/**
 * @brief balances the output interval [q_j, q_j + d_j) by inserting the newly created pair into
 * T_out[i_p] and Q[0..p-1][i_p]
 * @param Q reference to Q
 * @param pln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
 *                to [q_j, q_j + d_j) in the permutation graph
 * @param ptn_J (p_j,q_j), [q_j, q_j + d_j) must be the first unbalanced output interval starting
 *              in [s[i_p]..q_u]
 * @param ptn_J_next (p_j',q_j'), [q_j', q_j' + d_j') must be the first output interval starting
 *                   after [q_j, q_j + d_j)
 * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]]
 * @param p_cur starting position of the current input interval
 * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
 * @return the newly created pair (p_j+d,q_j+d)
 */
template <typename uint_t>
inline pair_tree_node<uint_t>* mds_builder<uint_t>::balance_upto_v4_par(
    qt_v4<uint_t> &Q,
    pair_list_node<uint_t> *pln_IpA,
    pair_tree_node<uint_t> *ptn_J,
    pair_tree_node<uint_t>* ptn_J_nxt,
    uint_t q_u,
    uint_t p_cur,
    uint_t *i_
) {
    // Index \in [0..p-1] of the current thread.
    uint16_t i_p = omp_get_thread_num();

    uint_t p_j = ptn_J->v.v.first;
    uint_t q_j = ptn_J->v.v.second;
    uint_t d_j = ptn_J_nxt->v.v.second - q_j;

    // d = p_{i+2a} - q_j is the maximum integer, so that [q_j, q_j + d) has a incoming edges in the permutation graph.
    uint_t d = pln_IpA->v.first - q_j;

    // Create the pair (p_j + d, q_j + d), which creates two new input intervals [p_j, p_j + d) and [p_j + d, p_j + d_j).
    pair_tree_node<uint_t> *ptn_NEW = new_nodes[i_p].emplace_back(pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{p_j + d, q_j + d})));
    T_out[i_p].insert_node_in(ptn_NEW,ptn_J);

    if (!(s[i_p] <= p_j + d && p_j + d < s[i_p+1])) {
        // If the new pair must be inserted in L_in[i_p_] of another thread i_p_ != i_p, find i_p_ with a binary search.
        uint16_t l = 0;
        uint16_t r = p-1;
        uint16_t m;
        while (l != r) {
            m = l+(r-l)/2+1;
            if (p_j + d >= s[m]) {
                l = m;
            } else {
                r = m-1;
            }
        }
        uint16_t i_p_ = l;

        Q[i_p_][i_p].enqueue(q_pair<uint_t>{&ptn_NEW->v,&ptn_J->v});
    } else {
        // Else insert it in L_in[i_p].
        L_in[i_p].insert_after_node(&ptn_NEW->v,&ptn_J->v);

        if (p_j + d < q_u) {
            if (p_j + d < q_j || q_j + d_j <= p_j + d) {
                pair_tree_node<uint_t> *ptn_Y = T_out[i_p].maximum_leq(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{0,p_j + d}));
                uint_t q_y = ptn_Y->v.v.second;

                // find the output interval starting after [q_y, q_y + d_y)
                pair_tree_node<uint_t> *ptn_Y_nxt = ptn_Y->nxt();

                // find the first input interval [p_z, p_z + d_z), that is connected to [q_y, q_y + d_y) in the permutation graph
                pair_list_node<uint_t> *pln_Z = &ptn_NEW->v;
                uint_t i__ = 1;
                while (pln_Z->pr != NULL && pln_Z->pr->v.first >= q_y) {
                    pln_Z = pln_Z->pr;
                    i__++;
                }
                pln_Z = &ptn_NEW->v;

                pair_list_node<uint_t> *pln_ZpA = is_unbalanced(&pln_Z,&i__,ptn_Y,ptn_Y_nxt);
                if (pln_ZpA != NULL) {
                    balance_upto_v4_par(Q,pln_ZpA,ptn_Y,ptn_Y_nxt,q_u,p_cur,i_);
                }
            }
        } else if (p_j + d < p_cur) {
            (*i_)++;
        }
    }
    return ptn_NEW;
}

/**
 * @brief balances the disjoint interval sequence in L_in[0..p-1] and T_out[0..p-1] in parallel
 */
template <typename uint_t>
void mds_builder<uint_t>::balance_v4_par() {
    if (log) log_message("balancing phase 1");
    
    /** @brief [0..p-1] stores queues with tuples (*p1,*p2);
     *        Q[i] stores the tuples to insert into thread i's section [s[i]..s[i+1] */
    qt_v4<uint_t> Q(p);

    /** @brief [0..p-1] D_done[i] stores whether thread i is done */
    bool D_done[p];

    for (uint16_t i=0; i<p; i++) {
        Q[i] = std::vector<moodycamel::ConcurrentQueue<q_pair<uint_t>>>(p);
        D_done[i] = false;
    }

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // points to to the pair (p_i,q_i).
        pair_list_node<uint_t> *pln_I = L_in[i_p].head();
        // points to the pair (p_j,q_j).
        typename pair_tree<uint_t>::avl_it it_outp_cur = T_out[i_p].iterator();
        // points to the pair (p_{j'},q_{j'}), where q_j + d_j = q_{j'}.
        typename pair_tree<uint_t>::avl_it it_outp_nxt = T_out[i_p].iterator(T_out[i_p].second_smallest());

        // temporary variables
        pair_list_node<uint_t> *pln_IpA;
        uint_t i_ = 1;

        // At the start of each iteration, [p_i, p_i + d_i) is the first input interval connected to [q_j, q_j + d_j) in the permutation graph
        bool stop = false;
        do {
            pln_IpA = is_unbalanced(&pln_I,&i_,it_outp_cur.current(),it_outp_nxt.current());

            // If [q_j, q_j + d_j) is unbalanced, balance it and all output intervals starting before it, that might get unbalanced in the process.
            if (pln_IpA != NULL) {
                it_outp_cur.set(balance_upto_v4_par(Q,pln_IpA,it_outp_cur.current(),it_outp_nxt.current(),it_outp_cur.current()->v.v.second,pln_I->v.first,&i_));
                continue;
            }

            // Find the next output interval with an incoming edge in the permutation graph and the first input interval connected to it.
            do {
                if (!it_outp_nxt.has_next()) {stop = true; break;}
                it_outp_cur.set(it_outp_nxt.current());
                it_outp_nxt.next();
                while (pln_I->v.first < it_outp_cur.current()->v.v.second) {
                    if (pln_I->sc == NULL) {stop = true; break;}
                    pln_I = pln_I->sc;
                }
            } while (!stop && pln_I->v.first >= it_outp_nxt.current()->v.v.second);
            i_ = 1;
        } while (!stop);
    }

    if (log) {
        if (mf.is_open()) mf << " time_balance_phase_1=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("balancing phase 2");
    }

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // temporary variables
        q_pair<uint_t> pair_ins;
        pair_list_node<uint_t> *pln_IpA,*pln_I,*pln_Im1,*pln_Z,*pln_ZpA;
        pair_tree_node<uint_t> *ptn_Y,*ptn_Y_nxt;
        bool done_this,done_other;
        uint_t i_ = 1;
        uint_t q_y;

        done_other = false;
        do {
            done_this = true;

            for (uint16_t i=0; i<p; i++) {
                while (Q[i_p][i].try_dequeue(pair_ins)) {
                    if (done_this) {
                        done_this = false;
                        if (D_done[i_p]) {
                            D_done[i_p] = false;
                        }
                    }

                    pln_I = pair_ins.first;
                    pln_Im1 = pair_ins.second;

                    L_in[i_p].insert_after_node(pln_I,pln_Im1);

                    // check if an output interval could have become unbalanced by inserting the new pair
                    ptn_Y = T_out[i_p].maximum_leq(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{0,pln_I->v.first}));
                    q_y = ptn_Y->v.v.second;

                    // find the output interval starting after [q_y, q_y + d_y)
                    ptn_Y_nxt = ptn_Y->nxt();

                    // find the first input interval [p_z, p_z + d_z), that is connected to [q_y, q_y + d_y) in the permutation graph
                    pln_Z = pln_I;
                    i_ = 1;
                    while (pln_Z->pr != NULL && pln_Z->pr->v.first >= q_y) {
                        pln_Z = pln_Z->pr;
                        i_++;
                    }
                    pln_Z = pln_I;

                    pln_ZpA = is_unbalanced(&pln_Z,&i_,ptn_Y,ptn_Y_nxt);
                    if (pln_ZpA != NULL) {
                        balance_upto_v4_par(Q,pln_ZpA,ptn_Y,ptn_Y_nxt,s[i_p+1],s[i_p+1],&i_);
                    }
                }
            }

            if (done_this) {
                done_other = true;
                for (uint16_t i=0; i<p; i++) {
                    if (!D_done[i] || Q[i][i_p].size_approx() != 0) {
                        done_other = false;
                    }
                }

                if (!D_done[i_p]) {
                    D_done[i_p] = true;
                }
            }
        } while (!done_other || !done_this);
    }

    if (log) {
        if (mf.is_open()) mf << " time_balance_phase_2=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}