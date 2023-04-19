/**
 * @brief balances the output interval [q_j, q_j + d_j) and all unbalanced output intervals
 *        starting before or at q_u that have become unbalanced in the process
 * @param pln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
 *                to [q_j, q_j + d_j) in the permutation graph
 * @param ptn_J (p_j,q_j), [q_j, q_j + d_j) must be the first unbalanced output interval
 * @param q_u starting position of an output interval
 * @param p_cur starting position of the current input interval
 * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
 * @return the newly created pair (p_j+d,q_j+d)
 */
template <typename uint_t>
inline pair_tree_node<uint_t>* mds_builder<uint_t>::balance_upto_v3_seq(
    pair_list_node<uint_t> *pln_IpA,
    pair_tree_node<uint_t> *ptn_J,
    uint_t q_u,
    uint_t p_cur,
    uint_t *i_
) {
    uint_t p_j = ptn_J->v.v.first;
    uint_t q_j = ptn_J->v.v.second;
    uint_t d_j = interval_length_seq(&ptn_J->v);
    
    // d = p_{i+2a} - q_j is the maximum integer, so that [q_j, q_j + d) has a incoming edges in the permutation graph.
    uint_t d = pln_IpA->v.first - q_j;

    // Create the pair (p_j + d, q_j + d), which creates two new input intervals [p_j, p_j + d) and [p_j + d, p_j + d_j).
    pair_tree_node<uint_t> *ptn_NEW = new_nodes[0].emplace_back(pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{p_j + d, q_j + d})));
    T_out[0].insert_node_in(ptn_NEW,ptn_J);
    L_in[0].insert_after_node(&ptn_NEW->v,&ptn_J->v);

    /*
        Case 1: pj + d ≥ qu + du
            It is irrelevant, if [qy, qy + qy) is unbalanced, since it starts after qu.
        Case 2: pj + d ∈ [qu, qu + du)
            Then [qy, qy +dy) = [qu, qu+du), so there are a+1 input intervals starting in [qu, qu+du),
            hence, it and all output intervals starting before qu are balanced.
        Case 3: pj + d < qu
            Case 3.1: pj + d ∈ [qj , qj + dj)
                Then [qy, qy + dy) equals either [qj , qj + d) or [qj + d, qj + dj). Since before inserting
                (pj + d, qj + d), there were < 2a input intervals starting in [qj , qj + dj), there are now
                ≤ a + 1 input intervals starting in [qj , qj + d) and [qj + d, qj + dj) each, hence, they
                are both balanced.
            Case 3.2: pj + d /∈ [qj , qj + dj)
                [qy, qy + dy) is possibly unbalanced.
    */

    // check if case 3.2 holds
    if (p_j + d < q_u) {
        if (p_j + d < q_j || q_j + d_j <= p_j + d) {
            // find [q_y, q_y + d_y)
            pair_tree_node<uint_t> *ptn_Y = T_out[0].maximum_leq(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{0,p_j + d}));

            // find the first input interval [p_z, p_z + d_z), that is connected to [q_y, q_y + d_y) in the permutation graph
            pair_list_node<uint_t> *pln_Z = &ptn_NEW->v;
            uint_t i__ = 1;
            while (pln_Z->pr != NULL && pln_Z->pr->v.first >= ptn_Y->v.v.second) {
                pln_Z = pln_Z->pr;
                i__++;
            }
            pln_Z = &ptn_NEW->v;

            pair_list_node<uint_t> *pln_ZpA = is_unbalanced(&pln_Z,&i__,ptn_Y);
            if (pln_ZpA != NULL) {
                balance_upto_v3_seq(pln_ZpA,ptn_Y,q_u,p_cur,i_);
            }
        }
    } else if (p_j + d < p_cur) {
        (*i_)++;
    }

    return ptn_NEW;
}

/**
 * @brief balances the disjoint interval sequence in L_in[0] and T_out[0] sequentially
 */
template <typename uint_t>
void mds_builder<uint_t>::balance_v3_seq() {
    if (log) log_message("balancing");

    // points to to the pair (p_i,q_i).
    pair_list_node<uint_t> *pln_I = L_in[0].head();
    // points to the pair (p_j,q_j).
    typename pair_tree<uint_t>::avl_it it_outp_cur = T_out[0].iterator();
    // points to the pair (p_{j'},q_{j'}), where q_j + d_j = q_{j'}.
    typename pair_tree<uint_t>::avl_it it_outp_nxt = T_out[0].iterator(T_out[0].second_smallest());

    // temporary variables
    pair_list_node<uint_t> *pln_IpA;
    uint_t i_ = 1;

    // At the start of each iteration, [p_i, p_i + d_i) is the first input interval connected to [q_j, q_j + d_j) in the permutation graph
    bool stop = false;
    do {
        pln_IpA = is_unbalanced(&pln_I,&i_,it_outp_cur.current(),it_outp_nxt.current());

        // If [q_j, q_j + d_j) is unbalanced, balance it and all output intervals starting before it, that might get unbalanced in the process.
        if (pln_IpA != NULL) {
            it_outp_cur.set(balance_upto_v3_seq(pln_IpA,it_outp_cur.current(),it_outp_cur.current()->v.v.second,pln_I->v.first,&i_));
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

    if (log) {
        if (mf.is_open()) {
            mf << " time_balance_phase_1=" << time_diff_ns(time)
                << " time_balance_phase_2=" << 0;
        }
        time = log_runtime(time);
    }
}