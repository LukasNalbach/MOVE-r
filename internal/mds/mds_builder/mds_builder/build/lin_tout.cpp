/**
 * @brief builds L_in[0..p-1] and T_out[0..p-1] out of the disjoint interval sequence I
 */
template <typename uint_t>
void mds_builder<uint_t>::build_lin_tout() {
    L_in = std::vector<pair_list<uint_t>>(p);

    T_out = std::vector<pair_tree<uint_t>>(p,
        pair_tree<uint_t>(
            [](auto n1, auto n2){return n1.v.second < n2.v.second;},
            [](auto n1, auto n2){return n1.v.second > n2.v.second;},
            [](auto n1, auto n2){return n1.v.second == n2.v.second;}
        )
    );

    // [0..p] section starting positions in the range [0..n], 0 = s[0] < s[1] < ... < s[p-1] = n
    s = std::vector<uint_t>(p+1);
    s[0] = 0;
    s[p] = n;

    // [0..p], x[i] stores the number of input intervals in I starting before s[i]
    std::vector<uint_t> x(p+1);
    x[0] = 0;
    x[p] = k;

    // [0..p], u[i] stores the number of output intervals in I starting before s[i]
    std::vector<uint_t> u(p+1);
    u[0] = 0;
    u[p] = k;

    if (log) log_message("building pi");

    // create identity permutation pi of [0..k-1]
    pi = (uint_t*) malloc(k*sizeof(uint_t));
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k; i++) {
        pi[i] = i;
    }

    // sort pi by q
    auto comp = [this](uint_t i, uint_t j){return I[i].second < I[j].second;};
    if (p > 1) {
        ips4o::parallel::sort(&pi[0],&pi[k],comp);
    } else {
        ips4o::sort(&pi[0],&pi[k],comp);
    }

    if (log) {
        if (mf.is_open()) mf << " time_build_pi=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building L_in");
    }

    // calculate seperation positions
    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * \lfloor 2k/p \rfloor for s[i_p].
        uint_t o = i_p*((2*k)/p);

        // Left interval limit of the binary search for s[i_p].
        uint_t l_s;
        // Left interval limit of the binary search for x[i_p].
        uint_t l_x;
        // Left interval limit of the binary search for u[i_p].
        uint_t l_u;
        // Candidate position in the binary search for s[i_p].
        uint_t m_s;
        // Candidate position in the binary search for x[i_p].
        uint_t m_x;
        // Candidate position in the binary search for u[i_p].
        uint_t m_u;
        // Right interval limit of the binary search for s[i_p].
        uint_t r_s;
        // Right interval limit of the binary search for x[i_p].
        uint_t r_x;
        // Right interval limit of the binary search for u[i_p].
        uint_t r_u;

        // Initialize the search range for s[i_p] to [0..n-1].
        l_s = 0;
        r_s = n-1;

        // Perform a binary search over [0..n-1], to find s[i_p].
        do {
            /* Set the Candidate position for s[i_p] to the position
            in the middle between l_s and r_s. */
            m_s = l_s+(r_s-l_s)/2;

            // Find the minimum x' \in [0,k-1], s.t. p_{x'} >= m_s.
            l_x = 0;
            r_x = k-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (I[m_x].first < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' \in [0,k-1], s.t. q_{\pi[u']} >= m_s.
            l_u = 0;
            r_u = k-1;
            while (l_u != r_u) {
                m_u = l_u+(r_u-l_u)/2;
                if (I[pi[m_u]].second < m_s) {
                    l_u = m_u+1;
                } else {
                    r_u = m_u;
                }
            }

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x
            and l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x+l_u < o) {
                l_s = m_s + 1;
            } else {
                r_s = m_s;
            }
        } while (true);

        // Store l_s, l_x and l_u in s[i_p], x[i_p] and u[i_p], respectively.
        s[i_p] = l_s;
        x[i_p] = l_x;
        u[i_p] = l_u;
    }

    // write I into nodes and build L_in[0..p-1]
    nodes = (pair_tree_node<uint_t>*) malloc(k*sizeof(pair_tree_node<uint_t>));

    std::vector<std::queue<pair_list_node<uint_t>*>> Q_i(2*p);

    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i_p=0; i_p<2*p; i_p++) {
                uint_t l2 = x[i_p/2];
                uint_t r2 = x[i_p/2+1]-1;

                if (r2 >= l2) {
                    uint_t m2 = l2+(r2-l2)/2;

                    uint_t l = i_p%2 == 0 ? l2 : (m2+1);
                    uint_t r = i_p%2 == 0 ? m2 : r2;

                    #pragma omp task
                    {
                        for (uint_t i=l; i<=r; i++) {
                            nodes[i].v.v = I[i];
                            nodes[i].v.pr = &nodes[i-1].v;
                            nodes[i].v.sc = &nodes[i+1].v;

                            if (I[i+1].first - I[i].first > l_max) {
                                Q_i[i_p].emplace(&nodes[i].v);
                            }
                        }
                    }
                }
            }

            #pragma omp taskwait
        }
    }
    
    // Adjust the sizes of the lists in L_in, and set their head and tail nodes.
    for (uint16_t i=0; i<p; i++) {
        L_in[i].set_size(x[i+1]-x[i]);

        if (!L_in[i].empty()) {
            L_in[i].set_head(&nodes[x[i]].v);
            L_in[i].set_tail(&nodes[x[i+1]-1].v);
            L_in[i].head()->pr = NULL;
            L_in[i].tail()->sc = NULL;
        }
    }

    // Now, we do not need I anymore.
    free(I);
    I = NULL;

    if (log) {
        if (mf.is_open()) mf << " time_build_lin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    /* This function returns for the value i \in [0,k-1] the node in nodes[0..k-1],
    that stores the pair, which creates the i-th output interval. */
    std::function<pair_tree_node<uint_t>*(uint_t)> at = [this](uint_t i){
        return &nodes[pi[i]];
    };

    // build T_out[0..p-1] from nodes
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i=0; i<p; i++) {
                
                if (u[i+1] > u[i]) {
                    #pragma omp task
                    {
                        /* Build T_out[i] out of the pairs, that create the output starting
                        in the range[s[i]..s[i+1]-1], which are located at the positions
                        at[u[i]], at[u[i+1]], ..., at[u[i+1]-1]. */
                        T_out[i].insert_array(u[i],u[i+1]-1,at,2);
                    }
                }
            }

            #pragma omp taskwait
        }
    }

    if (log) {
        if (mf.is_open()) mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("splitting too long intervals");
    }

    free(pi);
    pi = NULL;

    // build new_nodes[0..p-1]
    new_nodes.reserve(p);

    for (uint16_t i=0; i<p; i++) {
        new_nodes.emplace_back(dg_io_nc<pair_tree_node<uint_t>>(k/(double)(16*p*(a-1))));
    }

    /* make sure each avl tree T_out[i], with i in [0..p-1], contains a pair creating
    an output interval starting at s[i] */
    for (uint16_t i=1; i<p; i++) {
        if (T_out[i].empty() || T_out[i].minimum()->v.v.second != s[i]) {
            pair_list_node<uint_t> *pln = &T_out[i-1].maximum()->v;

            pair_tree_node<uint_t> *ptn = new_nodes[i].emplace_back(
                pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{
                    pln->v.first+s[i]-pln->v.second,s[i]
                }))
            );

            // find i_ \in [0,p-1], so that s[i_] <= ptn->v.first < s[i_+1]
            uint16_t l = 0;
            uint16_t r = p-1;
            uint16_t m;
            while (l != r) {
                m = l+(r-l)/2+1;
                if (s[m] > pln->v.first) {
                    r = m-1;
                } else {
                    l = m;
                }
            }
            uint16_t i_ = l;

            L_in[i_].insert_after_node(&ptn->v,pln);
            T_out[i].insert_node(ptn);

            if (ptn->v.sc == NULL) {
                if (i_ < p-1) {
                    uint16_t i__ = i_+1;

                    while (i__ <= p-1 && L_in[i__].empty()) {
                        i__++;
                    }

                    if (i__ <= p-1 && !L_in[i__].empty() && L_in[i__].head()->v.first - ptn->v.v.first > l_max) {
                        Q_i[2*i_].emplace(&ptn->v);
                    }
                }
            } else if (ptn->v.sc->v.first - ptn->v.v.first > l_max) {
                Q_i[2*i_].emplace(&ptn->v);
            }
        }
    }

    /* make sure each list L_in[i], with i in [0..p-1], contains a pair creating an
    input interval starting at s[i] */
    for (uint16_t i=1; i<p; i++) {
        if (L_in[i].empty() || L_in[i].head()->v.first != s[i]) {
            pair_list_node<uint_t> *pln = L_in[i-1].tail();

            pair_tree_node<uint_t> *ptn = new_nodes[i].emplace_back(
                pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{
                    s[i],pln->v.second+s[i]-pln->v.first
                }))
            );

            // find i_ \in [0,p-1], so that s[i_] <= ptn->v.v.second < s[i_+1]
            uint16_t l = 0;
            uint16_t r = p-1;
            uint16_t m;
            while (l != r) {
                m = l+(r-l)/2+1;
                if (s[m] > ptn->v.v.second) {
                    r = m-1;
                } else {
                    l = m;
                }
            }
            uint16_t i_ = l;

            T_out[i_].insert_node(ptn);
            L_in[i].push_front_node(&ptn->v);

            if (ptn->v.sc == NULL) {
                if (i < p-1) {
                    uint16_t i__ = i+1;

                    while (i__ <= p-1 && L_in[i__].empty()) {
                        i__++;
                    }

                    if (i__ <= p-1 && !L_in[i__].empty() && L_in[i__].head()->v.first - ptn->v.v.first > l_max) {
                        Q_i[2*i].emplace(&ptn->v);
                    }
                }
            } else if (ptn->v.sc->v.first - ptn->v.v.first > l_max) {
                Q_i[2*i].emplace(&ptn->v);
            }
        }
    }

    // insert the pair (s[i+1],s[i+1]) into each avl tree T_out[i], with i in [0..p-1]
    for (uint16_t i=0; i<p; i++) {
        pair_tree_node<uint_t>* ptn = new_nodes[i].emplace_back(
            pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{
                s[i+1],s[i+1]
            }))
        );
        
        L_in[i].push_back_node(&ptn->v);
        T_out[i].insert_node(ptn);
    }

    std::vector<std::vector<std::queue<pair_tree_node<uint_t>*>>> Q_o(p,std::vector<std::queue<pair_tree_node<uint_t>*>>(p));

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        pair_list_node<uint_t> *pln_cur,*pln_last;
        pair_tree_node<uint_t> *ptn_tmp;

        for (uint16_t i=2*i_p; i<2*(i_p+1); i++) {
            while (!Q_i[i].empty()) {
                pln_cur = Q_i[i].front();
                Q_i[i].pop();
                
                if (pln_cur->sc != NULL && pln_cur->sc->v.first - pln_cur->v.first > l_max) {
                    // find i_p' \in [0,p-1], so that s[i_p'] <= ptn->v.v.second < s[i_p'+1]
                    uint16_t l = 0;
                    uint16_t r = p-1;
                    uint16_t m;
                    while (l != r) {
                        m = l+(r-l)/2+1;
                        if (s[m] > pln_cur->v.second) {
                            r = m-1;
                        } else {
                            l = m;
                        }
                    }
                    uint16_t i_p_ = l;

                    do {
                        pln_last = pln_cur;

                        ptn_tmp = new_nodes[i_p].emplace_back(
                            pair_tree_node<uint_t>(pair_list_node<uint_t>(std::pair<uint_t,uint_t>{
                                pln_last->v.first+l_max,
                                pln_last->v.second+l_max
                            }))
                        );

                        pln_cur = &ptn_tmp->v;
                        L_in[i_p].insert_after_node(pln_cur,pln_last);

                        while (pln_cur->v.second >= s[i_p_+1]) {
                            i_p_++;
                        }
                        
                        Q_o[i_p_][i_p].emplace(ptn_tmp);
                    } while (pln_cur->sc != NULL && pln_cur->sc->v.first - pln_cur->v.first > l_max);
                }
            }
        }
    }

    for (uint16_t i=0; i<p; i++) {
        L_in[i].remove_node(L_in[i].tail());
    }

    Q_i.clear();

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        for (uint16_t i=0; i<p; i++) {
            while (!Q_o[i_p][i].empty()) {
                T_out[i_p].insert_node(Q_o[i_p][i].front());
                Q_o[i_p][i].pop();
            }
        }
    }

    Q_o.clear();

    if (log) {
        if (mf.is_open()) mf << " time_split_too_long_input_intervals=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}