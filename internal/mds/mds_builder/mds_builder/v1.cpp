/**
 * @brief builds the move datastructure md
 */
template <typename uint_t>
void mds_builder<uint_t>::balance_v1() {
    size_t baseline;
    std::chrono::steady_clock::time_point time, time_start;

    if (log) {
        baseline = malloc_count_current() - sizeof(I[0])*k;
        time = now();
        time_start = time;
        std::cout << std::endl;
    }

    if (log) log_message("building T_in");

    // stores the pairs in I sorted by p_i
    avl_tree<std::pair<uint_t,uint_t>> T_in(
        [](auto n1, auto n2){return n1.first < n2.first;},
        [](auto n1, auto n2){return n1.first > n2.first;},
        [](auto n1, auto n2){return n1.first == n2.first;}
    );

    // stores the pairs in I sorted by q_i
    avl_tree<std::pair<uint_t,uint_t>> T_out(
        [](auto n1, auto n2){return n1.second < n2.second;},
        [](auto n1, auto n2){return n1.second > n2.second;},
        [](auto n1, auto n2){return n1.second == n2.second;}
    );

    /* stores the pairs in I sorted by p_i, whiches output intervals have at least 4 incoming edges in
    the permutation graph */
    avl_tree<std::pair<uint_t,uint_t>> T_e(
        [](auto n1, auto n2){return n1.first < n2.first;},
        [](auto n1, auto n2){return n1.first > n2.first;},
        [](auto n1, auto n2){return n1.first == n2.first;}
    );

    // build T_in and T_out
    for (uint_t i=0; i<k; i++) {
        T_in.insert_or_update(I[i]);
    }

    T_in.insert_or_update(std::make_pair(n,n));

    if (log) {
        if (mf.is_open()) mf << " time_build_tin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    for (uint_t i=0; i<k; i++) {
        T_out.insert_or_update(I[i]);
    }

    T_out.insert_or_update(std::make_pair(n,n));

    avl_node<std::pair<uint_t,uint_t>> *node_cur = T_in.minimum();

    while (node_cur != T_in.maximum()) {
        while (node_cur->nxt()->v.first - node_cur->v.first > l_max) {
            T_out.insert_or_update(std::make_pair(node_cur->v.first+l_max,node_cur->v.second+l_max));
            node_cur = T_in.insert_or_update(std::make_pair(node_cur->v.first+l_max,node_cur->v.second+l_max));
        }
        node_cur = node_cur->nxt();
    }

    if (log) {
        if (mf.is_open()) mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_e");
    }

    // Now, we do not need I anymore.
    free(I);
    I = NULL;

    // build T_e
    uint_t q_i,q_next;
    node_cur = T_in.minimum();
    avl_node<std::pair<uint_t,uint_t>> *node_cur_2;
    uint32_t e;
    while (node_cur != T_in.maximum()) {
        /* For each output interval [q_i, q_i + d_i), find the first input interval connected
        to it in the permutation graph. */
        q_i = node_cur->v.second;
        q_next = q_i + node_cur->nxt()->v.first - node_cur->v.first;
        node_cur_2 = T_in.minimum_geq(std::pair<uint_t,uint_t>{q_i,0});
        // Count the number of input intervals connected to it in the permutation graph.
        e = 0;
        while (node_cur_2 != NULL) {
            if (node_cur_2->v.first < q_next) {
                e++;
                if (e == 2*a) {
                    // If there are at least 2a, insert it's corresponding pair into T_e.
                    T_e.insert_or_update(node_cur->v);
                    break;
                }
            } else {
                break;
            }
            node_cur_2 = node_cur_2->nxt();
        }
        node_cur = node_cur->nxt();
    }

    if (log) {
        if (mf.is_open()) mf << " time_build_te=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("balancing");
    }

    // balance the disjoint interval sequence
    uint_t d,q_j,p_j,q_y,d_j,d_y;
    std::pair<uint_t,uint_t> pair_NEW,pair_Y;
    avl_node<std::pair<uint_t,uint_t>> *node_Ipa,*min,*node_NEW,*node_Y;
    std::vector<std::tuple<uint_t,uint_t,std::pair<uint_t,uint_t>>> intervals_to_check;
    while (!T_e.empty()) {
        /* Find the pair creating the first unbalanced output interval [q_j, q_j + d_j)
        and remove it from T_e. */
        min = T_e.minimum();
        p_j = min->v.first;
        q_j = min->v.second;
        delete T_e.remove(min->v);

        // Find the a+1-st input interval in [q_j, q_j + d_j) and set d = p_{i+a} - q_j.
        /* d is the smallest integer, so that [q_j, q_j + d) has a incoming edges in the
        permutation graph. */
        node_Ipa = T_in.minimum_geq(std::make_pair(q_j,0));
        for (uint16_t i=0; i<a; i++) {
            node_Ipa = node_Ipa->nxt();
        }
        d = node_Ipa->v.first-q_j;

        // Create the new pair (p_j + d, q_j + d) and insert it into T_in and T_out.
        pair_NEW = std::make_pair(p_j+d,q_j+d);
        T_in.insert_or_update(pair_NEW);
        node_NEW = T_out.insert_or_update(pair_NEW);

        /* Find the output interval [q_y, q_y + d_y), [p_j + d, p_j + d_j) is connected
        to in the permutation graph. */
        node_Y = T_out.maximum_leq(std::make_pair(0,p_j+d));
        pair_Y = node_Y->v;
        q_y = pair_Y.second;

        /* The number of input intervals connected to [q_j + d, q_j + d_j) and [q_y, q_y + d_y)
        may have changed. For each, check if it has at least 2a incoming edges in the permutation
        graph and insert it into T_e, if it has. */
        d_j = node_NEW->nxt()->v.second-q_j;
        d_y = node_Y->nxt()->v.second-q_y;
        intervals_to_check = {
            std::make_tuple(q_j+d,q_j+d_j-1,pair_NEW),
            std::make_tuple(q_y,q_y+d_y-1,pair_Y)
        };

        for (auto tup : intervals_to_check) {
            /* Find the first input interval connected to the output interval in the
            permutation graph. */
            node_cur = T_in.minimum_geq(std::make_pair(std::get<0>(tup),0));

            // Count the number of input intervals connected to it in the permutation graph.
            e = 0;
            while (node_cur != NULL) {
                if (node_cur->v.first <= std::get<1>(tup)) {
                    e++;
                    // If there are at least 2a, insert it's corresponding pair into T_e.
                    if (e == 2*a) {
                        T_e.insert_or_update(std::get<2>(tup));
                        break;
                    }
                } else {
                    break;
                }
                node_cur = node_cur->nxt();
            }
        }
    }
    // Because T_e is empty, there are no unbalanced output intervals.

    k_ = T_in.size()-1;

    if (log) {
        time = log_runtime(time);
        float k__k = std::round(100.0*k_/k)/100.0;
        if (mf.is_open()) {
            mf << " k=" << k;
            mf << " k_=" << k_;
            mf << " time_balance=" << time_diff_ns(time);
        }
        std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
        log_message("writing B(I)");
    }
    
    D_p = (uint_t*) malloc((k_+1)*sizeof(uint_t));
    D_q = (uint_t*) malloc((k_+1)*sizeof(uint_t));

    D_p[k_] = n;
    D_q[k_] = n;

    auto it = T_in.iterator();

    for (uint_t i=0; i<=k_; i++) {
        D_p[i] = it.current()->v.first;
        D_q[i] = it.current()->v.second;
        it.next();
    }

    // Deconstruct the additional data structures.    
    T_in.delete_nodes();
    T_out.delete_nodes();

    if (log) {
        if (mf.is_open()) mf << " time_write_b_i=" << time_diff_ns(time);
        time = log_runtime(time);
    }

    // Build D_offs and D_idx
    build_doffsidx();

    if (log) {
        log_message("move datastructure built");
        if (mf.is_open()) {
            std::cout << ", peak additional memory allocation during build: ~ "
                    << format_size(malloc_count_peak()-baseline);
                    
            mf << " time_total=" << time_diff_ns(time_start);
        }
    }
}