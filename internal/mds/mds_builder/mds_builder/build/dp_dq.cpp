/**
 * @brief inserts the pairs in L_in[0..p-1] into D_pair and writes D_idx
 */
template <typename uint_t>
void mds_builder<uint_t>::build_dp_dq() {
    // [0..p-1], x[i] stores the number of input intervals in I starting before s[i]
    std::vector<uint_t> x(p+1);
    x[0] = 0;

    for (uint16_t i=0; i<p; i++) {
        // remove the pairs (s[i+1],s[i+1]) from T_out[i], for each i \in [0..p-1]
        T_out[i].remove_node(T_out[i].maximum());

        // calculate x and u
        x[i+1] = x[i]+L_in[i].size();
    }

    k_ = x[p];
    s.clear();

    if (log) {
        if (log) {
            float k__k = std::round(100.0*k_/k)/100.0;
            if (mf.is_open()) {
                mf << " k=" << k;
                mf << " k_=" << k_;
            }
            std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
        }
        log_message("writing B_I");
    }

    D_p = (uint_t*) malloc((k_+1)*sizeof(uint_t));
    D_q = (uint_t*) malloc((k_+1)*sizeof(uint_t));

    D_p[k_] = n;
    D_q[k_] = n;

    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i_p=0; i_p<p; i_p++) {
                uint_t l = x[i_p];
                uint_t r = x[i_p+1]-1;
                uint_t m = l+(r-l)/2;

                #pragma omp task
                {
                    auto pln_cur = L_in[i_p].head();

                    for (uint_t i=l; i<=m; i++) {
                        D_p[i] = pln_cur->v.first;
                        D_q[i] = pln_cur->v.second;
                        pln_cur = pln_cur->sc;
                    }
                }

                #pragma omp task
                {
                    auto pln_cur = L_in[i_p].tail();

                    for (uint_t i=r; i>m; i--) {
                        D_p[i] = pln_cur->v.first;
                        D_q[i] = pln_cur->v.second;
                        pln_cur = pln_cur->pr;
                    }
                }
            }

            #pragma omp taskwait
        }
    }

    // Deconstruct the additional data structures.
    x.clear();

    for (uint16_t i=0; i<p; i++) {
        L_in[i].disconnect_nodes();
        T_out[i].disconnect_nodes();
        new_nodes[i].clear();
    }
    
    L_in.clear();
    T_out.clear();
    new_nodes.clear();

    free(nodes);
    nodes = NULL;

    if (log) {
        if (mf.is_open()) mf << " time_write_b_i=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}