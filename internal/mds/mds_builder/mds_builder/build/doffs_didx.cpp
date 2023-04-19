/**
 * @brief builds D_offs and D_idx
 */
template <typename uint_t>
void mds_builder<uint_t>::build_doffsidx() {
    if (log)  log_message("building D_offs and D_idx");

    pi = (uint_t*) malloc((k_+1)*sizeof(uint_t));
    pi[k_] = k_;

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        pi[i] = i;
    }

    auto comp_pi = [this](auto i, auto j) {return D_q[i] < D_q[j];};
    if (p > 1) {
        ips4o::parallel::sort(&pi[0],&pi[k_],comp_pi);
    } else {
        ips4o::sort(&pi[0],&pi[k_],comp_pi);
    }

    D_offsidx = (std::pair<uint_t,uint_t>*) malloc(k_*sizeof(std::pair<uint_t,uint_t>));

    std::vector<uint_t> x(p+1);
    x[0] = 0;
    x[p] = k_;
    
    std::vector<uint_t> u(p+1);
    u[0] = 0;
    u[p] = k_;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * \lfloor (r'+r'')/p \rfloor for s[i_p].
        uint_t o = i_p*((2*k_)/p);

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
            /* Set the Candidate position for s[i_p] to the position in the middle
            between l_s and r_s. */
            m_s = l_s+(r_s-l_s)/2;

            // Find the minimum x' \in [0,r''-1], s.t. M_Phi.p(x') >= m_s.
            l_x = 0;
            r_x = k_-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (D_p[m_x] < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' \in [0,r'-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = 0;
            r_u = k_-1;
            while (l_u != r_u) {
                m_u = l_u+(r_u-l_u)/2;
                if (D_q[pi[m_u]] < m_s) {
                    l_u = m_u+1;
                } else {
                    r_u = m_u;
                }
            }

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x and
            l_u are valid values for x[i_p] and u[i_p], respectively. */
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

        // Store l_x and l_u in x[i_p] and u[i_p], respectively.
        x[i_p] = l_x;
        u[i_p] = l_u;
        
        #pragma omp barrier

        /* Check if the range [u[i_p]..u[i_p+1]-1], over which the thread i_p
        has to iterate in SA_s, is empty. */
        if (u[i_p+1] > u[i_p]) {
            // Iteration range start position in M_Phi.
            uint_t i = x[i_p];
            // Iteration range start position in SA_s.
            uint_t j = u[i_p];

            // Iteration range end position in M_Phi.
            uint_t i_ = x[i_p+1]-1;
            // Iteration range end position in SA_s.
            uint_t j_ = u[i_p+1]-1;

            /* Check if the first suffix array sample lies before x[i_p]-th 
            input interval of M_Phi. */
            if (D_q[pi[j]] < D_p[i]) {
                i--;
            }

            /* Iterate until one of the iteration end positions i_ and j_ has
            been reached. */
            while (i <= i_ && j <= j_) {

                /* Iterate over the suffix array samples that lie in the current
                i-th input interval of M_Phi. */
                while (j <= j_ && D_q[pi[j]] < D_p[i+1]) {
                    
                    /* Because each of those j-th largest suffix array samples lie in the i-th
                    input interval of M_Phi, we can set SA_idx[\pi'[j]] = i for each of them. */
                    D_offsidx[pi[j]].first = D_q[pi[j]]-D_p[i];
                    D_offsidx[pi[j]].second = i;

                    j++;
                }

                i++;
            };
        }
    }

    if (log) {
        if (mf.is_open()) mf << " time_build_doffs_didx=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}