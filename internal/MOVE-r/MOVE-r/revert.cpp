/**
 * @brief writes the original text to the string T
 * @param p the number of threads to use, must not exceed the maximum
 * number ofthreads to use that has been specified when the index was built
 * @return T
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
uint8_t* MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::revert(uint16_t p) {
    uint8_t* T = (uint8_t*) malloc(n-1);

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        uint16_t b_pr = (i_p*p_r)/p;
        uint16_t e_pr = i_p == p-1 ? p_r-1 : ((i_p+1)*p_r)/p-1;

        // Iteration range start position of thread i_p.
        int64_t b = b_pr == 0 ? 0 : b_pr*((n-1)/p_r)+1;
        // Iteration range end position of thread i_p.
        int64_t e = e_pr == p_r-1 ? n-2 : (e_pr+1)*((n-1)/p_r);

        // The position in the bwt of the current character in T.
        pos_t i = D_e[e_pr];
        /* The position of the character in L' (and therefore the index of the input interval
        in M_LF), in which i lies. */
        idx_t x;

        // Calculate x with a binary search over the input intervals of M_LF.

        // Left interval limit of the binary search.
        idx_t l = 0;
        // Left interval limit of the binary search.
        idx_t r = r_-1;
        // Candidate position for the binary search.
        idx_t m;

        while (l != r) {
            m = l+(r-l)/2+1;
            if (M_LF.p(m) <= i) {
                l = m;
            } else {
                r = m-1;
            }
        }

        x = l;

        if (chars_mapped) {
            // Set T[e] <- L[i] = L'[x]
            T[e] = map_char_rev[L_(x)];

            // Build T from right to left starting at the second to the last psition.
            for (int64_t j=e-1; j>=b; j--) {
                // Set i <- LF(i).
                M_LF.move(i,x);
                // Set T[j] <- L[i] = L'[x].
                T[j] = map_char_rev[L_(x)];
            }
        } else {
            T[e] = L_(x);

            for (int64_t j=e-1; j>=b; j--) {
                M_LF.move(i,x);
                T[j] = L_(x);
            }
        }
    }

    return T;
}