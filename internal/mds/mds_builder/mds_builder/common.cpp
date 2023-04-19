/**
 * @brief returns the length of the input/output interval starting at p_i/q_i
 * @param pln (p_i,q_i)
 * @return length of the input/output interval starting at p_i/q_i
 */
template <typename uint_t>
uint_t mds_builder<uint_t>::interval_length_seq(pair_list_node<uint_t> *pln) {
    return (pln->sc != NULL ? pln->sc->v.first : n) - pln->v.first;
}

/**
 * @brief checks if the output interval [q_j, q_j + d_j) is unbalanced and iterates
 *        pln_IpI_ tothe last output interval connected to it in the permutation graph
 * @param pln_IpI_ (p_{i+i_},q_{i+i_}), [p_i, p_i + d_i) must be the first input
 *                 interval connected to [q_j, q_j + d_j) in the permutation graph
 * @param i_ 1 <= i_ <= 2a
 * @param ptn_J (p_j,q_j)
 * @param ptn_J_nxt (p_{j'},q_{j'}), with q_j + d_j = q_{j'}
 * @return (p_{i+a},q_{i+a}) if [q_j, q_j + d_j) is unbalanced, else NULL
 */
template <typename uint_t>
pair_list_node<uint_t>* mds_builder<uint_t>::is_unbalanced(
    pair_list_node<uint_t> **pln_IpI_,
    uint_t* i_,
    pair_tree_node<uint_t> *ptn_J,
    pair_tree_node<uint_t> *ptn_J_nxt
) {
    // [l,r] = [q_j, q_j + d_j)

    // r + 1
    uint_t rp1 = ptn_J_nxt == NULL ?
        ptn_J->v.v.second + interval_length_seq(&ptn_J->v)
        : ptn_J_nxt->v.v.second;

    /* If |[l,r]| < 2a, there cannot be at least 2a input intervals
    connected to [l,r] in the permutation graph. */
    if (rp1-ptn_J->v.v.second < 2*a) return NULL;

    uint_t i_start = *i_;

    /* Count the number i of input intervals connected to [l,r]
    in the permutation graph and stop as soon as i > a. */
    do {
        if (*i_ > a) break;
        if ((*pln_IpI_)->sc == NULL) return NULL;
        *pln_IpI_ = (*pln_IpI_)->sc;
        (*i_)++;
        if ((*pln_IpI_)->v.first >= rp1) return NULL;
    } while (true);

    pair_list_node<uint_t> *pln_IpA = *pln_IpI_;

    // Count further and stop as soon as i >= 2a.
    do {
        if (*i_ >= 2*a) {
            // if i_ > a+1 held in the beginning, correct pln_IpA
            while (i_start > a+1) {
                pln_IpA = pln_IpA->pr;
                i_start--;
            }
            *i_ = a;
            return pln_IpA;
        }
        if ((*pln_IpI_)->sc == NULL) return NULL;
        *pln_IpI_ = (*pln_IpI_)->sc;
        (*i_)++;
        if ((*pln_IpI_)->v.first >= rp1) return NULL;
    } while (true);
}

/**
 * @brief verifies the correctness of the resulting interval sequence
 */
template <typename uint_t>
void mds_builder<uint_t>::verify_correctness() {
    std::cout << std::endl << "verifying correctness of the interval sequence:" << std::endl;
    bool correct = true;

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (D_p[i+1] <= D_p[i]) {
            #pragma omp critical
            {
                std::cout << "input interval starting positions are not ascending:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << D_p[i] << std::endl;
                std::cout << "p_{i+1} = " << D_p[i+1] << std::endl << std::endl;
                correct = false;
            }
        }
    }

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (D_p[i+1] - D_p[i] > l_max) {
            #pragma omp critical
            {
                std::cout << "too long input interval (> l_max = " << l_max << "):" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << D_p[i] << std::endl;
                std::cout << "p_{i+1} = " << D_p[i+1] << std::endl << std::endl;
                correct = false;
            }
        }
    }
    
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k_; i++) {
        if (D_q[pi[i+1]] - D_q[pi[i]] != D_p[pi[i]+1] - D_p[pi[i]]) {
            #pragma omp critical
            {
                std::cout << "input- and output interval lengths do not match:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_{pi[i]} = " << D_p[pi[i]] << std::endl;
                std::cout << "p_{pi[i]+1} = " << D_p[pi[i]+1] << std::endl;
                std::cout << "q_{pi[i]} = " << D_q[pi[i]] << std::endl;
                std::cout << "q_{pi[i+1]} = " << D_q[pi[i+1]] << std::endl << std::endl;
                correct = false;
            }
        }
    }

    uint_t i = 0;
    uint_t j = 0;
    uint32_t e;

    while (i < k && j < k) {
        while (i < k && D_p[i] < D_q[pi[j]]) {
            i++;
        }

        if (D_p[i] < D_q[pi[j+1]]) {
            e = 1;
            
            while (i+1 < k && D_p[i+1] < D_q[pi[j+1]]) {
                i++;
                e++;

                if (e >= 2*a) {
                    std::cout << "a-heavy output interval:" << std::endl;
                    std::cout << "i = " << i-2*a+1 << std::endl;
                    std::cout << "j = " << j << std::endl;
                    std::cout << "q_{pi[j]} = " << D_q[pi[j]] << std::endl;
                    std::cout << "p_{i} = " << D_p[i-2*a+1] << std::endl;
                    std::cout << "p_{i+2a} = " << D_p[i] << std::endl;
                    std::cout << "q_{pi[j+1]} = " << D_q[pi[j+1]] << std::endl << std::endl;
                    correct = false;
                    j++;
                    break;
                }
            }

            i++;
        } else {
            j++;
        }
    }

    std::cout << "the interval sequence has " << (correct ? "" : "not ") << "been built correctly";
}