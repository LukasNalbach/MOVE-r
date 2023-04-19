/**
 * @brief locates the pattern P in T
 * @param P the pattern to locate in T
 * @param m the length of the pattern
 * @param Occ an empty array to write the positions of occurences of P in T to
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
void MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::locate(uint8_t* P, pos_t m, std::vector<pos_t>& Occ) {
    // Left interval limit of the suffix array interval.
    pos_t i_l = 0;
    // Right interval limit of the suffix array interval.
    pos_t i_r = n-1;

    /* position of the character in L' (and therefore the index of the input interval in M_LF),
    in which i_l lies. */
    idx_t x_l = 0;
    /* position of the character in L' (and therefore the index of the input interval in M_LF),
    in which i_r lies. */
    idx_t x_r = r_-1;

    /* m-j', if the j'-th iteration of the backward search is the last iteration, at whiches beginning
    P[j] != L'[x_r] holds. */
    pos_t m_m_j_ = m-1;
    /* x'_{r,j'}, that is the value of x_r before the LF queries and after the rank-select queries in the j'-th iteration
    of the backward search. */
    idx_t x__r_j_ = r_-1;

    // Temporary variable for P[j].
    uint8_t P_j;

    // Read P backwards from right to left.
    for (int64_t j=m-1; j>=0; j--) {
        // If the characters have been remapped internally, the pattern also has to be remapped.
        P_j = chars_mapped ? map_char[P[j]] : P[j];

        /* If P[j] != L[i_l] = L'[x_l], set i_l <- select(L,P[j],rank(L,P[j],i_l-1)+1) and adjust x_l,
        s.t. M_LF.p(x_l) <= i_l < M_LF.p(x_l+1). */

        if (P_j != L_(x_l)) {
            /* To do so, we can at first find the first run with character P[j] after the x_l-th run, save
            its index in x_l and set i_l <- M_LF.p(x_l). */
            x_l = RS_L_select_1[P_j].select(RS_L_rank_1[P_j].rank(x_l)+1);
            i_l = M_LF.p(x_l);
        }

        /* If P[j] != L[i_r] = L'[x_r], set i_r <- select(L,P[j],rank(L,P[j],i_r)) and adjust x_r, s.t.
        M_LF.p(x_r) <= i_r < M_LF.p(x_r+1). */

        if (P_j != L_(x_r)) {
            /* To do so, we can at first find the last run with character P[j] before the x_r-th run, save
            its index in x_r and set i_r <- M_LF.p(x_r+1)-1. */
            x_r = RS_L_select_1[P_j].select(RS_L_rank_1[P_j].rank(x_r));
            i_r = M_LF.p(x_r+1)-1;
            // update m-j, because P[j] != L'[x_r] held in this iteration.
            m_m_j_ = j;
            // update x__r_j_, because P[j] != L'[x_r] held in this iteration.
            x__r_j_ = x_r;
        }

        /* If the next suffix array interval [LF(i_l),LF(i_r)] is empty, then i_l > i_r,
        because LF(j) is monotonic for a fixed L[j], hence it suffices to check, whether
        i_l <= i_r holds. */

        // If the suffix array interval is empty, P does not occur in T, so return an empty array.
        if (i_l > i_r) {
            return;
        }
        
        // Set i_l <- LF(i_l) and i_r <- LF(i_r).
        M_LF.move(i_l,x_l);
        M_LF.move(i_r,x_r);
    }
    
    /* Initially, this is the value in the suffix array at the right interval limit of the suffix array
    interval of P. This variable is then used to iterate over the values of the suffix array in the suffix
    array interval of P from right (i_r) to left (i_l). */
    pos_t i_s;
    /* The index of the input interval in M_Phi, in which i_s currently lies, that is M_Phi.p(x_s) <= i_s
    < M_Phi.p(x_s+1) holds. */
    idx_t x_s;

    x_s = SA_idx(x__r_j_);
    i_s = M_Phi.p(x_s)+SA_offs(x__r_j_)-m_m_j_-1;

    /* More precisely, if SA_s[x'_{r,j'}]-(m-j+1) < M_Phi.p(SA_idx[x'_{r,j'}]) holds, i_s now lies in a input interval of
    M_Phi before the SA_idx[x'_{r,j'}]-th one, so we have to decrease x_s. To find the correct value for x_s, we perform an
    exponential search to the left over the input interval starting positions of M_Phi starting at SA_idx[x'_{r,j'}] = x_s. */
    if (i_s < M_Phi.p(x_s)) {
        // Current step size of the exponential search.
        idx_t s = 1;
        // Perform the first step.
        x_s -= 1;

        // Perform the exponential search.
        while (i_s < M_Phi.p(x_s)) {
            s *= 2;
            if (s > x_s) {
                x_s = 0;
            } else {
                x_s -= s;
            }
        }

        // Left limit of the binary search range.
        idx_t b = x_s;
        // Right limit of the binary search range.
        idx_t e = x_s + s - 1;
        // Candidate positon in the middle of the binary search range.
        idx_t c;

        // Perform the binary search.
        while (b != e) {
            c = b+(e-b)/2+1;
            if (i_s < M_Phi.p(c)) {
                e = c-1;
            } else {
                b = c;
            }
        }

        // Store the found posiiton in x_s.
        x_s = b;
    }

    /* Because the suffix array interval [i_l,i_r] of P has length i_r-i_l+1, there are i_r-i_l+1
    occurences of P in T. */
    Occ.reserve(i_r-i_l+1);

    /* Because i_s contains the value in the suffix array of the right limit of the suffix array
    interval of P, i_s is an occurence position of P in T. */
    Occ.emplace_back(i_s);
    
    // Perform i_r-i_l Phi queries with (i_s,x_s) and at each step, store i_s in Occ.
    for (pos_t i=1; i<=i_r-i_l; i++) {
        M_Phi.move(i_s,x_s);
        Occ.emplace_back(i_s);
    }
}