/**
 * @brief returns the number of occurences of P in T
 * @param P the pattern to locate in T
 * @param m the length of the pattern
 * @return the number of occurences of P in T
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
pos_t MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::count(uint8_t* P, pos_t m) {
    pos_t i_l = 0;
    pos_t i_r = n-1;

    idx_t x_l = 0;
    idx_t x_r = r_-1;

    uint8_t P_j;

    for (int64_t j=m-1; j>=0; j--) {
        P_j = chars_mapped ? map_char[P[j]] : P[j];

        if (P_j != L_(x_l)) {
            x_l = RS_L_select_1[P_j].select(RS_L_rank_1[P_j].rank(x_l)+1);
            i_l = M_LF.p(x_l);
        }

        if (P_j != L_(x_r)) {
            x_r = RS_L_select_1[P_j].select(RS_L_rank_1[P_j].rank(x_r));
            i_r = M_LF.p(x_r+1)-1;
        }

        // If the suffix array interval is empty, return 0.
        if (i_l > i_r) {
            return 0;
        }
        
        M_LF.move(i_l,x_l);
        M_LF.move(i_r,x_r);
    }

    // Return the size i_r-i_l+1 of the suffix array interval [i_l,i_r] of P.
    return i_r-i_l+1;
}