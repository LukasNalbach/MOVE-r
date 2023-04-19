#pragma once

/**
 * @brief stores a bijective function f_I : [0..n-1] -> [0..n-1] as a balanced disjoint interval
 *        sequence B(I)[0..k'] and a string c[0..k'], supports calculation of f_I(i) = i', with i
 *        in [0..n-1], by calculating Move(i,x) = (i',x'), stores the balanced disjoint inteval
 *        sequence B(I) interleaved with D_idx = ((p_0,q_0,D_idx[0],c[1]),(p_1,q_1,D_idx[1],c[1])
 *        ,...,(p_{k'-1},q_{k'-1},D_idx[k'-1],c[k'-1])), (n,n,-1,0), where D_idx[j] = i <=> q_j
 *        in [p_i, p_i + d_i - 1], with i,j in [0..k'-1]
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
class mds_lf : public mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> {
    protected:
    static constexpr uint8_t pos_char = bytes_pos + bytes_idx + bytes_offs;
    uint8_t* base_c = NULL;

    inline constexpr uint64_t size_entry() override { return bytes_pos + bytes_idx + bytes_offs + 1; }

    public:
    /**
     * @brief creates an empty move datastructure
     */
    mds_lf() : mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>() {};

    mds_lf(const mds_lf&) = delete;

    /**
     * @brief resizes the move data structure to size k
     * @param n maximum value
     * @param k size
     */
    virtual void resize(uint64_t n, uint64_t k) override {
        mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::resize(n,k);
        base_c = (uint8_t*) mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::data + pos_char;
        c(k) = 0;
    }

    virtual uint64_t get_size() override {
        return
            67+ // variables
            (mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::k+1)*size_entry(); // data
    }

    /**
     * @brief returns c[i]
     * @param i in [0..k']
     * @return c[i]
     */
    inline uint8_t& c(idx_t i) {
        return *(base_c + i * size_entry());
    }
};