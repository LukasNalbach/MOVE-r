#pragma once

#include <coding/par_delta_vector.hpp>
#include <coding/par_var_width_vector.hpp>

/**
 * @brief stores a bijective function f_I : [0..n-1] -> [0..n-1] as a balanced disjoint
 *        interval sequence B(I)[0..k'], supports calculation of f_I(i) = i', with i
 *        in [0..n-1], by calculating Move(i,x) = (i',x')
 * @tparam 
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
class mds_phi {
    protected:
    static constexpr uint8_t pos_idx = bytes_pos;
    static constexpr uint8_t pos_offs = bytes_pos + bytes_idx;
    uint8_t* base_pos = NULL;
    uint8_t* base_idx = NULL;
    uint8_t* base_offs = NULL;
    
    virtual inline constexpr uint64_t size_entry() { return bytes_pos + bytes_idx + bytes_offs; }

    uint64_t n = 0; // maximum value, n = p_{k'-1} + d_{k'-1}, k' <= n
    uint64_t k = 0; // number of intervals in the balanced disjoint inteval sequence B(I), 0 < k'

    /** @brief stores the balanced disjoint inteval sequence B(I) interleaved with D_idx
     *  = ((p_0,q_0,D_idx[0],S[1]),(p_1,q_1,D_idx[1],S[1]),...,(p_{k'-1},q_{k'-1},D_idx[k'-1])),
     *  (n,n,-1,0), where D_idx[j] = i <=> q_j in [p_i, p_i + d_i - 1], with i,j in [0..k'-1] */
    uint8_t* data = NULL;

    public:
    mds_phi() = default;

    ~mds_phi() {
        if (data != NULL) {
            free(data);
            data = NULL;
        }
    }

    mds_phi(const mds_phi&) = delete;

    virtual uint64_t get_size() {
        return
            58+ // variables
            (k+1)*size_entry(); // data
    }

    /**
     * @brief returns n
     * @return n = p_{k'-1} + d_{k'-1}
     */
    inline uint64_t get_n() {
        return n;
    }

    /**
     * @brief returns the number k' of intervals in the disjoint interval sequnece
     * @return number k' of intervals in the disjoint interval sequnece
     */
    inline uint64_t get_k() {
        return k;
    }

    /**
     * @brief resizes the move data structure to size k
     * @param n maximum value
     * @param k size
     */
    virtual void resize(uint64_t n, uint64_t k) {
        this->n = n;
        this->k = k;
        data = (uint8_t*) malloc((k+1)*size_entry());
        base_pos = data;
        base_idx = data + pos_idx;
        base_offs = data + pos_offs;
        set_p(k,n);
        set_offs(k,0);
        set_idx(k,k);
    }

    /**
     * @brief returns p_i
     * @param i in [0..k']
     * @return p_i
     */
    inline pos_t p(idx_t i) {
        if constexpr (sizeof(pos_t) == bytes_pos) {
            return *(reinterpret_cast<pos_t*>(base_pos + i * size_entry()));
        } else {
            if constexpr (bytes_pos == 3) {
                return (*(reinterpret_cast<pos_t*>(base_pos + i * size_entry()))) & 0xFFFF00;
            } else {
                return (*(reinterpret_cast<pos_t*>(base_pos + i * size_entry()))) & 0x000000FFFFFFFFFF;
            }
        }
    }

    /**
     * @brief 
     * @param i 
     * @return 
     */
    inline pos_t offs(idx_t i) {
        if constexpr (sizeof(pos_t) == bytes_offs) {
			return *(reinterpret_cast<pos_t*>(base_offs + i * size_entry()));
		} else if constexpr (sizeof(pos_t) % bytes_offs == 0) {
			if constexpr (bytes_offs == 1) {
				return static_cast<pos_t>(*(reinterpret_cast<uint8_t*>(base_offs + i * size_entry())));
			} else if constexpr (bytes_offs == 2) {
				return static_cast<pos_t>(*(reinterpret_cast<uint16_t*>(base_offs + i * size_entry())));
			} else {
				return static_cast<pos_t>(*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry())));
			}
        } else {
            if constexpr (bytes_offs == 3) {
                if constexpr (sizeof(pos_t) == 4) {
                    return (*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))) & 0x00FFFFFF;
                } else {
                    return static_cast<uint64_t>((*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))) & 0x00FFFFFF);
                }
            } else {
                return (*(reinterpret_cast<uint64_t*>(base_offs + i * size_entry()))) & 0x000000FFFFFFFFFF;
            }
        }
    }

    /**
     * @brief returns q_i
     * @param i in [0..k']
     * @return q_i
     */
    inline pos_t q(idx_t i) {
        return p(idx(i))+offs(i);
    }

    /**
     * @brief returns D_idx[i]
     * @param i in [0..k']
     * @return D_idx[i]
     */
    inline idx_t idx(idx_t i) {
        if constexpr (sizeof(idx_t) == bytes_idx) {
            return *(reinterpret_cast<idx_t*>(base_idx + i * size_entry()));
        } else {
            if constexpr (bytes_idx == 3) {
                return (*(reinterpret_cast<idx_t*>(base_idx + i * size_entry()))) & 0x00FFFFFF;
            } else {
                return (*(reinterpret_cast<idx_t*>(base_idx + i * size_entry()))) & 0x000000FFFFFFFFFF;
            }
        }
    }

    inline void set_p(idx_t i, pos_t p) {
        if constexpr (sizeof(pos_t) == bytes_pos) {
            *(reinterpret_cast<pos_t*>(base_pos + i * size_entry())) = p;
        } else {
            if constexpr (bytes_pos == 3) {
                *(reinterpret_cast<pos_t*>(base_pos + i * size_entry()))
                = (*(reinterpret_cast<pos_t*>(base_pos + i * size_entry())) & 0xFF000000) | p;
            } else {
                *(reinterpret_cast<pos_t*>(base_pos + i * size_entry()))
                = (*(reinterpret_cast<pos_t*>(base_pos + i * size_entry())) & 0xFFFFFF0000000000) | p;
            }
        }
    }

    inline void set_offs(idx_t i, pos_t offs) {
        if constexpr (sizeof(pos_t) == bytes_offs) {
			*(reinterpret_cast<pos_t*>(base_offs + i * size_entry())) = offs;
		} else if constexpr (sizeof(pos_t) % bytes_offs == 0) {
			if constexpr (bytes_offs == 1) {
				*(reinterpret_cast<uint8_t*>(base_offs + i * size_entry())) = static_cast<uint8_t>(offs);
			} else if constexpr (bytes_offs == 2) {
				*(reinterpret_cast<uint16_t*>(base_offs + i * size_entry())) = static_cast<uint16_t>(offs);
			} else {
				*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry())) = static_cast<uint32_t>(offs);
			}
        } else {
            if constexpr (bytes_offs == 3) {
                if constexpr (sizeof(pos_t) == 4) {
                    *(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))
                    = ((*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))) & 0xFF000000) | offs;
                } else {
                    *(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))
                    = ((*(reinterpret_cast<uint32_t*>(base_offs + i * size_entry()))) & 0xFF000000) | static_cast<uint32_t>(offs);
                }
            } else {
                *(reinterpret_cast<uint64_t*>(base_offs + i * size_entry()))
                = ((*(reinterpret_cast<uint64_t*>(base_offs + i * size_entry()))) & 0xFFFFFF0000000000) | offs;
            }
        }
    }

    inline void set_idx(idx_t i, idx_t idx) {
        if constexpr (sizeof(idx_t) == bytes_idx) {
            *(reinterpret_cast<idx_t*>(base_idx + i * size_entry())) = idx;
        } else {
            if constexpr (bytes_idx == 3) {
                *(reinterpret_cast<idx_t*>(base_idx + i * size_entry()))
                = (*(reinterpret_cast<idx_t*>(base_idx + i * size_entry())) & 0xFF000000) | idx;
            } else {
                *(reinterpret_cast<idx_t*>(base_idx + i * size_entry()))
                = (*(reinterpret_cast<idx_t*>(base_idx + i * size_entry())) & 0xFFFFFF0000000000) | idx;
            }
        }
    }

    /**
     * @brief calculates the move query Move(I,i,x) = (i',x') by changing ix = (i,x)
     *        to ix' = (i',x'), with i' = f_I(i) and i' in [p_x', p_x' + d_x' - 1]
     * @param i in [0,n-1]
     * @param x in [0..k'-1], where i in [p_x, p_x + d_x - 1]
     */
    inline void move(pos_t& i, idx_t& x) {
        i = q(x)+(i-p(x));
        x = idx(x);
        while (i >= p(x+1)) {
            x++;
        }
    }

    /**
     * @brief creates a move datastructure from an input stream
     * @param in input stream
     * @param p_ the number of threads to use
     */
    void load(std::ifstream &in, bool log, uint16_t p_) {
        auto time = now();
        in.read((char*)&n,sizeof(uint64_t));
        in.read((char*)&k,sizeof(uint64_t));
        uint16_t p_d;
        in.read((char*)&p_d,sizeof(uint16_t));
        std::vector<uint64_t> D_ps(p_d);
        in.read((char*)&D_ps[0],p_d*sizeof(uint64_t));
        resize(n,k);
        if (log) std::cout << std::endl << "reconstructing D_p" << std::flush;

        bool var_width_encode;
        in >> var_width_encode;
        if (var_width_encode) {
            par_var_width_vector<pos_t> D_p_pvwv;
            D_p_pvwv.load_and_decode_buffered(in,[this](uint64_t i, pos_t p){set_p(i,p);},0,k-1,p_);
        } else {
            par_delta_vector<pos_t> D_p_pdv;
            D_p_pdv.load_and_decode_buffered(in,[this](uint64_t i, pos_t p){set_p(i,p);},0,k-1,p_);
        }

        #pragma omp parallel for num_threads(p_)
        for (uint16_t ip_d=0; ip_d<p_d; ip_d++) {
            idx_t b = ip_d*(k/p_d);
            idx_t e = ip_d == p_d-1 ? k-1 : (ip_d+1)*(k/p_d)-1;
            pos_t p_cur = D_ps[ip_d];
            pos_t tmp;

            for (idx_t i=b; i<=e; i++) {
                tmp = p(i);
                set_p(i,p_cur);
                p_cur += tmp;
            }
        }

        D_ps.clear();
        D_ps.shrink_to_fit();

        if (log) {
			time = log_runtime(time);
            std::cout << "reconstructing D_offs"  << std::flush;
		}

        in >> var_width_encode;
        if (var_width_encode) {
            par_var_width_vector<pos_t> D_offs_pvwv;
            D_offs_pvwv.load_and_decode_buffered(in,[this](uint64_t i, pos_t offs){set_offs(i,offs);},0,k-1,p_);
        } else {
            par_delta_vector<pos_t> D_offs_pdv;
            D_offs_pdv.load_and_decode_buffered(in,[this](uint64_t i, pos_t offs){set_offs(i,offs);},0,k-1,p_);
        }

        if (log) {
			time = log_runtime(time);
            std::cout << "reconstructing D_idx"  << std::flush;
		}

        par_var_width_vector<idx_t> D_idx_pvwv;
        D_idx_pvwv.load_and_decode_buffered(in,[this](uint64_t i, idx_t idx){set_idx(i,idx);},0,k-1,p_);

        if (log) {
            time = log_runtime(time);
            std::cout << "move data structure reconstructed" << std::flush;
		}
    }
};

/**
 * @brief writes a move_datastructure to ab output stream
 * @param n 
 * @param k 
 * @param D_p 
 * @param pi 
 * @param out 
 * @param p 
 * @param p_d 
 */
template <typename uint_t>
static void serialize_mds(
    uint64_t n, uint64_t k,
    uint_t*& D_p,
    std::pair<uint_t,uint_t>*& D_offsidx,
    uint8_t bytes_pos,
    uint8_t bytes_idx,
    uint8_t bytes_offs,
    std::ofstream &out,
    uint64_t& time_write_index_file,
    uint64_t& time_encode,
    uint16_t p,
    uint16_t p_d,
    bool log,
    std::ofstream &mf_idx,
    std::ofstream &mf_mds
) {
    out.write((char*)&n,sizeof(uint64_t));
    out.write((char*)&k,sizeof(uint64_t));
    out.write((char*)&p_d,sizeof(uint16_t));

    std::vector<uint64_t> D_ps(p_d);

    for(uint16_t i=0; i<p_d; i++) {
        D_ps[i] = D_p[i*(k/p_d)];
    }

    out.write((char*)&D_ps[0],p_d*sizeof(uint64_t));
    D_ps.clear();
    auto time = now();
    if (log) std::cout << std::endl << "delta-encoding input interval lengths";

    par_delta_vector<uint_t> D_p_pdv;
    D_p_pdv.build([&D_p](uint64_t i){return D_p[i+1]-D_p[i];},0,k-1,p_d,p);
    float comp_ratio = (k*bytes_offs)/(float)(D_p_pdv.get_bit_size()/8);
    bool var_width_encode = comp_ratio < 1;
    out << var_width_encode;

    if (log) {
        time_encode += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "compression ratio: " << comp_ratio;
        if (comp_ratio < 1) {
            std::cout << " < 1" << std::endl << "=> variable-width-encoding input interval lengths";
        } else {
            std::cout << ", " << format_size(k*bytes_offs) << " -> " << format_size(D_p_pdv.get_bit_size()/8);
        };
    }

    par_var_width_vector<uint_t> D_p_pvwv;
    if (var_width_encode) {
        D_p_pdv.clear();
        D_p_pvwv.build([&D_p](uint64_t i){return D_p[i+1]-D_p[i];},0,k-1,bytes_offs*8,p_d,p);
        comp_ratio = (k*bytes_offs)/(float)(D_p_pvwv.get_bit_size()/8);
    }

    if (log) {
        time_encode += time_diff_ns(time,now());
        if (var_width_encode) {
            time = log_runtime(time);
            std::cout << "compression ratio: " << comp_ratio;
            std::cout << ", " << format_size(k*bytes_offs) << " -> " << format_size(D_p_pvwv.get_bit_size()/8);
        }
        std::cout << std::endl << "serializing input interval lengths";
    }

    if (var_width_encode) {
        D_p_pvwv.serialize(out);
        D_p_pvwv.clear();
    } else {
        D_p_pdv.serialize(out);
        D_p_pdv.clear();
    }

    if (log) {
        time_write_index_file += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "delta-encoding D_offs";
    }

    par_delta_vector<uint_t> D_offs_pdv;
    D_offs_pdv.build([&D_offsidx](uint64_t i){return D_offsidx[i].first;},0,k-1,p_d,p);
    comp_ratio = (k*bytes_offs)/(float)(D_offs_pdv.get_bit_size()/8);
    var_width_encode = comp_ratio < 1;
    out << var_width_encode;
    
    if (log) {
        time_encode += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "compression ratio: " << comp_ratio;
        if (var_width_encode) {
            std::cout << " < 1" << std::endl << "=> variable-width-encoding D_offs";
        } else {
            std::cout << ", " << format_size(k*bytes_offs) << " -> " << format_size(D_offs_pdv.get_bit_size()/8);
        };
    }

    par_var_width_vector<uint_t> D_offs_pvwv;
    if (var_width_encode) {
        D_offs_pdv.clear();
        D_offs_pvwv.build([&D_offsidx](uint64_t i){return D_offsidx[i].first;},0,k-1,bytes_offs*8,p_d,p);
        comp_ratio = (k*bytes_offs)/(float)(D_offs_pvwv.get_bit_size()/8);
    }

    if (log) {
        time_encode += time_diff_ns(time,now());
        if (var_width_encode) {
            time = log_runtime(time);
            std::cout << "compression ratio: " << comp_ratio;
            std::cout << ", " << format_size(k*bytes_offs) << " -> " << format_size(D_offs_pvwv.get_bit_size()/8);
        }
        std::cout << std::endl << "serializing D_offs";
    }

    if (var_width_encode) {
        D_offs_pvwv.serialize(out);
        D_offs_pvwv.clear();
    } else {
        D_offs_pdv.serialize(out);
        D_offs_pdv.clear();
    }
    
    if (log) {
        time_write_index_file += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "variable-width-encoding D_idx";
    }

    par_var_width_vector<uint_t> D_idx_pvwv;
    D_idx_pvwv.build([&D_offsidx](uint64_t i){return D_offsidx[i].second;},0,k-1,(uint8_t)std::ceil(log2(k-1)),p_d,p);
    comp_ratio = (k*bytes_idx)/(float)(D_idx_pvwv.get_bit_size()/8);

    if (log) {
        time_encode += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "compression ratio: " << comp_ratio;
        std::cout << ", " << format_size(k*bytes_idx) << " -> " << format_size(D_idx_pvwv.get_bit_size()/8);
        std::cout << std::endl << "serializing D_idx";
    }

    D_idx_pvwv.serialize(out);
    D_idx_pvwv.clear();

    if (log) {
        time_write_index_file += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "move data structure serialized";
    }
}