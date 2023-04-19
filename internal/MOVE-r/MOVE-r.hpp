#pragma once

#include <iostream>
#include <omp.h>
#include <sdsl/bit_vectors.hpp>
#include <misc/utils.hpp>
#include <mds/mds_phi.hpp>
#include <mds/mds_lf.hpp>
#include <coding/par_delta_vector.hpp>
#include <coding/par_huff_string.hpp>

enum support_mode {
	revert,
	revert_count,
	revert_count_locate
};

template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
class MOVE_r {
	static constexpr uint64_t size_entry_saidxoffs = bytes_idx + bytes_offs;
	static constexpr uint8_t pos_saoffs = bytes_idx;
	uint8_t* base_saidx;
	uint8_t* base_saoffs;

	private:
	support_mode support;
	uint64_t n = 0; // the length of T
	uint8_t sigma = 0;
	uint64_t r = 0; // r, the number of runs in L
	uint64_t r_ = 0; // r', the number of input/output intervals in M_LF.
	uint64_t r__ = 0; // r'', the number of input/output intervals in M_Phi.
	uint16_t a = 0; // Balancing parameter, restricts size to O(r*(a/(a-1))), 2 <= a.
	/* Maximum possible number of threads to use while reverting. If an index is built
	with some value p_r, then p threads can be used while reverting, wehere p <= p_r. */
	uint16_t p_r = 0;
	// Maximum number of threads to use during decoding and encoding
	uint16_t p_d = 0;
	bool chars_mapped = false; // Stores, whether the characters have been remapped internally.

	/* The Move Data Structure for LF. It also stores L', which can be accessed at
	position i with M_LF.S(i). */
	mds_lf<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> M_LF;
	mds_phi<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs> M_Phi; // The Move Data Structure for Phi.
	/* Stores SA_offs and SA_idx interleaved with each other,
	so SA_idxoffs[i] = (SA_offs[i],SA_idx[i]), for all i \in [0..r'-1] */
	uint8_t* SA_idxoffs = NULL;
	std::vector<sdsl::sd_vector<>> RS_L_bvs;
	std::vector<sdsl::sd_vector<>::rank_1_type> RS_L_rank_1;
	std::vector<sdsl::sd_vector<>::select_1_type> RS_L_select_1;
	/* Maps a character c in T to a value in the range [2..2+\sigma], where
	map_char[c] < map_char[c'] <=> c < c'holds for each two characters c,c'
	in T, to preserve the order of the characters among T. */
	std::vector<uint8_t> map_char;
	/* Reverse mapping function map_char_rev = map_char^{-1}. */
	std::vector<uint8_t> map_char_rev;
	/* [0..p_r-1] D_e[j] = i, s.t. L[i] = T[(j+1) * \lfoor (n-1)/p_r \rfloor], for j \in
	[0..p_r-2] and D_e[p_r-1] = 0, because L[n-2] = T[0]. */
	std::vector<uint64_t> D_e;

	inline void set_SA_idx(idx_t x, idx_t idx) {
		if constexpr (sizeof(idx_t) == bytes_idx) {
            *(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs)) = idx;
        } else {
            if constexpr (bytes_idx == 3) {
                *(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs))
                = (*(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs)) & 0xFF000000) | idx;
            } else {
                *(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs))
                = (*(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs)) & 0xFFFFFF0000000000) | idx;
            }
        }
	}

	inline void set_SA_offs(idx_t x, pos_t offs) {
		if constexpr (sizeof(pos_t) == bytes_offs) {
			*(reinterpret_cast<pos_t*>(base_saoffs + x * size_entry_saidxoffs)) = offs;
		} else if constexpr (sizeof(pos_t) % bytes_offs == 0) {
			if constexpr (bytes_offs == 1) {
				*(reinterpret_cast<uint8_t*>(base_saoffs + x * size_entry_saidxoffs)) = static_cast<uint8_t>(offs);
			} else if constexpr (bytes_offs == 2) {
				*(reinterpret_cast<uint16_t*>(base_saoffs + x * size_entry_saidxoffs)) = static_cast<uint16_t>(offs);
			} else {
				*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs)) = static_cast<uint32_t>(offs);
			}
        } else {
            if constexpr (bytes_offs == 3) {
                if constexpr (sizeof(pos_t) == 4) {
                    *(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))
                    = ((*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0xFF000000) | offs;
                } else {
                    *(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))
                    = ((*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0xFF000000) | static_cast<uint32_t>(offs);
                }
            } else {
                *(reinterpret_cast<uint64_t*>(base_saoffs + x * size_entry_saidxoffs))
                = ((*(reinterpret_cast<uint64_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0xFFFFFF0000000000) | offs;
            }
        }
	}

	/**
	 * @brief 
	 * @param x
	 * @return 
	 */
	inline idx_t SA_idx(idx_t x) {
		if constexpr (sizeof(idx_t) == bytes_idx) {
            return *(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs));
        } else {
            if constexpr (bytes_idx == 3) {
                return (*(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs))) & 0x00FFFFFF;
            } else {
                return (*(reinterpret_cast<idx_t*>(base_saidx + x * size_entry_saidxoffs))) & 0x000000FFFFFFFFFF;
            }
        }
	}

	/**
	 * @brief 
	 * @param x
	 * @return 
	 */
	inline pos_t SA_offs(idx_t x) {
		if constexpr (sizeof(pos_t) == bytes_offs) {
			return *(reinterpret_cast<pos_t*>(base_saoffs + x * size_entry_saidxoffs));
		} else if constexpr (sizeof(pos_t) % bytes_offs == 0) {
			if constexpr (bytes_offs == 1) {
				return static_cast<pos_t>(*(reinterpret_cast<uint8_t*>(base_saoffs + x * size_entry_saidxoffs)));
			} else if constexpr (bytes_offs == 2) {
				return static_cast<pos_t>(*(reinterpret_cast<uint16_t*>(base_saoffs + x * size_entry_saidxoffs)));
			} else {
				return static_cast<pos_t>(*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs)));
			}
        } else {
            if constexpr (bytes_offs == 3) {
				if constexpr (sizeof(pos_t) == 4) {
                    return (*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0x00FFFFFF;
                } else {
                    return static_cast<uint64_t>((*(reinterpret_cast<uint32_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0x00FFFFFF);
                }
            } else {
                return (*(reinterpret_cast<uint64_t*>(base_saoffs + x * size_entry_saidxoffs))) & 0x000000FFFFFFFFFF;
            }
        }
	}

	/**
	 * @brief 
	 * @param x
	 */
	inline uint8_t L_(idx_t x) {
		return M_LF.c(x);
	}

	uint64_t get_size_rs_l_() {
		uint64_t size_rs_l_ = 0;

		if (support != support_mode::revert) {
			for (uint16_t i=0; i<256; i++) {
				size_rs_l_ += sdsl::size_in_bytes(RS_L_bvs[i]);
				size_rs_l_ += sdsl::size_in_bytes(RS_L_rank_1[i]);
				size_rs_l_ += sdsl::size_in_bytes(RS_L_select_1[i]);
			}
		}

		return size_rs_l_;
	}

	public:
	/**
	 * @brief returns n
	 * @return n
	 */
	inline uint64_t get_n() {
		return n;
	}

	inline uint8_t get_sigma() {
		return sigma;
	}

	/**
	 * @brief returns r
	 * @return r
	 */
	inline uint64_t get_r() {
		return r;
	}

	/**
	 * @brief returns r'
	 * @return r'
	 */
	inline uint64_t get_r_() {
		return r_;
	}

	/**
	 * @brief returns r''
	 * @return r''
	 */
	inline uint64_t get_r__() {
		return r__;
	}

	/**
	 * @brief returns a
	 * @return a
	 */
	inline uint16_t get_a() {
		return a;
	}

	/**
	 * @brief returns p_r
	 * @return p_r
	 */
	inline uint16_t get_pd() {
		return p_r;
	}

	uint64_t get_size() {
		return
			73+ // variables
			(chars_mapped?2*256:0)+ // map_char & map_char^-1
			p_r*8+ // D_e
			M_LF.get_size()-(r_+1)+ // M_LF
			r_+1+ // L'
			get_size_rs_l_()+ // RS_L'
			M_Phi.get_size()+ // M_Phi
			r_*bytes_idx+ // SA_idx
			r_*bytes_offs; // SA_offs
	}

	void log_data_structure_sizes() {
		std::cout << "M_LF: " << format_size(M_LF.get_size()-(r_+1)) << std::endl;
		std::cout << "L': " << format_size(r_+1) << std::endl;
		
		if (support != support_mode::revert) {
			std::cout << "RS_L': " << format_size(get_size_rs_l_()) << std::endl;

			if (support == revert_count_locate) {
				std::cout << "M_Phi: " << format_size(M_Phi.get_size()) << std::endl;
				std::cout << "SA_idx: " << format_size(r_*bytes_idx) << std::endl;
				std::cout << "SA_offs: " << format_size(r_*bytes_offs) << std::endl;
			}
		}
	}

	void log_data_structure_sizes(std::ofstream &mf_idx) {
		mf_idx << " size_m_lf=" << M_LF.get_size()-(r_+1);
		mf_idx << " size_l_=" << r_+1;
		
		if (support != support_mode::revert) {
			mf_idx << " size_rs_l_=" << get_size_rs_l_();

			if (support == revert_count_locate) {
				mf_idx << " size_m_phi=" << M_Phi.get_size();
				mf_idx << " size_sa_idx=" << r_*bytes_idx;
				mf_idx << " size_sa_offs=" << r_*bytes_offs;
			}
		}
	}
	
	/**
	 * @brief 
	 * @param x 
	 * @return 
	 */
	inline uint8_t get_L_(idx_t x) {
		return L_(x);
	}

	/**
	 * @brief 
	 * @param x 
	 * @return 
	 */
	uint8_t get_L(pos_t i) {
		// x = max x' \in [1,r']: M_LF.p(x') <= i
		idx_t x;

		// Left interval limit of the binary search.
		idx_t l = 0;
		// Right interval limit of the binary search.
		idx_t r = r_-1;
		// Candidate position for the binary search.
		idx_t m;

		while (l != r) {
			m = l+(r-l)/2;
			if (M_LF.p(m) < i) {
				l = m+1;
			} else {
				r = m;
			}
		}

		x = l;

		// M_LF.p(x) <= i < M_LF.p(x+1) => L[i] = L'[x]
		return L_(x);
	}

	pos_t get_SA(pos_t i) {
		// x = max x' \in [1,r']: M_LF.p(x') <= i
		idx_t x;

		// Left interval limit of the binary search.
		idx_t l = 0;
		// Right interval limit of the binary search.
		idx_t r = r_-1;
		// Candidate position for the binary search.
		idx_t m;

		while (l != r) {
			m = l+(r-l)/2;
			if (M_LF.p(m) < i) {
				l = m+1;
			} else {
				r = m;
			}
		}

		x = l;

		if (i == M_LF.p(x)) {
			return M_Phi.p(SA_idx(x))+SA_offs(x);
		}

		x++;
		pos_t j = M_LF.p(x);

		// M_Phi.p(x_s) <= i_s < M_Phi.p(x_s+1).
		idx_t x_s = SA_idx(x);
		// i_s = SA[j]
		pos_t i_s = M_Phi.p(x_s)+SA_offs(x);

		// In each iteration, i_s = SA[j] = \Phi^{i-j}(SA[i]) holds.
		while (j > i) {
			// Set i_s = \Phi(i_s)
			M_Phi.move(i_s,x_s);
			j--;
		}

		// Since j = i, now i_s = SA[i] holds.
		return i_s;
	}

	MOVE_r() = default;

	MOVE_r(const MOVE_r&) = delete;

	~MOVE_r() {
		if (SA_idxoffs == NULL) {
			free(SA_idxoffs);
			SA_idxoffs = NULL;
		}
	}

	/**
	 * @brief reads a serialized r-index from an input stream
	 * @param index_file an input stream storing a serialized r-index
	 * @param p the number of threads to use (default: all threads)
	 */
	void load(std::ifstream& index_file, support_mode support, bool log, uint16_t p, std::ofstream &mf);

	/**
	 * @brief writes the original text to the string T
	 * @param p the number of threads to use, must not exceed the maximum
	 * number ofthreads to use that has been specified when the index was built
	 * @return T
	 */
	uint8_t* revert(uint16_t p);

	/**
	 * @brief returns the number of occurences of P in T
	 * @param P the pattern to locate in T
	 * @param m the length of the pattern
	 * @return the number of occurences of P in T
	 */
	pos_t count(uint8_t* P, pos_t m);

	/**
	 * @brief locates the pattern P in T
	 * @param P the pattern to locate in T
	 * @param m the length of the pattern
	 * @param Occ an empty array to write the positions of occurences of P in T to
	 */
	void locate(uint8_t* P, pos_t m, std::vector<pos_t>& Occ);
};

#include "MOVE-r/load.cpp"
#include "MOVE-r/revert.cpp"
#include "MOVE-r/count.cpp"
#include "MOVE-r/locate.cpp"