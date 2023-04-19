#pragma once

#include <malloc_count.h>
#include <libsais.h>
#include <libsais64.h>
#include <ips4o.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/io.hpp>
#include <mds/mds_builder/mds_builder.hpp>
#include <MOVE-r/MOVE-r.hpp>
#include <coding/par_delta_vector.hpp>
#include <coding/par_huff_string.hpp>
#include <coding/par_var_width_vector.hpp>

enum memory_usage_mode {
	low,
	high,
	automatic
};

/**
 * @brief Builds an r-index
 * @param T the text, to build the r-index of, has to be of size at most
 * INT_MAX-1 bytes for 32-bit, and of at least INT_MAX bytes for 64-bit
 * @param p the number of threads to use (default: all threads)
 * @param index_file the file to write the index to
 * @param support support mode
 * @param memory_usage_mode
 * @param p_r maximum possible number of threads to use while reverting. If an
 * index is built with some value p_r, then p' threads can be used while reverting,
 * wehere p' <= p_r.
 * @param p_r maximum possible number of threads to use while encoding and reconstructing the index
 * @param a balancing parameter, restricts size to O(r*(a/(a-1))), 2 <= a (default: 8)
 * @param epsilon
 * @param log determines, whether to print log messages (default: false)
 * @param mf_idx file to write runtime and peak memory usage data of the index construction
 * to (default: NULL)
 * @param mf_mds file to write runtime and peak memory usage data of the construction of the
 * move-datastructures to (default: NULL)
 * @param name_textfile name of the text file to build the index of
 */
template <typename uint_t, typename int_t>
static void MOVE_r_build(
    uint8_t*& T,
    uint64_t n,
    std::ofstream &index_file,
    support_mode support,
    memory_usage_mode memory_usage_mode,
    uint16_t p,
    uint16_t p_r,
    uint16_t p_d,
    uint16_t a,
    double epsilon,
    bool log,
    std::ofstream &mf_idx,
    std::ofstream &mf_mds,
    std::string name_textfile
) {
    auto time = now();
    auto time_start = time;
    uint64_t time_store_on_disk = 0;
    uint64_t time_write_index_file = 0;
    uint64_t time_encode = 0;
    uint64_t size_index_in_ram = 73;
    
    if (p > 1 && 1000*p > n) {
        p = std::max((uint64_t)1,n/1000);
        if (log) std::cout << "warning: p > n/1000, setting p to n/1000 ~ " << std::to_string(p) << std::endl;
    }

    if (p > 1 && 1000*p_r > n) {
        p_r = std::max((uint64_t)1,n/1000);
        if (log) std::cout << "warning: p_r > n/1000, setting p_r to n/1000 ~ " << std::to_string(p_r) << std::endl;
    }

    if (p > 1 && 1000*p_d > n) {
        p_d = std::max((uint64_t)1,n/1000);
        if (log) std::cout << "warning: p_d > n/1000, setting p_d to n/1000 ~ " << std::to_string(p_d) << std::endl;
    }

    omp_set_num_threads(p);
    if (log) std::cout << "preprocessing T" << std::flush;

    bool chars_mapped = false;
    // contains_char_thr[i_p][c] = true <=> thread i_p found the character c in its section of T.
    std::vector<std::vector<bool>> contains_char_thr(p,std::vector<bool>(256,false));

    // Iterate over T and report the occurence of each found character in T in contains_char_thr.
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<n-1; i++) {
        contains_char_thr[omp_get_thread_num()][T[i]] = true;
    }

    std::vector<bool> contains_char(256,false);

    /* Combine the results of each thread's sections in contains_char_thr[0..i_p-1][0..255]
    into contains_char[0..255]. */
    for (uint16_t i=0; i<256; i++) {
        for (uint16_t i_p=0; i_p<p; i_p++) {
            if (contains_char_thr[i_p][i]) {
                contains_char[i] = true;
                break;
            }
        }
    }

    contains_char_thr.clear();
    std::vector<uint8_t> map_char;

    // The number of distinct characters in T.
    uint8_t sigma = 1;

    // Count the number of ones in contains_char[0..255] in sigma.
    for (uint16_t i=0; i<256; i++) {
        if (contains_char[i]) {
            sigma++;
        }
    }

    // If T contains 0 or 1, we have to remap the characters in T.
    if (contains_char[0] || contains_char[1]) {
        if (sigma > 253) {
            /* If T contains more than 253 distinct characters, we cannot remap them into the 
            range [0..255] without using 0 or 1, hence we cannot build the r-index for T. */
            std::cout << "Error: the input contains more than 253 distinct characters" << std::endl;
            return;
        }

        // Else, we build the mapping function map_char

        chars_mapped = true;
        map_char.resize(256,0);

        // the character, to map the currently next largest character in T to
        uint16_t j = 2;

        /* To preserve the order among characters in T, we start by mapping smallest
        character in T to 2, the second smallest to 3, ... . */
        for (uint16_t i=0; i<256; i++) {
            if (contains_char[i]) {
                map_char[i] = j;
                j++;
            }
        }

        // Apply the mapping function map_char to T.
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T[i] = map_char[T[i]];
        }
    }

    if (chars_mapped) {
        size_index_in_ram += 2*256;
    }

    if (log) {
        if (mf_idx.is_open()) mf_idx << " time_preprocess_t=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "building SA" << std::flush;
    }

    int32_t fs_SA = 6*sigma;

    // The suffix array.
    int_t* SA = (int_t*) malloc((n+fs_SA)*sizeof(int_t));

    // Choose the correct suffix srray construction algorithm.
    if constexpr (std::is_same<int_t,int32_t>::value) {
        if (p == 1) {
            libsais(&T[0],&SA[0],n,fs_SA,NULL);
        } else {
            libsais_omp(&T[0],&SA[0],n,fs_SA,NULL,p);
        }
    } else {
        if (p == 1) {
            libsais64(&T[0],&SA[0],n,fs_SA,NULL);
        } else {
            libsais64_omp(&T[0],&SA[0],n,fs_SA,NULL,p);
        }
    }
    
    if (log) {
        if (mf_idx.is_open()) mf_idx << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "building L and C" << std::flush;
    }

    // The BWT.
    uint8_t* L = (uint8_t*) malloc(n);

    // [0..p][0..255]
    std::vector<std::vector<uint_t>> C(p+1);
    // [0..p]
    std::vector<uint_t> r_p(p+1);
    // \lfloor (n-1)/p_r \rfloor
    uint_t nm1pr = (n-1)/p_r;

    // For now, we want D_e to have size p_r+1, which is one more than its final size.
    std::vector<uint64_t> D_e(p_r+1);

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        C[i_p].resize(256,0);

        // Iteration range start position of thread i_p.
        uint_t b = i_p*(n/p);
        // Iteration range end position of thread i_p.
        uint_t e = i_p == p-1 ? n-1 : (i_p+1)*(n/p)-1;

        // Build L in the range [b..e].
        for (uint_t i=b; i<=e; i++) {
            L[i] = SA[i] == 0 ? T[n-1] : T[SA[i]-1];
        }

        /* Store in C[i_p][c] the number of occurences of c in L[b..e], for each c \in [0..255].
        Store in r_p[i_p] the number of runs starting in L[b..e]. */

        // We add an artificial run starting at b.
        C[i_p][L[b]]++;
        r_p[i_p] = 1;

        // Temporary variable for the position of L[i] in T.
        uint_t j;
        
        // Iterate over the range L[b..e]
        for (uint_t i=b+1; i<=e; i++) {

            // Calculate the position j of L[i] in T.
            j = SA[i] == 0 ? n-1 : SA[i]-1;
            
            /* If j is a multiple of (n-1)/p_r, store store i at D_e[j/((n-1)/p_r)].
            Because we may find a j = p_r*(n-1)/p_r, D_e must have size p_r+1. */
            if (j%nm1pr == 0) {
                D_e[j/nm1pr] = i;
            }

            // Count the occurence of L[i] in C[i_p][L[i]].
            C[i_p][L[i]]++;

            // If there is a run starting at position i, count it in r_p[i_p].
            if (L[i] != L[i-1]) {
                r_p[i_p]++;
            }
        }
    }

    if (support != revert_count_locate) {
        free(SA);
        SA = NULL;
    }

    /* Now we have D_e[i] = T[j * \lfoor (n-1)/p_r \rfloor], for all i \in [0..p_r]. We
    want to have D_e[i] = T[(j+1) * \lfoor (n-1)/p_r \rfloor], for all i \in [0..p_r-1]. */

    for (uint_t i=0; i<p_r; i++) {
        D_e[i] = D_e[i+1];
    }

    D_e.resize(p_r);
    D_e[p_r-1] = 0;
    size_index_in_ram += p_r*8;

    // Now we are done with D_e.
    
    free(T);
    T = NULL;
    C[p].resize(256,0);

    /* Now, C[i_p][c] is the number of occurences of c in L[b..e], where [b..e] is the range of the
    thread i_p \in [0..p-1]. Also, we have C[p][0..255] = 0. */

    /* We want to have C[i_p][c] = rank(L,c,b-1), where b is the iteration range start position of
    thread i_p \in [0..p-1]. Also, we want C[p][0..255] to be the C-array, that is C[p][c] stores
    the number of occurences of all smaller characters c' < c in L[0..n-1], for c \in [0..255]. */

    for (uint16_t i=1; i<p; i++) {
        for (uint16_t j=0; j<256; j++) {
            C[i][j] += C[i-1][j];
        }
    }

    /* Now, we have C[i_p][c] = rank(L,c,e), for each c \in [0..255] and i_p \in [0..p-1],
    where e is the iteration range end position of thread i_p. */

    for (uint16_t i=p; i>0; i--) {
        for (uint16_t j=0; j<256; j++) {
            C[i][j] = C[i-1][j];
        }
    }
    
    for (uint16_t j=0; j<256; j++) {
        C[0][j] = 0;
    }

    /* Now, we have C[i_p][c] = rank(L,c,b-1), for each c \in [0..255] and i_p \in [0..p-1],
    where b is the iteration range start position of thread i_p, so we are done with C[0..p-1][0..255].
    Also, we have C[p][c] = rank(L,c,n-1), for c \in [0..255]. */

    for (uint16_t i=255; i>0; i--) {
        C[p][i] = C[p][i-1];
    }

    C[p][0] = 0;

    // Now, we have C[p][c] = rank(L,c-1,n-1), for c \in [1..255], and C[p][0] = 0.

    for (uint16_t i=2; i<256; i++) {
        C[p][i] += C[p][i-1];
    }

    // Now we are done with C, since C[p] is the C-array.

    /* Now, r_p[i_p] stores the number of runs starting in the iteration range L[b..e] of thread
    i_p \in [0..p-1], and r_p[p] = 0. We want r_p[i_p] to store the number of runs starting before
    the iteration range start position b of thread i_p \in [0..p-1]. Also, we want r_p[p] to store
    the number of all runs. This way, we can build I_LF[r_p[i_p]..r_p[i_p+1]-1] with thread i_p \in [0..p-1]
    using the C-array in C[p] while using and updating the rank-function in C[i_p] on-the-fly. */

    for (uint16_t i=p; i>0; i--) {
        r_p[i] = r_p[i-1];
    }

    r_p[0] = 0;

    for (uint16_t i=2; i<=p; i++) {
        r_p[i] += r_p[i-1];
    }

    /* Now, r_p[i_p] stores the number of runs starting before the iteration range start position b of
    thread i_p \in [0..p-1] and r_p[p] stores the number of all runs, so we are done with r_p[0..p] */
    
    uint64_t r = r_p[p];

    if (log) {
        double n_r = std::round(100.0*(n/(double)r))/100.0;
        if (mf_idx.is_open()) {
            mf_idx << " time_build_l_c=" << time_diff_ns(time,now());
            mf_idx << " n=" << n;
            mf_idx << " sigma=" << std::to_string(sigma);
            mf_idx << " r=" << r;
        }
        time = log_runtime(time);
        std::cout << "n = " << n << ", sigma = " << std::to_string(sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }

    uint8_t support_uint;

    switch (support) {
        case revert:
            support_uint = 0;
            break;
        case revert_count:
            support_uint = 1;
            break;
        case revert_count_locate:
            support_uint = 2;
    }

    uint8_t bytes_pos = 0;
    uint8_t bytes_idx = 0;
    uint8_t bytes_offs = 0;

    for (uint8_t bytes=1; bytes<=5; bytes++) {
        if (n <= pow(2,8*bytes)-1) {
            bytes_pos = bytes;
            break;
        }
    }

    for (uint8_t bytes=1; bytes<=5; bytes++) {
        if (n/((pow(2,8*bytes)-1)*(double)r) <= epsilon) {
            bytes_offs = bytes;
            break;
        }
    }

    if (log) {
        std::cout << "building I_LF";
        if (support != revert_count_locate) std::cout << " and I_Phi";
        std::cout << std::flush;
    }

    std::pair<uint_t,uint_t>* I_LF = NULL;
    std::pair<uint_t,uint_t>* I_Phi = NULL;
    
    I_LF = (std::pair<uint_t,uint_t>*) malloc((r+1)*sizeof(std::pair<uint_t,uint_t>));
    I_LF[r] = std::pair<uint_t,uint_t>(n,n);

    if (support == revert_count_locate) {
        I_Phi = (std::pair<uint_t,uint_t>*) malloc((r+1)*sizeof(std::pair<uint_t,uint_t>));
        I_Phi[r] = std::pair<uint_t,uint_t>(n,n);
    }

    #pragma omp parallel num_threads(p)
    {
        // Index \in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b = i_p*(n/p);
        // Iteration range end position of thread i_p.
        uint_t e = i_p == p-1 ? n-1 : (i_p+1)*(n/p)-1;

        // Start position of the last-seen run.
        uint_t i_ = b;
        // Position in I_LF and M_Phi to write the next pairs to.
        uint_t j = r_p[i_p];

        // Build I_LF[r_p[i_p]..r_p[i_p+1]-1] and I_Phi[r_p[i_p]..r_p[i_p+1]-1]

        /* Each thread creates one input interval starting at r_p[b] and SA[b], respectively,
        whiches pairs have to be placed at LF[r_p[i_p]] and I_Phi[r_p[i_p]]. */

        I_LF[j] = std::pair<uint_t,uint_t>{i_,C[p][L[b]]+C[i_p][L[b]]};

        if (support == revert_count_locate) {
            I_Phi[j] = std::pair<uint_t,uint_t>{SA[i_],SA[i_ == 0 ? n-1 : i_-1]};
        }

        j++;

        // Iterate over the range L[b+1..e]
        for (uint_t i=b+1; i<=e; i++) {

            // Check, if a run starts at a position i
            if (L[i] != L[i-1]) {

                /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i-1),
                for each c \in [0..255]. */
                C[i_p][L[i-1]] += i-i_;

                // Update the position i_ of the last-seen run to i.
                i_ = i;

                /* Write the pair (i,LF(i)) to the next position j in I_LF, where
                LF(i) = C[L[i]] + rank(L,L[i],i-1) = C[p][L[i]] + C[i_p][L[i]]. */
                I_LF[j] = std::pair<uint_t,uint_t>{i,C[p][L[i]]+C[i_p][L[i]]};

                if (support == revert_count_locate) {
                    // Write the pair (SA[i],Phi(SA[i])) to the next position j in I_Phi.
                    I_Phi[j] = std::pair<uint_t,uint_t>{SA[i],SA[i-1]};
                }

                j++;
            }
        }
    }

    r_p.clear();
    C.clear();

    if (log) {
        std::string str;

        if (support == revert_count_locate) {
            str = " time_build_ilf_iphi=";
        } else {
            str = " time_build_ilf=";
        }

        if (mf_idx.is_open()) mf_idx << str << time_diff_ns(time,now());
    }

    if (log) time = log_runtime(time);

    if (support == revert_count_locate) {
        if (log) std::cout << "sorting I_Phi" << std::flush;

        // Sort I_Phi by the starting positions of its output intervals.
        auto comp_I_Phi = [](auto p1, auto p2) {return p1.first < p2.first;};
        // Choose the correct sorting algorithm.
        if (p > 1) {
            ips4o::parallel::sort(&I_Phi[0],&I_Phi[r],comp_I_Phi);
        } else {
            ips4o::sort(&I_Phi[0],&I_Phi[r],comp_I_Phi);
        }

        if (log) {
            if (mf_idx.is_open()) mf_idx << " time_sort_iphi=" << time_diff_ns(time,now());
            time = log_runtime(time);
        }
    }

    if (p > 1 && 1000*p > r) {
        p = std::max((uint64_t)1,r/1000);
        if (log) std::cout << "warning: p > r/1000, setting p to r/1000 ~ " << std::to_string(p) << std::endl;
    }

    if (p_d > 1 && 1000*p_d > r) {
        p_d = std::max((uint64_t)1,r/1000);
        if (log) std::cout << "warning: p_d > r/1000, setting p_d to r/1000 ~ " << std::to_string(p_d) << std::endl;
    }

    index_file.write((char*)&support_uint,sizeof(uint8_t));
    index_file.write((char*)&bytes_pos,sizeof(uint8_t));
    std::streampos pos_bytes_idx = index_file.tellp();
    index_file.write((char*)&bytes_idx,sizeof(uint8_t));
    index_file.write((char*)&bytes_offs,sizeof(uint8_t));
    index_file.write((char*)&p_r,sizeof(uint16_t));
    index_file.write((char*)&p_d,sizeof(uint16_t));
    index_file.write((char*)&n,sizeof(uint64_t));
    index_file.write((char*)&sigma,sizeof(uint8_t));
    index_file.write((char*)&a,sizeof(uint16_t));
    index_file.write((char*)&r,sizeof(uint64_t));
    index_file.write((char*)&chars_mapped,sizeof(bool));
    if (chars_mapped) index_file.write((char*)&map_char[0],256*sizeof(uint8_t));
    index_file.write((char*)&D_e[0],p_r*sizeof(uint64_t));
    D_e.clear();

    if (log && mf_mds.is_open()) {
        mf_mds << "RESULT"
               << " type=build_mlf"
               << " text=" << name_textfile
               << " p=" << p
               << " a=" << a;
    }

    bool store_on_disk;
    
    switch (memory_usage_mode) {
        case low: store_on_disk = true; break;
        case high: store_on_disk = false; break;
        case automatic: store_on_disk = malloc_count_current()+(r*(a/(float)(a-1))*(1+epsilon)+3*p)*(2*sizeof(uint_t)+sizeof(pair_tree_node<uint_t>)) > ram_size();
    }

    if (store_on_disk) {
        if (log) std::cout << "storing L" << (support == revert_count_locate ? ", SA and I_Phi" : "") << " on the disk to reduce ram usage" << std::flush;

        if (support == revert_count_locate) {
            std::ofstream SA_file("tmp_SA");
            write_to_file(SA_file,(uint8_t*)&SA[0],n*sizeof(int_t));
            SA_file.close();
            free(SA);
            SA = NULL;

            std::ofstream I_Phi_file("tmp_I_Phi");
            write_to_file(I_Phi_file,(uint8_t*)&I_Phi[0],(r+1)*sizeof(std::pair<uint_t,uint_t>));
            I_Phi_file.close();
            free(I_Phi);
            I_Phi = NULL;
        }
                    
        std::ofstream L_file("tmp_L");
        write_to_file(L_file,&L[0],n);
        L_file.close();
        free(L);
        L = NULL;

        if (log) {
            time_store_on_disk += time_diff_ns(time,now());
            time = log_runtime(time);
        }
    }

    if (log) {
        std::cout << std::endl << "building M_LF:" << std::flush;
    }

    // Build the move data structure M_LF out of I_LF.
    uint_t* D_p_LF = NULL;
    std::pair<uint_t,uint_t>* D_offsidx_LF = NULL;
    uint64_t r_;
    mds_builder<uint_t>::build_mds(I_LF,n,r,D_p_LF,D_offsidx_LF,r_,pow(2,8*bytes_offs)-1,a,p,log,mf_mds);
    index_file.write((char*)&r_,sizeof(uint64_t));

    if (log) {
        if (mf_mds.is_open()) mf_mds << std::endl;
        if (mf_idx.is_open()) {
            mf_idx << " time_build_mlf=" << time_diff_ns(time,now())
                    << " r_=" << r_;
        }
        time = log_runtime(time);
        std::cout << std::endl << "serializing M_LF:" << std::flush;
    }

    for (uint8_t bytes=1; bytes<=5; bytes++) {
        if (r_ <= pow(2,8*bytes)-1) {
            bytes_idx = bytes;
            break;
        }
    }

    size_index_in_ram += 67+(r_+1)*(bytes_pos+bytes_idx+bytes_offs);

    if (mf_idx.is_open()) {
        mf_idx << " omega_p=" << std::to_string(bytes_pos*8);
        mf_idx << " omega_idx=" << std::to_string(bytes_idx*8);
        mf_idx << " omega_offs=" << std::to_string(bytes_offs*8);
    }

    serialize_mds<uint_t>(
        n,
        r_,
        D_p_LF,
        D_offsidx_LF,
        bytes_pos,
        bytes_idx,
        bytes_offs,
        index_file,
        time_write_index_file,
        time_encode,
        p,
        p_d,
        log,
        mf_idx,
        mf_mds
    );

    free(D_offsidx_LF);
    D_offsidx_LF = NULL;

    if (log) {
        time = log_runtime(time);
        std::cout << std::endl;
    }

    if (store_on_disk) {
        if (log) std::cout << "reading L" << (support == revert_count_locate ? " and SA" : "") << " from disk" << std::flush;

        if (support == revert_count_locate) {
            SA = (int_t*) malloc(n*sizeof(int_t));
            std::ifstream SA_file("tmp_SA");
            read_from_file(SA_file,(uint8_t*)&SA[0],n*sizeof(int_t));
            SA_file.close();
            std::remove("tmp_SA");
        }

        L = (uint8_t*) malloc(n);
        std::ifstream L_file("tmp_L");
        read_from_file(L_file,&L[0],n);
        L_file.close();
        std::remove("tmp_L");

        if (log) {
            time_store_on_disk += time_diff_ns(time,now());
            time = log_runtime(time);
        }
    }

    uint_t* SA_s = NULL;

    if (support == revert_count_locate) {
        if (log) std::cout << "building SA_s" << std::flush;

        SA_s = (uint_t*) malloc(r_*sizeof(uint_t));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r_; i++) {
            SA_s[i] = SA[D_p_LF[i+1]-1];
        }

        free(SA);
        SA = NULL;

        if (log) {
            if (mf_idx.is_open()) mf_idx << " time_build_sas=" << time_diff_ns(time,now());
            time = log_runtime(time);
        }

        if (store_on_disk) {
            if (log) std::cout << "storing SA_s on the disk to reduce memory usage" << std::flush;
            std::ofstream SA_s_file("tmp_SA_s");
            write_to_file(SA_s_file,(uint8_t*)&SA_s[0],r_*sizeof(uint_t));
            SA_s_file.close();
            free(SA_s);
            SA_s = NULL;

            if (log) {
                time_store_on_disk += time_diff_ns(time,now());
                time = log_runtime(time);
            }
        }
    }

    if (log) std::cout << "building L'" << std::flush;
    uint8_t* L_ = (uint8_t*) malloc(r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<r_; i++) {
        L_[i] = L[D_p_LF[i]];
    }

    size_index_in_ram += r_+1;
    free(D_p_LF);
    D_p_LF = NULL;
    free(L);
    L = NULL;

    if (log) {
        if (mf_idx.is_open()) mf_idx << " time_build_l_=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "huffman-encoding L'" << std::flush;
    }

    par_huff_string L_phs;
    L_phs.build([&L_](uint64_t i){return L_[i];},0,r_-1,p_d,p);

    if (log) {
        time_encode += time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "compression ratio: " << L_phs.get_compression_ratio();
        std::cout << ", " << format_size(r_) << " -> " << format_size(L_phs.get_bit_size()/8);
        std::cout << std::endl << "serializing L'" << std::flush;
    }

    L_phs.serialize(index_file);
    L_phs.clear();
    if (log) time_write_index_file += time_diff_ns(time,now());
    uint64_t size_rs_l_ = 0;

    if (support != revert) {
        if (log) {
            time = log_runtime(time);
            std::cout << "building RS_L'" << std::flush;
        }
        
        std::vector<uint8_t> chars;
        chars.emplace_back(1);
        
        for (uint16_t i=0; i<256; i++) {
            if (contains_char[i]) {
                chars.emplace_back(chars_mapped ? map_char[i] : i);
            }
        }

        contains_char.clear();
        if (chars_mapped) map_char.clear();
        index_file.write((char*)&chars[0],sigma);
        std::vector<sdsl::bit_vector> RS_L_plain_bvs(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<sigma; i++) {
            RS_L_plain_bvs[chars[i]] = sdsl::bit_vector(r_);
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r_; i++) {
            RS_L_plain_bvs[L_[i]][i] = 1;
        }

        free(L_);
        L_ = NULL;
        std::vector<sdsl::sd_vector<>> RS_L_bvs(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<sigma; i++) {
            RS_L_bvs[chars[i]] = sdsl::sd_vector<>(RS_L_plain_bvs[chars[i]]);
        }

        RS_L_plain_bvs.clear();

        std::vector<sdsl::sd_vector<>::rank_1_type> RS_L_rank_1(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<sigma; i++) {
            RS_L_rank_1[chars[i]] = sdsl::sd_vector<>::rank_1_type(&RS_L_bvs[chars[i]]);
        }

        std::vector<sdsl::select_support_sd<>> RS_L_select_1(256);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i=0; i<sigma; i++) {
            RS_L_select_1[chars[i]] = sdsl::sd_vector<>::select_1_type(&RS_L_bvs[chars[i]]);
        }

        size_rs_l_ = 0;

        for (uint16_t i=0; i<256; i++) {
            size_rs_l_ += sdsl::size_in_bytes(RS_L_bvs[i]);
            size_rs_l_ += sdsl::size_in_bytes(RS_L_rank_1[i]);
            size_rs_l_ += sdsl::size_in_bytes(RS_L_select_1[i]);
        }

        size_index_in_ram += size_rs_l_;

        if (log) {
            std::cout << " (" << format_size(size_rs_l_) << ")";
            if (mf_idx.is_open()) mf_idx << " time_build_rs_l_=" << time_diff_ns(time,now());
            time = log_runtime(time);
            std::cout << "serializing RS_L'"<< std::flush;
        }

        for (uint16_t i=0; i<sigma; i++) {
            RS_L_bvs[chars[i]].serialize(index_file);
        }

        RS_L_bvs.clear();
        
        for (uint16_t i=0; i<sigma; i++) {
            RS_L_rank_1[chars[i]].serialize(index_file);
        }

        RS_L_rank_1.clear();

        for (uint16_t i=0; i<sigma; i++) {
            RS_L_select_1[chars[i]].serialize(index_file);
        }

        RS_L_select_1.clear();
        chars.clear();

        if (log) time_write_index_file += time_diff_ns(time,now());
    }

    if (log) time = log_runtime(time);

    uint64_t r__;
    uint_t* D_p_Phi = NULL;
    std::pair<uint_t,uint_t>* D_offsidx_Phi = NULL;

    if (support == revert_count_locate) {
        if (store_on_disk) {
            if (log) std::cout << "reading I_Phi from disk" << std::flush;
            I_Phi = (std::pair<uint_t,uint_t>*) malloc((r+1)*sizeof(std::pair<uint_t,uint_t>));
            std::ifstream I_Phi_file("tmp_I_Phi");
            read_from_file(I_Phi_file,(uint8_t*)&I_Phi[0],(r+1)*sizeof(std::pair<uint_t,uint_t>));
            I_Phi_file.close();
            std::remove("tmp_I_Phi");

            if (log) {
                time_store_on_disk += time_diff_ns(time,now());
                time = log_runtime(time);
            }
        }

        if (log) {
            if (mf_mds.is_open()) {
                mf_mds << "RESULT"
                        << " type=build_mphi"
                        << " text=" << name_textfile
                        << " p=" << p
                        << " a=" << a;
            }
            std::cout << std::endl << "building M_Phi:" << std::flush;
        }

        mds_builder<uint_t>::build_mds(I_Phi,n,r,D_p_Phi,D_offsidx_Phi,r__,pow(2,8*bytes_offs)-1,a,p,log,mf_mds);
        index_file.write((char*)&r__,sizeof(uint64_t));
        size_index_in_ram += (r__+1)*(bytes_pos+bytes_idx+bytes_offs);

        if (log) {
            if (mf_mds.is_open()) mf_mds << std::endl;
            if (mf_idx.is_open()) {
                mf_idx << " time_build_mphi=" << time_diff_ns(time,now());
                mf_idx << " r__=" << r__;
            }
            time = log_runtime(time);
            std::cout << std::endl;
        }

        for (uint8_t bytes=1; bytes<=5; bytes++) {
            if (std::max(r_,r__) <= pow(2,8*bytes)-1) {
                bytes_idx = bytes;
                break;
            }
        }
    } else {
        for (uint8_t bytes=1; bytes<=5; bytes++) {
            if (r_ <= pow(2,8*bytes)-1) {
                bytes_idx = bytes;
                break;
            }
        }
    }

    std::streampos int_tmp = index_file.tellp();
    index_file.seekp(pos_bytes_idx);
    index_file.write((char*)&bytes_idx,sizeof(uint8_t));

    if (support == revert_count_locate) {
        index_file.seekp(int_tmp);
        if (log) std::cout << "serializing M_Phi:" << std::flush;

        serialize_mds<uint_t>(
            n,
            r__,
            D_p_Phi,
            D_offsidx_Phi,
            bytes_pos,
            bytes_idx,
            bytes_offs,
            index_file,
            time_write_index_file,
            time_encode,
            p,
            p_d,
            log,
            mf_idx,
            mf_mds
        );

        free(D_offsidx_Phi);
        D_offsidx_Phi = NULL;

        if (log) {
            time = log_runtime(time);
            std::cout << std::endl;
        }

        if (store_on_disk) {
            if (log) std::cout << "reading SA_s from disk" << std::flush;
            SA_s = (uint_t*) malloc(r_*sizeof(uint_t));
            std::ifstream SA_s_file("tmp_SA_s");
            read_from_file(SA_s_file,(uint8_t*)&SA_s[0],r_*sizeof(uint_t));
            SA_s_file.close();
            std::remove("tmp_SA_s");

            if (log) {
                time_store_on_disk += time_diff_ns(time,now());
                time = log_runtime(time);
            }
        }

        if (log) std::cout << "building SA_offs and SA_idx" << std::flush;
        uint_t* pi_ = (uint_t*) malloc(r_*sizeof(uint_t));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r_; i++) {
            pi_[i] = i;
        }

        auto comp_pi_ = [&SA_s](auto i, auto j) {return SA_s[i] < SA_s[j];};
        if (p > 1) {
            ips4o::parallel::sort(&pi_[0],&pi_[r_],comp_pi_);
        } else {
            ips4o::sort(&pi_[0],&pi_[r_],comp_pi_);
        }

        std::pair<uint_t,uint_t>* SA_offs_idx = (std::pair<uint_t,uint_t>*) malloc(r_*sizeof(std::pair<uint_t,uint_t>));

        /* Now we will divide the range [0..n-1] up into p non-overlapping sub-ranges [s[i_p]..s[i_p+1]-1],
        for each i_p \in [0..p-1], with 0 = s[0] < s[1] < ... < s[p] = n, where
        s[i_p] = min {s' \in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * \lfloor (r'+r'')/p \rfloor, where 
                        x[i_p] = min {x' \in [0,r''-1], s.t. M_Phi.p(x') >= s'} and
                        u[i_p] = min {u' \in [0,r'-1], s.t. SA_s[u'] >= s'}
        }.
        By doing so, we ensure, that the number of the input intervals of M_Phi starting in the range
        [s[i_p]..s[i_p+1]-1] plus the number of suffix array samples in SA_s lying in the range
        [s[i_p]..s[i_p+1]-1] is \lfloor (r'+r'')/p \rfloor +- 1. This property is useful, because it
        ensures that if with each thread i_p, we simultaneously iterate over those, then each thread
        iterates over almost exactly the same number \lfloor (r'+r'')/p \rfloor +- 1 of entries in M_Phi
        and SA_s combined. This way, we can acheive good load-balancing. Because we do not have to access
        s[0..p] later, we will not store those values in an array. */

        /* [0..p], x[i_p] = min {x' \in [0,r''-1], s.t. M_Phi.p(x') >= s'} stores the number of input
        intervals in M_Phi starting before s[i_p]. */
        std::vector<uint_t> x(p+1);
        x[0] = 0;
        x[p] = r__;

        /* [0..p], u[i_p] = min {u' \in [0,r'-1], s.t. SA_s[u'] >= s'} stores the number of suffix array
        samples in SA_s that are smaller than s[i_p]. */
        std::vector<uint_t> u(p+1);
        u[0] = 0;
        u[p] = r_;

        // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
        #pragma omp parallel num_threads(p)
        {
            // Index \in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            // The optimal value i_p * \lfloor (r'+r'')/p \rfloor for s[i_p].
            uint_t o = i_p*((r_+r__)/p);

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
                r_x = r__-1;
                while (l_x != r_x) {
                    m_x = l_x+(r_x-l_x)/2;
                    if (D_p_Phi[m_x] < m_s) {
                        l_x = m_x+1;
                    } else {
                        r_x = m_x;
                    }
                }

                // Find the minimum u' \in [0,r'-1], s.t. SA_s[pi'[u']] >= m_s.
                l_u = 0;
                r_u = r_-1;
                while (l_u != r_u) {
                    m_u = l_u+(r_u-l_u)/2;
                    if (SA_s[pi_[m_u]] < m_s) {
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
                if (SA_s[pi_[j]] < D_p_Phi[i]) {
                    i--;
                }

                /* Iterate until one of the iteration end positions i_ and j_ has
                been reached. */
                while (i <= i_ && j <= j_) {

                    /* Iterate over the suffix array samples that lie in the current
                    i-th input interval of M_Phi. */
                    while (j <= j_ && SA_s[pi_[j]] < D_p_Phi[i+1]) {
                        
                        /* Because each of those j-th largest suffix array samples lie in the i-th
                        input interval of M_Phi, we can set SA_idx[\pi'[j]] = i for each of them. */
                        SA_offs_idx[pi_[j]].first = SA_s[pi_[j]]-D_p_Phi[i];
                        SA_offs_idx[pi_[j]].second = i;

                        j++;
                    }

                    i++;
                };
            }
        }

        size_index_in_ram += r_*bytes_idx;
        size_index_in_ram += r_*bytes_offs;

        free(SA_s);
        SA_s = NULL;
        free(D_p_Phi);
        D_p_Phi = NULL;
        free(pi_);
        pi_ = NULL;
        x.clear();
        u.clear();

        if (log) {
            if (mf_idx.is_open()) mf_idx << " time_build_sa_idx_offs=" << time_diff_ns(time,now());
            time = log_runtime(time);
            std::cout << "delta-encoding SA_offs" << std::flush;
        }

        par_delta_vector<uint_t> SA_offs_pdv;
        SA_offs_pdv.build([&SA_offs_idx](uint64_t i){return SA_offs_idx[i].first;},0,r_-1,p_d,p);
        float comp_ratio = (r_*bytes_offs)/(float)(SA_offs_pdv.get_bit_size()/8);
        bool var_width_encode = comp_ratio < 1;
        index_file << var_width_encode;

        if (log) {
            time_encode += time_diff_ns(time,now());
            time = log_runtime(time);
            std::cout << "compression ratio: " << comp_ratio;
            if (var_width_encode) {
                std::cout << " < 1" << std::endl << "=> variable-width-encoding D_offs";
            } else {
                std::cout << ", " << format_size(r_*bytes_offs) << " -> " << format_size(SA_offs_pdv.get_bit_size()/8);
            }
        }

        par_var_width_vector<uint_t> SA_offs_pvwv;
        if (var_width_encode) {
            SA_offs_pdv.clear();
            SA_offs_pvwv.build([&SA_offs_idx](uint64_t i){return SA_offs_idx[i].first;},0,r_-1,bytes_offs*8,p_d,p);
            comp_ratio = SA_offs_pvwv.get_compression_ratio();
        }

        if (log) {
            time_encode += time_diff_ns(time,now());
            if (var_width_encode) {
                time = log_runtime(time);
                std::cout << "compression ratio: " << comp_ratio;
                std::cout << ", " << format_size(r_*bytes_offs) << " -> " << format_size(SA_offs_pvwv.get_bit_size()/8);
            }
            std::cout << std::endl << "serializing SA_offs" << std::flush;
        }

        if (var_width_encode) {
            SA_offs_pvwv.serialize(index_file);
            SA_offs_pvwv.clear();
        } else {
            SA_offs_pdv.serialize(index_file);
            SA_offs_pdv.clear();
        }

        if (log) {
            time_write_index_file += time_diff_ns(time,now());
            time = log_runtime(time);
            std::cout << "encoding SA_idx" << std::flush;
        }

        par_var_width_vector<uint_t> SA_idx_pvwv;
        SA_idx_pvwv.build([&SA_offs_idx](uint64_t i){return SA_offs_idx[i].second;},0,r_-1,(uint8_t)std::ceil(log2(r__-1)),p_d,p);
        comp_ratio = (r_*bytes_idx)/(float)(SA_idx_pvwv.get_bit_size()/8);

        free(SA_offs_idx);
        SA_offs_idx = NULL;

        if (log) {
            time_encode += time_diff_ns(time,now());
            time = log_runtime(time);
            std::cout << "compression ratio: " << comp_ratio;
            std::cout << ", " << format_size(r_*bytes_idx) << " -> " << format_size(SA_idx_pvwv.get_bit_size()/8);
            std::cout << std::endl << "serializing SA_idx" << std::flush;
        }
        
        SA_idx_pvwv.serialize(index_file);
        SA_idx_pvwv.clear();
        
        if (log) {
            time_write_index_file += time_diff_ns(time,now());
            time = log_runtime(time);
        }
    }

    if (log) {
        std::cout << std::endl;
        std::cout << "bytes_pos = " << std::to_string(bytes_pos);
        std::cout << ", bytes_idx = " << std::to_string(bytes_idx);
        std::cout << ", bytes_offs = " << std::to_string(bytes_offs);
        std::cout << std::endl;

        uint64_t time_total = time_diff_ns(time_start,now());
        uint64_t time_build_index_data_structures = time_total - time_write_index_file - time_encode;
        uint64_t size_index_file = index_file.tellp();

        std::cout << "index file size: " << format_size(size_index_file) << std::endl;
        std::cout << "index size in ram: " << format_size(size_index_in_ram) << std::endl;
        std::cout << "M_LF: " << format_size(67+(r_+1)*(bytes_pos+bytes_idx+bytes_offs)) << std::endl;
		std::cout << "L': " << format_size(r_+1) << std::endl;
		
		if (support != support_mode::revert) {
			std::cout << "RS_L': " << format_size(size_rs_l_) << std::endl;

			if (support == revert_count_locate) {
				std::cout << "M_Phi: " << format_size(58+(r__+1)*(bytes_pos+bytes_idx+bytes_offs)) << std::endl;
				std::cout << "SA_idx: " << format_size(r_*bytes_idx) << std::endl;
				std::cout << "SA_offs: " << format_size(r_*bytes_offs) << std::endl;
			}
		}

		std::cout << std::endl;
        if (store_on_disk) std::cout << "overall time to temporarily store index data structures on the disk: " << format_time(time_store_on_disk) << std::endl;
        std::cout << "overall time to build the index data structures: " << format_time(time_build_index_data_structures) << std::endl;
        std::cout << "overall time to encode the index data structures: " << format_time(time_encode) << std::endl;
        std::cout << "overall time to write the index file: " << format_time(time_write_index_file) << std::endl;
        std::cout << "r-index built";

        if (mf_idx.is_open()) {
            mf_idx << " time_total=" << time_total;
            mf_idx << " time_build_index_data_structures=" << time_build_index_data_structures;
            mf_idx << " time_encode=" << time_encode;
            mf_idx << " time_write_index_file=" << time_write_index_file;
            mf_idx << " peak_memory_usage=" << malloc_count_peak();
            mf_idx << " size_index_file=" << size_index_file;
            mf_idx << " size_index_in_ram=" << size_index_in_ram;

            mf_idx << " size_m_lf=" << 67+(r_+1)*(bytes_pos+bytes_idx+bytes_offs);
            mf_idx << " size_l_=" << r_+1;
            
            if (support != support_mode::revert) {
                mf_idx << " size_rs_l_=" << size_rs_l_;

                if (support == revert_count_locate) {
                    mf_idx << " size_m_phi=" << 58+(r__+1)*(bytes_pos+bytes_idx+bytes_offs);
                    mf_idx << " size_sa_idx=" << r_*bytes_idx;
                    mf_idx << " size_sa_offs=" << r_*bytes_offs;
                }
            }

            mf_idx << std::endl;
        }

        std::cout << ", peak additional memory allocation during build: ~ " << format_size(malloc_count_peak());
        std::cout << ", in " << format_time(time_total) << std::endl;
    }
}