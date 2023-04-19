/**
 * @brief reads a serialized r-index from an input stream
 * @param index_file an input stream storing a serialized r-index
 * @param p the number of threads to use (default: all threads)
 */
template <typename pos_t, typename idx_t, unsigned int bytes_pos, unsigned int bytes_idx, unsigned int bytes_offs>
void MOVE_r<pos_t,idx_t,bytes_pos,bytes_idx,bytes_offs>::load(std::ifstream& index_file, support_mode support, bool log, uint16_t p, std::ofstream &mf) {
    malloc_count_reset_peak();
    auto time = now();
    auto time_start = time;
    this->support = support;

    std::cout << std::endl;
    std::cout << "bytes_pos = " << std::to_string(bytes_pos);
    std::cout << ", bytes_idx = " << std::to_string(bytes_idx);
    std::cout << ", bytes_offs = " << std::to_string(bytes_offs);
    std::cout << std::endl;

    if (log) std::cout << std::endl << "reconstructing M_LF:" << std::flush;

    index_file.seekg(4);
    index_file.read((char*)&p_r,sizeof(uint16_t));
    index_file.read((char*)&p_d,sizeof(uint16_t));
    index_file.read((char*)&n,sizeof(uint64_t));
    index_file.read((char*)&sigma,sizeof(uint8_t));
    index_file.read((char*)&a,sizeof(uint16_t));
    index_file.read((char*)&r,sizeof(uint64_t));
    index_file.read((char*)&chars_mapped,sizeof(bool));

    if (chars_mapped) {
        map_char.resize(256);
        index_file.read((char*)&map_char[0],256*sizeof(unsigned char));
        map_char_rev.resize(256);

        for (uint16_t i=0; i<256; i++) {
            map_char_rev[map_char[i]] = i;
        }
    }

    D_e.resize(p_r);
    index_file.read((char*)&D_e[0],p_r*sizeof(uint64_t));
    index_file.read((char*)&r_,sizeof(uint64_t));
    M_LF.load(index_file,log,p);

    if (log) {
        time = log_runtime(time);
        std::cout << std::endl << "reconstructing L'" << std::flush;
    }

    par_huff_string L_phs;
    L_phs.load_and_decode_buffered(index_file,[this](uint64_t i, uint8_t c){M_LF.c(i) = c;},0,r-1,p);

    if (support != support_mode::revert) {
        if (log) {
            time = log_runtime(time);
            std::cout << "reconstructing RS_L'" << std::flush;
        }

        std::vector<uint8_t> chars(sigma,0);
        index_file.read((char*)&chars[0],sigma);

        RS_L_bvs.resize(256);

        for (uint16_t i=0; i<sigma; i++) {
            RS_L_bvs[chars[i]].load(index_file);
        }

        RS_L_rank_1.resize(256);

        for (uint16_t i=0; i<sigma; i++) {
            RS_L_rank_1[chars[i]].load(index_file);
            RS_L_rank_1[chars[i]].set_vector(&RS_L_bvs[chars[i]]);
        }

        RS_L_select_1.resize(256);

        for (uint16_t i=0; i<sigma; i++) {
            RS_L_select_1[chars[i]].load(index_file);
            RS_L_select_1[chars[i]].set_vector(&RS_L_bvs[chars[i]]);
        }
    }
    
    if (support == revert_count_locate) {
        if (log) {
            time = log_runtime(time);
            std::cout << std::endl << "reconstructing M_Phi:" << std::flush;
        }

        index_file.read((char*)&r__,sizeof(uint64_t));
        M_Phi.load(index_file,log,p);

        if (log) {
            time = log_runtime(time);
            std::cout << std::endl << "reconstructing SA_offs" << std::flush;
        }

        SA_idxoffs = (uint8_t*) malloc(r_*size_entry_saidxoffs);
        base_saidx = SA_idxoffs;
        base_saoffs = SA_idxoffs + pos_saoffs;

        bool var_width_encode;
        index_file >> var_width_encode;
        if (var_width_encode) {
            par_var_width_vector<pos_t> SA_offs_pvwv;
            SA_offs_pvwv.load_and_decode_buffered(index_file,[this](uint64_t i, pos_t offs){set_SA_offs(i,offs);},0,r_-1,p);
        } else {
            par_delta_vector<pos_t> SA_offs_pdv;
            SA_offs_pdv.load_and_decode_buffered(index_file,[this](uint64_t i, pos_t offs){set_SA_offs(i,offs);},0,r_-1,p);
        }

        if (log) {
            time = log_runtime(time);
            std::cout << "reconstructing SA_idx" << std::flush;
        }

        par_var_width_vector<idx_t> SA_idx_pdv;
        SA_idx_pdv.load_and_decode_buffered(index_file,[this](uint64_t i, idx_t idx){set_SA_idx(i,idx);},0,r_-1,p);
    }

    if (log) {
        time = log_runtime(time);
        
        std::cout << std::endl;
        std::cout << "r-index reconstructed, peak memory consumption: " << format_size(malloc_count_peak());
        std::cout << ", in " << format_time(time_diff_ns(time_start,now())) << std::endl;
    }
}