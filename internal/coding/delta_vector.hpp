#pragma once

#include <misc/utils.hpp>

template <typename uint_t>
class delta_vector {
    protected:
    uint64_t size = 0;
    uint64_t delta_data_bit_size = 0;
    std::vector<uint8_t,default_init_allocator<uint8_t>> delta_data;

    inline static std::pair<uint_t,uint64_t> encode_delta(uint_t value) {
        if (value == 0) return std::pair<uint_t,uint64_t>{1, (const uint64_t) (((uint64_t) 1) << (64 - 1))};
        uint_t log2_value = log2_clz<uint64_t>(value);
        uint_t _2_log2_log2_value_p1_p1 = 2 * (log2_clz<uint_t>(log2_value + 1) + 1);
        return std::pair<uint_t,uint64_t>{
            _2_log2_log2_value_p1_p1 + log2_value,
            ((uint64_t) (log2_value + 1)) << (64 - _2_log2_log2_value_p1_p1) | (((uint64_t) (value & ~(((uint_t) 1) << log2_value))) << (64 - log2_value - _2_log2_log2_value_p1_p1))
        };
    }

    inline static std::pair<uint_t,uint_t> decode_delta(uint64_t delta) {
        uint_t _2_leading_zeros_delta = 2 * __builtin_clzl(delta);
        if (_2_leading_zeros_delta == 0) return std::pair<uint_t,uint_t>{1,0};
        uint_t gamma_m1 = ((uint_t) (delta >> (64 - _2_leading_zeros_delta))) - 1;
        uint_t code_length_delta = _2_leading_zeros_delta + gamma_m1;
        return std::pair<uint_t,uint_t>{
            code_length_delta,
            ((uint_t) ((((uint64_t) -1) >> (64 - gamma_m1)) & (delta >> (64 - code_length_delta)))) | (((uint_t) 1) << gamma_m1)
        };
    }

    public:
    delta_vector() = default;
    delta_vector(const delta_vector&) = default;
    delta_vector(std::vector<uint_t> &vector, uint64_t l, uint64_t r) {build([&vector](uint64_t i){return vector[i];},l,r);}
    delta_vector(std::vector<uint_t> &vector) : delta_vector(vector,0,vector.size()-1) {}

    void build(std::function<uint_t(uint64_t)> read, uint64_t l, uint64_t r) {
        size = r-l+1;
        delta_data_bit_size = 0.75*(r-l+1)*sizeof(uint_t);
        delta_data.resize(delta_data_bit_size/8+1);
        std::memset(&delta_data[0],0,delta_data.size());
        std::vector<std::pair<uint_t,uint64_t>> precomputed_codes;
        uint64_t block_size = std::min((uint64_t) 256, r-l+1);
        precomputed_codes.reserve(block_size);
        std::pair<uint_t,uint64_t> delta_pair;
        uint64_t idx_in_delta_data = 0;
        uint64_t offs_in_delta_data = 0;
        uint64_t offs_in_code = 0;
        uint64_t written_bits;
        uint64_t current_block_length;
        uint64_t block_start_pos_in_vector = l;
        uint64_t block_delta_bit_size;

        while (block_start_pos_in_vector <= r) {
            current_block_length = std::min(block_size, r - block_start_pos_in_vector + 1);
            block_delta_bit_size = 0;
            
            for (uint64_t i=0; i<current_block_length; i++) {
                delta_pair = encode_delta(read(block_start_pos_in_vector + i));
                block_delta_bit_size += delta_pair.first;
                precomputed_codes.emplace_back(delta_pair);
            }

            delta_data_bit_size += block_delta_bit_size;

            while (delta_data.size()*8 < delta_data_bit_size) {
                std::vector<uint8_t,default_init_allocator<uint8_t>> tmp(delta_data.size()*2);
                std::copy(delta_data.begin(),delta_data.end(),tmp.begin());
                delta_data.swap(tmp);
                std::memset(&delta_data[delta_data.size()/2],0,delta_data.size()/2);
            }

            for (uint64_t i=0; i<current_block_length; i++) {
                delta_pair = precomputed_codes[i];

                while (offs_in_code < delta_pair.first) {
                    delta_data[idx_in_delta_data] |= (uint8_t) (delta_pair.second >> (56 + offs_in_delta_data - offs_in_code));
                    written_bits = std::min((uint64_t) 8 - offs_in_delta_data, (uint64_t) delta_pair.first - offs_in_code);
                    offs_in_code += written_bits;
                    offs_in_delta_data += written_bits;

                    if (offs_in_delta_data >= 8) {
                        idx_in_delta_data++;
                        offs_in_delta_data -= 8;
                    }
                }

                offs_in_code = 0;
            }
            
            precomputed_codes.clear();
            block_start_pos_in_vector += current_block_length;
        }
    }

    void decode(std::function<void(uint64_t,uint_t)> write, uint64_t l, uint64_t r) {
        std::pair<uint_t,uint_t> delta_pair;
        uint64_t idx_in_delta_data = 0;
        uint64_t offs_in_delta_data = 0;
        uint64_t code = 0;
        uint64_t offs_in_code = 0;
        int8_t shift_left_delta_data;

        for (uint64_t i=l; i<=r; i++) {
            while (offs_in_code < 64) {
                shift_left_delta_data = 56 + offs_in_delta_data - offs_in_code;

                if (shift_left_delta_data >= 0) {
                    code |= ((uint64_t) delta_data[idx_in_delta_data]) << shift_left_delta_data;
                    offs_in_code += 8 - offs_in_delta_data;
                    idx_in_delta_data++;
                    offs_in_delta_data = 0;
                } else {
                    code |= ((uint64_t) delta_data[idx_in_delta_data]) >> - shift_left_delta_data;
                    offs_in_delta_data += 64 - offs_in_code;
                    offs_in_code = 64;
                }
            }

            delta_pair = decode_delta(code);
            write(i, delta_pair.second);
            code <<= delta_pair.first;
            offs_in_code -= delta_pair.first;
        }
    }

    void decode(std::vector<uint_t> &vector, uint64_t l, uint64_t r) {
        decode([&vector](uint64_t i, uint_t val){vector[i] = val;},l,r);
    }

    void decode(std::vector<uint_t> &vector) {
        decode(vector,0,vector.size()-1);
    }

    std::vector<uint_t> decode() {
        std::vector<uint_t> vector;
        vector.resize(size);
        decode(vector);
        return vector;
    }

    void clear() {
        size = 0;
        delta_data_bit_size = 0;
        delta_data.clear();
    }

    uint64_t get_size() {
        return size;
    }

    uint64_t get_bit_size() {
        return 16+delta_data_bit_size;
    }

    float get_compression_ratio() {
        return (size*sizeof(uint_t)*8)/(float)get_bit_size();
    }

    void load(std::ifstream &in) {
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&delta_data_bit_size,sizeof(uint64_t));
        delta_data.resize(delta_data_bit_size/8+1);
        in.read((char*)&delta_data[0],delta_data_bit_size/8+1);
    }

    void serialize(std::ofstream &out) {
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&delta_data_bit_size,sizeof(uint64_t));
        out.write((char*)&delta_data[0],delta_data_bit_size/8+1);
    }
};