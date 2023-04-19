#pragma once

#include <omp.h>
#include "huff_string.hpp"

class par_huff_string {
    protected:
    uint16_t num_huff_strings = 0;
    uint64_t size = 0;
    uint64_t bit_size = 0;
    std::vector<huff_string> huff_strings;

    public:
    par_huff_string() = default;
    par_huff_string(const par_huff_string&) = default;
    par_huff_string(std::vector<uint8_t> &string, uint64_t l, uint64_t r, uint16_t num_huff_strings = 256, uint16_t p = omp_get_max_threads()) {build([&string](uint64_t i){return string[i];},l,r,num_huff_strings,p);}
    par_huff_string(std::vector<uint8_t> &string, uint16_t num_huff_strings = 256, uint16_t p = omp_get_max_threads()) : par_huff_string(string,0,string.size()-1,num_huff_strings,p) {}

    void build(std::function<uint8_t(uint64_t)> read, uint64_t l, uint64_t r, uint16_t num_huff_strings = 256, uint16_t p = omp_get_max_threads()) {
        size = r-l+1;
        num_huff_strings = std::max((uint64_t)1,std::min((uint64_t)num_huff_strings,size/1000));
        this->num_huff_strings = num_huff_strings;
        huff_strings.resize(num_huff_strings);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_huff_strings; i++) {
			uint64_t b = l + i * (size / num_huff_strings);
            uint64_t e = l + (i == num_huff_strings - 1 ? size - 1 : (i + 1) * (size / num_huff_strings) - 1);
            huff_strings[i].build(read,b,e);
        }

        bit_size = 18;

        for (uint16_t i = 0; i < num_huff_strings; i++) {
            bit_size += huff_strings[i].get_bit_size();
        }
    }

    void decode(std::function<void(uint64_t,uint8_t)> write, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_huff_strings; i++) {
			uint64_t b = l + i * (size / num_huff_strings);
            uint64_t e = l + (i == num_huff_strings - 1 ? size - 1 : (i + 1) * (size / num_huff_strings) - 1);
            huff_strings[i].decode(write,b,e);
        }
    }

    void decode(std::vector<uint8_t> &string, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        decode([&string](uint64_t i, uint8_t c){string[i] = c;},l,r,p);
    }

    void decode(std::vector<uint8_t> &string, uint16_t p = omp_get_max_threads()) {
        decode(string,0,string.size()-1,p);
    }

    std::vector<uint8_t> decode(uint16_t p = omp_get_max_threads()) {
        std::vector<uint8_t> string(size);
        decode(string,p);
        return string;
    }

    void load_and_decode_buffered(std::ifstream &in, std::function<void(uint64_t,uint8_t)> write, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        in.read((char*)&num_huff_strings,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        huff_strings.resize(num_huff_strings);
        uint16_t num_huff_strings_left = num_huff_strings;
        uint16_t num_huff_strings_to_read;

        while (num_huff_strings_left > 0) {
            num_huff_strings_to_read = std::min(num_huff_strings_left,(uint16_t)p);

            for (uint16_t i = num_huff_strings - num_huff_strings_left; i < num_huff_strings - num_huff_strings_left + num_huff_strings_to_read; i++) {
                huff_strings[i].load(in);
            }

            #pragma omp parallel for num_threads(p)
            for (uint16_t i = num_huff_strings - num_huff_strings_left; i < num_huff_strings - num_huff_strings_left + num_huff_strings_to_read; i++) {
                uint64_t b = l + i * (size / num_huff_strings);
                uint64_t e = l + (i == num_huff_strings - 1 ? size - 1 : (i + 1) * (size / num_huff_strings) - 1);

                huff_strings[i].decode(write,b,e);
                huff_strings[i].clear();
            }

            num_huff_strings_left -= num_huff_strings_to_read;
        }

        clear();
    }

    void clear() {
        num_huff_strings = 0;
        size = 0;
        bit_size = 0;
        huff_strings.clear();
    }

    uint64_t get_size() {
        return size;
    }

    uint64_t get_bit_size() {
        return bit_size;
    }

    float get_compression_ratio() {
        return (size*8)/(float)bit_size;
    }

    void load(std::ifstream &in) {
        in.read((char*)&num_huff_strings,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        huff_strings.resize(num_huff_strings);

        for (uint16_t i=0; i<num_huff_strings; i++) {
            huff_strings[i].load(in);
        }
    }

    void serialize(std::ofstream &out) {
        out.write((char*)&num_huff_strings,sizeof(uint16_t));
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&bit_size,sizeof(uint64_t));

        for (uint16_t i=0; i<num_huff_strings; i++) {
            huff_strings[i].serialize(out);
        }
    }
};