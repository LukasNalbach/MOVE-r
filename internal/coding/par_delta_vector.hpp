#pragma once

#include <omp.h>
#include "delta_vector.hpp"

template <typename uint_t>
class par_delta_vector {
    protected:
    uint16_t num_delta_vectors = 0;
    uint64_t size = 0;
    uint64_t bit_size = 0;
    std::vector<delta_vector<uint_t>> delta_vectors;

    public:
    par_delta_vector() = default;
    par_delta_vector(const par_delta_vector&) = default;
    par_delta_vector(std::vector<uint_t> &vector, uint64_t l, uint64_t r, uint16_t num_delta_vectors = 256, uint16_t p = omp_get_max_threads()) {build([&vector](uint64_t i){return vector[i];},l,r,num_delta_vectors,p);}
    par_delta_vector(std::vector<uint_t> &vector, uint16_t num_delta_vectors = 256, uint16_t p = omp_get_max_threads()) : par_delta_vector(vector,0,vector.size()-1,num_delta_vectors,p) {}

    void build(std::function<uint_t(uint64_t)> read, uint64_t l, uint64_t r, uint16_t num_delta_vectors = 256, uint16_t p = omp_get_max_threads()) {
        size = r-l+1;
        num_delta_vectors = std::max((uint64_t)1,std::min((uint64_t)num_delta_vectors,size/1000));
        this->num_delta_vectors = num_delta_vectors;
        delta_vectors.resize(num_delta_vectors);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_delta_vectors; i++) {
			uint64_t b = l + i * (size / num_delta_vectors);
            uint64_t e = l + (i == num_delta_vectors - 1 ? size - 1 : (i + 1)*(size / num_delta_vectors ) - 1);
            delta_vectors[i].build(read,b,e);
        }

        bit_size = 18;

        for (uint16_t i = 0; i < num_delta_vectors; i++) {
            bit_size += delta_vectors[i].get_bit_size();
        }
    }

    void decode(std::function<void(uint64_t,uint_t)> write, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_delta_vectors; i++) {
			uint64_t b = l + i * (size / num_delta_vectors);
            uint64_t e = l + (i == num_delta_vectors - 1 ? size - 1 : (i + 1) * (size / num_delta_vectors) - 1);
            delta_vectors[i].decode(write,b,e);
        }
    }

    void decode(std::vector<uint_t> &vector, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        decode([&vector](uint64_t i, uint_t c){vector[i] = c;},l,r,p);
    }

    void decode(std::vector<uint_t> &vector, uint16_t p = omp_get_max_threads()) {
        decode(vector,0,vector.size()-1,p);
    }

    std::vector<uint_t> decode(uint16_t p = omp_get_max_threads()) {
        std::vector<uint_t> vector(size);
        decode(vector,p);
        return vector;
    }

    void load_and_decode_buffered(std::ifstream &in, std::function<void(uint64_t,uint_t)> write, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        in.read((char*)&num_delta_vectors,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        delta_vectors.resize(num_delta_vectors);
        uint16_t num_delta_vectors_left = num_delta_vectors;
        uint16_t num_delta_vectors_to_read;

        while (num_delta_vectors_left > 0) {
            num_delta_vectors_to_read = std::min(num_delta_vectors_left,(uint16_t)p);

            for (uint16_t i = num_delta_vectors - num_delta_vectors_left; i < num_delta_vectors - num_delta_vectors_left + num_delta_vectors_to_read; i++) {
                delta_vectors[i].load(in);
            }

            #pragma omp parallel for num_threads(p)
            for (uint16_t i = num_delta_vectors - num_delta_vectors_left; i < num_delta_vectors - num_delta_vectors_left + num_delta_vectors_to_read; i++) {
                uint64_t b = l + i * (size / num_delta_vectors);
                uint64_t e = l + (i == num_delta_vectors - 1 ? size - 1 : (i + 1) * (size / num_delta_vectors) - 1);

                delta_vectors[i].decode(write,b,e);
                delta_vectors[i].clear();
            }

            num_delta_vectors_left -= num_delta_vectors_to_read;
        }

        clear();
    }

    void clear() {
        num_delta_vectors = 0;
        size = 0;
        bit_size = 0;
        delta_vectors.clear();
    }

    uint64_t get_size() {
        return size;
    }

    uint64_t get_bit_size() {
        return bit_size;
    }

    float get_compression_ratio() {
        return (size*sizeof(uint_t)*8)/(float)bit_size;
    }

    void load(std::ifstream &in) {
        in.read((char*)&num_delta_vectors,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        delta_vectors.resize(num_delta_vectors);

        for (uint16_t i=0; i<num_delta_vectors; i++) {
            delta_vectors[i].load(in);
        }
    }

    void serialize(std::ofstream &out) {
        out.write((char*)&num_delta_vectors,sizeof(uint16_t));
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&bit_size,sizeof(uint64_t));

        for (uint16_t i=0; i<num_delta_vectors; i++) {
            delta_vectors[i].serialize(out);
        }
    }
};