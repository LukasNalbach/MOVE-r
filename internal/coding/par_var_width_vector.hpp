#pragma once

#include <omp.h>
#include <sdsl/vectors.hpp>

template <typename uint_t>
class par_var_width_vector {
    protected:
    uint16_t num_vectors;
    uint64_t size;
    uint64_t bit_size;
    std::vector<sdsl::int_vector<0>> vectors;

    public:
    par_var_width_vector() = default;
    par_var_width_vector(const par_var_width_vector&) = default;
    par_var_width_vector(std::vector<uint_t> &vector, uint64_t l, uint64_t r, uint8_t word_width = sizeof(uint_t)*8, uint16_t num_vectors = 256, uint16_t p = omp_get_max_threads()) {build([&vector](uint64_t i){return vector[i];},l,r,word_width,num_vectors,p);}
    par_var_width_vector(std::vector<uint_t> &vector, uint8_t word_width = sizeof(uint_t)*8, uint16_t num_vectors = 256, uint16_t p = omp_get_max_threads()) : par_var_width_vector(vector,0,vector.size()-1,word_width,num_vectors,p) {}

    void build(std::function<uint_t(uint64_t)> read, uint64_t l, uint64_t r, uint8_t word_width = sizeof(uint_t)*8, uint16_t num_vectors = 256, uint16_t p = omp_get_max_threads()) {
        size = r-l+1;
        num_vectors = std::max((uint64_t)1,std::min((uint64_t)num_vectors,size/1000));
        this->num_vectors = num_vectors;
        vectors.resize(num_vectors);

        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_vectors; i++) {
			uint64_t b = l + i * (size / num_vectors);
            uint64_t e = l + (i == num_vectors - 1 ? size - 1 : (i + 1) * (size / num_vectors) - 1);
            uint64_t pos_in_vector = 0;
            vectors[i].width(word_width);
            vectors[i].resize(e-b+1);

            for (uint64_t j=b; j<=e; j++) {
                vectors[i][pos_in_vector++] = read(j);
            }
        }

        bit_size = 18;

        for (uint16_t i=0; i<num_vectors; i++) {
            bit_size += vectors[i].bit_size();
        }
    }

    void decode(std::function<void(uint64_t,uint_t)> write, uint64_t l, uint64_t r, uint16_t p = omp_get_max_threads()) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t i = 0; i < num_vectors; i++) {
			uint64_t b = l + i * (size / num_vectors);
            uint64_t e = l + (i == num_vectors - 1 ? size - 1 : (i + 1) * (size / num_vectors) - 1);
            uint64_t pos_in_vector = 0;

            for (uint64_t j=b; j<=e; j++) {
                write(j,vectors[i][pos_in_vector++]);
            }
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
        in.read((char*)&num_vectors,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        vectors.resize(num_vectors);
        uint16_t num_vectors_left = num_vectors;
        uint16_t num_vectors_to_read;

        while (num_vectors_left > 0) {
            num_vectors_to_read = std::min(num_vectors_left,(uint16_t)p);

            for (uint16_t i = num_vectors - num_vectors_left; i < num_vectors - num_vectors_left + num_vectors_to_read; i++) {
                vectors[i].load(in);
            }

            #pragma omp parallel for num_threads(p)
            for (uint16_t i = num_vectors - num_vectors_left; i < num_vectors - num_vectors_left + num_vectors_to_read; i++) {
                uint64_t b = l + i * (size / num_vectors);
                uint64_t e = l + (i == num_vectors - 1 ? size - 1 : (i + 1) * (size / num_vectors) - 1);
                uint64_t pos_in_vector = 0;

                for (uint64_t j=b; j<=e; j++) {
                    write(j,vectors[i][pos_in_vector++]);
                }
                
                vectors[i].resize(0);
            }

            num_vectors_left -= num_vectors_to_read;
        }

        clear();
    }

    void clear() {
        num_vectors = 0;
        size = 0;
        vectors.clear();
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
        in.read((char*)&num_vectors,sizeof(uint16_t));
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&bit_size,sizeof(uint64_t));
        vectors.resize(num_vectors);

        for (uint16_t i=0; i<num_vectors; i++) {
            vectors[i].load(in);
        }
    }

    void serialize(std::ofstream &out) {
        out.write((char*)&num_vectors,sizeof(uint16_t));
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&bit_size,sizeof(uint64_t));

        for (uint16_t i=0; i<num_vectors; i++) {
            vectors[i].serialize(out);
        }
    }
};