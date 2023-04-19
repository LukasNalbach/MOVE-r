#pragma once

#include <sdsl/uint256_t.hpp>
#include <sdsl/bit_vectors.hpp>
#include <misc/avl_tree.hpp>
#include <misc/utils.hpp>

class huff_string {
    protected:
    uint64_t size = 0;
    uint16_t num_dist_chars = 0;
    uint64_t huff_data_bit_size = 0;
    sdsl::bit_vector huff_tree_bps;
    std::vector<uint8_t> leaf_char_order;
    std::vector<uint8_t,default_init_allocator<uint8_t>> huff_data;

    struct huff_tree_node_encode {
        uint16_t character;
        uint64_t frequency;
        huff_tree_node_encode* parent = NULL;
        huff_tree_node_encode* left_child = NULL;
        huff_tree_node_encode* right_child = NULL;
    };

    struct huff_tree_node_decode {
        uint16_t character = 256;
        uint16_t left_child = 0;
        uint16_t right_child = 0;
    };

    public:
    huff_string() = default;
    huff_string(const huff_string&) = default;
    huff_string(std::vector<uint8_t> &string, uint64_t l, uint64_t r) {build([&string](uint64_t i){return string[i];},l,r);}
    huff_string(std::vector<uint8_t> &string) : huff_string(string,0,string.size()-1) {}

    void build(std::function<uint8_t(uint64_t)> read, uint64_t l, uint64_t r) {
        std::vector<uint64_t> char_freq(256,0);
        std::vector<avl_node<huff_tree_node_encode>*> avl_nodes;

        for (uint64_t i=l; i<=r; i++) {
            char_freq[read(i)]++;
        }

        avl_tree<huff_tree_node_encode> huff_tree(
            [](auto n1, auto n2){return n1.frequency < n2.frequency || (n1.frequency == n2.frequency && n1.character < n2.character);},
            [](auto n1, auto n2){return n1.frequency > n2.frequency || (n1.frequency == n2.frequency && n1.character > n2.character);},
            [](auto n1, auto n2){return n1.frequency == n2.frequency && n1.character == n2.character;}
        );;

        for (uint16_t i=0; i<256; i++) {
            if (char_freq[i] != 0) {
                avl_nodes.emplace_back(huff_tree.insert_or_update(huff_tree_node_encode{(uint8_t)i,char_freq[i]}));
            }
        }

        char_freq.clear();
        num_dist_chars = huff_tree.size();
        avl_node<huff_tree_node_encode> *tmp_left,*tmp_right,*tmp_new;
        uint16_t next_huff_node_id = 256;

        while (huff_tree.size() > 1) {
            tmp_left = huff_tree.minimum();
            huff_tree.remove_node(huff_tree.minimum());
            tmp_right = huff_tree.minimum();
            huff_tree.remove_node(huff_tree.minimum());

            tmp_new = huff_tree.insert_or_update(huff_tree_node_encode{
                next_huff_node_id++,
                tmp_left->v.frequency+tmp_right->v.frequency,
                NULL,
                &tmp_left->v,
                &tmp_right->v
            });
            
            avl_nodes.emplace_back(tmp_new);
            tmp_left->v.parent = &tmp_new->v;
            tmp_right->v.parent = &tmp_new->v;
        }

        leaf_char_order.resize(num_dist_chars,0);
        huff_tree_bps.bit_resize(2*avl_nodes.size());
        huff_tree_node_encode* cur_huff_node = &huff_tree.minimum()->v;
        uint8_t cur_tree_depth = 0;
        uint16_t pos_in_huff_tree_bps = 0;

        while (cur_huff_node->right_child != NULL) {
            cur_huff_node = cur_huff_node->right_child;
            cur_tree_depth++;
            huff_tree_bps[pos_in_huff_tree_bps++] = 1;
        }

        sdsl::uint256_t cur_huff_code = ((sdsl::uint256_t) -1) << (256 - cur_tree_depth);
        std::vector<std::pair<uint8_t,sdsl::uint256_t>> map_to_huff_code(256);
        map_to_huff_code[cur_huff_node->character] = std::pair<uint8_t,sdsl::uint256_t>{cur_tree_depth,cur_huff_code};
        leaf_char_order[0] = cur_huff_node->character;
        huff_data_bit_size += cur_tree_depth*cur_huff_node->frequency;

        for (uint16_t i=1; i<num_dist_chars; i++) {
            do {
                if (cur_huff_node->left_child != NULL) {
                    cur_huff_node = cur_huff_node->left_child;
                    cur_tree_depth++;
                    huff_tree_bps[pos_in_huff_tree_bps++] = 1;

                    while (cur_huff_node->right_child != NULL) {
                        cur_huff_node = cur_huff_node->right_child;
                        cur_tree_depth++;
                        cur_huff_code |= ((sdsl::uint256_t) 1) << (256 - cur_tree_depth);
                        huff_tree_bps[pos_in_huff_tree_bps++] = 1;
                    }
                } else {
                    while (cur_huff_node->parent != NULL && cur_huff_node == cur_huff_node->parent->left_child) {
                        cur_huff_node = cur_huff_node->parent;
                        cur_tree_depth--;
                        huff_tree_bps[pos_in_huff_tree_bps++] = 0;
                    }

                    cur_huff_node = cur_huff_node->parent;
                    cur_tree_depth--;
                    cur_huff_code = cur_huff_code & (((sdsl::uint256_t) -1) << (256 - cur_tree_depth));
                    huff_tree_bps[pos_in_huff_tree_bps++] = 0;
                }
            } while (cur_huff_node->character >= 256);

            map_to_huff_code[cur_huff_node->character] = std::pair<uint8_t,sdsl::uint256_t>{cur_tree_depth,cur_huff_code};
            leaf_char_order[i] = cur_huff_node->character;
            huff_data_bit_size += cur_tree_depth*cur_huff_node->frequency;
        }

        huff_tree_bps[pos_in_huff_tree_bps++] = 0;
        huff_tree.disconnect_nodes();
        avl_nodes.clear();
        size = r-l+1;
        huff_tree_bps.bit_resize(pos_in_huff_tree_bps);
        huff_data.resize(huff_data_bit_size/8+1);
        std::memset(&huff_data[0],0,huff_data.size());
        uint64_t idx_in_huff_data = 0;
        uint64_t offs_in_huff_data = 0;
        uint64_t offs_in_code = 0;
        uint64_t current_block_length;
        uint64_t block_start_pos_in_string = l;
        uint64_t written_bits;
        uint64_t block_size = std::min((uint64_t) 256, r-l+1);
        std::pair<uint8_t,sdsl::uint256_t> huff_pair;
        std::vector<std::pair<uint8_t,sdsl::uint256_t>> precomputed_codes;
        precomputed_codes.reserve(block_size);

        while (block_start_pos_in_string <= r) {
            current_block_length = std::min(block_size, r - block_start_pos_in_string + 1);
            
            for (uint64_t i=0; i<current_block_length; i++) {
                precomputed_codes.emplace_back(map_to_huff_code[read(block_start_pos_in_string + i)]);
            }

            for (uint64_t i=0; i<current_block_length; i++) {
                huff_pair = precomputed_codes[i];
                offs_in_code = 0;

                while (offs_in_code < huff_pair.first) {
                    huff_data[idx_in_huff_data] |= (uint8_t) (huff_pair.second >> (uint8_t) ((248 + offs_in_huff_data) - offs_in_code));
                    written_bits = std::min((uint64_t) 8 - offs_in_huff_data, (uint64_t) huff_pair.first - offs_in_code);
                    offs_in_code += written_bits;
                    offs_in_huff_data += written_bits;

                    if (offs_in_huff_data >= 8) {
                        idx_in_huff_data++;
                        offs_in_huff_data -= 8;
                    }
                }
            }
            
            precomputed_codes.clear();
            block_start_pos_in_string += current_block_length;
        }
    }

    void decode(std::function<void(uint64_t,uint8_t)> write, uint64_t l, uint64_t r) {
        std::vector<huff_tree_node_decode> huff_tree_nodes;
        std::vector<uint16_t> path_to_cur_node;
        huff_tree_nodes.reserve(num_dist_chars+(num_dist_chars/2+1)*log2_clz<uint32_t>((num_dist_chars/2)+1)-1);
        path_to_cur_node.reserve(64);
        huff_tree_nodes.emplace_back(huff_tree_node_decode{});
        uint32_t huff_tree_root = 0;
        uint16_t cur_node = huff_tree_root;
        uint16_t pos_in_huff_tree_bps = 0;

        while (huff_tree_bps[pos_in_huff_tree_bps] == 1) {
            pos_in_huff_tree_bps++;
            path_to_cur_node.emplace_back(cur_node);
            huff_tree_nodes.emplace_back(huff_tree_node_decode{});
            huff_tree_nodes[cur_node].right_child = huff_tree_nodes.size()-1;
            cur_node = huff_tree_nodes[cur_node].right_child;
        }

        huff_tree_nodes[cur_node].character = leaf_char_order[0];

        for (uint16_t i=1; i<num_dist_chars; i++) {
            while (huff_tree_bps[pos_in_huff_tree_bps] == 0) {
                pos_in_huff_tree_bps++;
                cur_node = path_to_cur_node.back();
                path_to_cur_node.pop_back();
            }

            path_to_cur_node.emplace_back(cur_node);
            huff_tree_nodes.emplace_back(huff_tree_node_decode{});
            huff_tree_nodes[cur_node].left_child = huff_tree_nodes.size()-1;
            cur_node = huff_tree_nodes[cur_node].left_child;
            pos_in_huff_tree_bps++;

            while (huff_tree_bps[pos_in_huff_tree_bps] == 1) {
                pos_in_huff_tree_bps++;
                path_to_cur_node.emplace_back(cur_node);
                huff_tree_nodes.emplace_back(huff_tree_node_decode{});
                huff_tree_nodes[cur_node].right_child = huff_tree_nodes.size()-1;
                cur_node = huff_tree_nodes[cur_node].right_child;
            }

            huff_tree_nodes[cur_node].character = leaf_char_order[i];
        }

        path_to_cur_node.clear();
        uint64_t idx_in_huff_data = 0;
        uint64_t offs_in_huff_data = 0;
        sdsl::uint256_t code = 0;
        uint64_t offs_in_code = 0;
        int16_t shift_left_huff_data;
        uint8_t offs_from_right_of_code;

        for (uint64_t i=l; i<=r; i++) {
            while (offs_in_code < 256) {
                shift_left_huff_data = (248 - offs_in_code) + offs_in_huff_data;

                if (shift_left_huff_data >= 0) {
                    code |= ((sdsl::uint256_t) huff_data[idx_in_huff_data]) << shift_left_huff_data;
                    offs_in_code += 8 - offs_in_huff_data;
                    idx_in_huff_data++;
                    offs_in_huff_data = 0;
                } else {
                    code |= ((sdsl::uint256_t) huff_data[idx_in_huff_data]) >> - shift_left_huff_data;
                    offs_in_huff_data += 256 - offs_in_code;
                    offs_in_code = 256;
                }
            }

            offs_from_right_of_code = 255;
            cur_node = huff_tree_root;
            
            while (huff_tree_nodes[cur_node].character == 256) {
                if (((uint8_t) (code >> offs_from_right_of_code)) & ((uint8_t) 1)) {
                    cur_node = huff_tree_nodes[cur_node].right_child;
                } else {
                    cur_node = huff_tree_nodes[cur_node].left_child;
                }
                
                offs_from_right_of_code--;
            }

            write(i,huff_tree_nodes[cur_node].character);
            code = code << 255 - offs_from_right_of_code;
            offs_in_code = offs_from_right_of_code + 1;
        }
    }

    void decode(std::vector<uint8_t> &string, uint64_t l, uint64_t r) {
        decode([&string](uint64_t i, uint8_t c){string[i] = c;},l,r);
    }

    void decode(std::vector<uint8_t> &string) {
        decode(string,0,string.size()-1);
    }

    std::vector<uint8_t> decode() {
        std::vector<uint8_t> string(size);
        decode(string);
        return string;
    }

    void clear() {
        size = 0;
        num_dist_chars = 0;
        huff_data_bit_size = 0;
        huff_tree_bps.bit_resize(0);
        leaf_char_order.clear();
        huff_data.clear();
    }

    uint64_t get_size() {
        return size;
    }

    uint64_t get_bit_size() {
        return (18+num_dist_chars+huff_data_bit_size/8+1)*8+huff_tree_bps.bit_size();
    }

    float get_compression_ratio() {
        return (size*8)/(float)get_bit_size();
    }

    void load(std::ifstream &in) {
        in.read((char*)&size,sizeof(uint64_t));
        in.read((char*)&num_dist_chars,sizeof(uint16_t));
        in.read((char*)&huff_data_bit_size,sizeof(uint64_t));
        huff_tree_bps.load(in);
        leaf_char_order.resize(num_dist_chars);
        in.read((char*)&leaf_char_order[0],num_dist_chars);
        huff_data.resize(huff_data_bit_size/8+1);
        in.read((char*)&huff_data[0],huff_data_bit_size/8+1);
    }

    void serialize(std::ofstream &out) {
        out.write((char*)&size,sizeof(uint64_t));
        out.write((char*)&num_dist_chars,sizeof(uint16_t));
        out.write((char*)&huff_data_bit_size,sizeof(uint64_t));
        huff_tree_bps.serialize(out);
        out.write((char*)&leaf_char_order[0],num_dist_chars);
        out.write((char*)&huff_data[0],huff_data_bit_size/8+1);
    }
};