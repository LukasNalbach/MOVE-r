#pragma once

/**
 * @brief dynamically-growing insert-only no-copy datastructure
 * @tparam T value type
 */
template <typename T>
class dg_io_nc {
    protected:
    std::vector<std::vector<T>> vectors; // vectors that store the elements

    public:
    /**
     * @brief creates an empty datastructure
     */
    dg_io_nc() {}

    /**
     * @brief creates an empty datastructure with a certain amount of elements reserved
     * @param size initially reserved number of elements
     */
    dg_io_nc(uint64_t size) {
        vectors = std::vector<std::vector<T>>();
        vectors.reserve(50);
        vectors.emplace_back(std::vector<T>());
        vectors.back().reserve(size);
    }

    ~dg_io_nc() {
        clear();
    }

    /**
     * @brief clears all vectors
     */
    inline void clear() {
        vectors.clear();
        vectors.shrink_to_fit();
    }

    /**
     * @brief inserts an element into the datastructure, doubles the number of elements reserved
     *        reserved by the datastructure if it is full
     * @param v element
     * @return pointer to the element in the datastructure
     */
    inline T* emplace_back(T &&v) {
        if (vectors.back().size() == vectors.back().capacity()) {
            size_t new_capacity = 2*vectors.back().capacity();
            vectors.emplace_back(std::vector<T>());
            vectors.back().reserve(new_capacity);
        }

        vectors.back().emplace_back(v);
        return &(vectors.back().back());
    }

    /**
     * @brief inserts an element into the datastructure, doubles the number of elements reserved
     *        reserved by the datastructure if it is full
     * @param v element
     * @return pointer to the element in the datastructure
     */
    inline T* emplace_back(T &v) {
        return emplace_back(std::move(v));
    }
};