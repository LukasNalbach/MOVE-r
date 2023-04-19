#pragma once

/**
 * @brief node in a dll_list
 * @tparam T value type
 */
template <typename T>
struct dll_node {
    T v; // value
    dll_node<T> *pr; // predecesor
    dll_node<T> *sc; // successor

    /**
     * @brief creates an empty dll_node
     */
    dll_node() {};

    /**
     * @brief creates a dll_node with value v
     */
    dll_node(T v) {
        this->v = v;
        pr = sc = NULL;
    };

    /**
     * @brief creates a dll_node with value v, predecessor pr and successor sc
     * @param pr a dll_node
     * @param sc a dll_node
     */
    dll_node(T v, dll_node<T>* pr, dll_node<T>* sc) {
        this->v = v;
        this->pr = pr;
        this->sc = sc;
    }

    /**
     * @brief deletes the dll_node
     */
    ~dll_node() {
        pr = sc = NULL;
    }
};

/**
 * @brief doubly linked list
 * @tparam T 
 */
template <typename T>
class dl_list {
    protected:
    dll_node<T> *hd; // first node
    dll_node<T> *tl; // last node
    uint64_t s; // size
    
    public:
    /**
     * @brief creates an empty dl_list
     */
    dl_list() {
        hd = tl = NULL;
        s = 0;
    }

    /**
     * @brief deletes all nodes of the dl_list and the lsit
     */
    ~dl_list() {
        delete_nodes();
        hd = tl = NULL;
        s = 0;
    }

    /**
     * @brief checks whether the dl_list is empty
     * @return whether the dl_list is empty
     */
    inline bool empty() {
        return s == 0;
    }

    /**
     * @brief returns the number of nodes in the dl_list
     * @return number of nodes in the dl_list
     */
    inline uint64_t size() {
        return s;
    }

    /**
     * @brief adjusts the size of the dl_list
     * @param s new size
     */
    inline void set_size(uint64_t s) {
        this->s = s;
    }

    /**
     * @brief returns the head of the dl_list
     * @return the head of the dl_list
     */
    inline dll_node<T>* head() {
        return hd;
    }

    /**
     * @brief adjusts the head of the dl_list
     * @param n new head
     */
    inline void set_head(dll_node<T> *n) {
        this->hd = n;
    }

    /**
     * @brief returns the tail of the dl_list
     * @return the tail of the dl_list
     */
    inline dll_node<T>* tail() {
        return tl;
    }

    /**
     * @brief adjusts the head of the dl_list
     * @param n new tail
     */
    inline void set_tail(dll_node<T> *n) {
        this->tl = n;
    }

    /**
     * @brief inserts a node with value v before the head of the dl_list
     * @param v value
     * @return the newly created node
     */
    inline dll_node<T>* push_front(T &&v) {
        return push_front_node(new dll_node<T>(v));
    }

    /**
     * @brief inserts a node with value v before the head of the dl_list
     * @param v value
     * @return the newly created node
     */
    inline dll_node<T>* push_front(T &v) {
        return push_front_node(new dll_node<T>(std::move(v)));
    }

    /**
     * @brief inserts a node with value v after the head of the dl_list
     * @param v value
     * @return the newly created node
     */
    inline dll_node<T>* push_back(T &&v) {
        return push_back_node(new dll_node<T>(v));
    }

    /**
     * @brief inserts a node with value v after the head of the dl_list
     * @param v value
     * @return the newly created node
     */
    inline dll_node<T>* push_back(T &v) {
        return push_back_node(new dll_node<T>(std::move(v)));
    }

    /**
     * @brief inserts a node with value v before the node n2
     * @param v value
     * @param n2 a dll_node in the dl_list
     * @return the newly created node
     */
    inline dll_node<T>* insert_before(T &&v, dll_node<T> *n) {
        return insert_before_node(new dll_node<T>(v),n);
    }

    /**
     * @brief inserts a node with value v before the node n2
     * @param v value
     * @param n2 a dll_node in the dl_list
     * @return the newly created node
     */
    inline dll_node<T>* insert_before(T &v, dll_node<T> *n) {
        return insert_before_node(new dll_node<T>(std::move(v)),n);
    }

    /**
     * @brief inserts a node with value v after the node n2
     * @param v value
     * @param n2 a dll_node in the dl_list
     * @return the newly created node
     */
    inline dll_node<T>* insert_after(T &&v, dll_node<T> *n) {
        return insert_after_node(new dll_node<T>(v),n);
    }

    /**
     * @brief inserts a node with value v after the node n2
     * @param v value
     * @param n2 a dll_node in the dl_list
     * @return the newly created node
     */
    inline dll_node<T>* insert_after(T &v, dll_node<T> *n) {
        return insert_after_node(new dll_node<T>(std::move(v)),n);
    }

    /**
     * @brief inserts the node n1 before the node n2
     * @param n1 a dll_node, that is not in the dl_list
     * @param n2 a dll_node in the dl_list, n1 != n2
     */
    inline void insert_before_node(dll_node<T> *n1, dll_node<T> *n2) {
        if (n2 == hd) {
            hd = n1;
        } else {
            n2->pr->sc = n1;
            n1->pr = n2->pr;
        }
        n1->sc = n2;
        n2->pr = n1;
        s++;
    }

    /**
     * @brief inserts the node n1 after the node n2
     * @param n1 a dll_node, that is not in the dl_list
     * @param n2 a dll_node in the dl_list, n1 != n2
     */
    inline void insert_after_node(dll_node<T> *n1, dll_node<T> *n2) {
        if (n2 == tl) {
            tl = n1;
        } else {
            n2->sc->pr = n1;
            n1->sc = n2->sc;
        }
        n1->pr = n2;
        n2->sc = n1;
        s++;
    }

    /**
     * @brief inserts the node n before the head of the dl_list
     * @param n a dll_node, that is not in the dl_list
     */
    inline void push_front_node(dll_node<T> *n) {
        if (empty()) {
            hd = tl = n;
            s = 1;
        } else {
            insert_before_node(n,hd);
        }
    }

    /**
     * @brief inserts the node n after the tail of the dl_list
     * @param n a dll_node, that is not in the dl_list
     */
    inline void push_back_node(dll_node<T> *n) {
        if (empty()) {
            hd = tl = n;
            s = 1;
        } else {
            insert_after_node(n,tl);
        }
    }

    /**
     * @brief concatenates the dl_list l to the end of the dl_list
     * @param l another dl_list
     */
    inline void concat(dl_list<T> *l) {
        if (empty()) {
            if (!l->empty()) {
                hd = l->hd;
                tl = l->tl;
                s = l->s;
                l->disconnect_nodes();
            }
        } else if (!l->empty()) {
            tl->sc = l->hd;
            l->hd->pr = tl;
            tl = l->tail();
            s += l->s;
            l->disconnect_nodes();
        }
    }

    /**
     * @brief removes the node n from the dl_list
     * @param n a dll_node in the dl_list
     */
    inline void remove_node(dll_node<T> *n) {
        s--;
        if (n == hd) {
            hd = n->sc;
        } else if (n == tl) {
            tl = n->pr;
        }
        if (n->pr != NULL) {
            n->pr->sc = n->sc;
        }
        if (n->sc != NULL) {
            n->sc->pr = n->pr;
        }
    }

    /**
     * @brief disconnects the dl_list from it's nodes
     */
    inline void disconnect_nodes() {
        hd = tl = NULL;
        s = 0;
    }

    /**
     * @brief deletes all nodes in the dl_list
     */
    void delete_nodes() {
        if (!empty()) {
            dll_node<T> *n = hd;
            for (uint64_t i=1; i<s; i++) {
                n = n->sc;
                delete n->pr;
            }
            hd = tl = NULL;
            s = 0;
        }
    }

    /**
     * @brief iterator for a dl_list
     */
    class dll_it {
        protected:
        dl_list<T> *l; // the dl_list, the iterator iterates through
        dll_node<T> *cur; // the node the iterator points to

        public:
        /**
         * @brief creates an dl_it pointing to the node n in the dl_list l
         * @param l a dl_list
         * @param n a dll_node in l
         */
        dll_it(dl_list<T> *l, dll_node<T> *n) {
            this->l = l;
            this->cur = n;
        }

        /**
         * @brief deletes the iterator
         */
        ~dll_it() {
            l = NULL;
            cur = NULL;
        }

        /**
         * @brief checks whether the iterator can iterate forward
         * @return whether it can iterate forward
         */
        inline bool has_next() {
            return cur->sc != NULL;
        }

        /**
         * @brief checks whether the iterator can iterate backward
         * @return whether it can iterate backward
         */
        inline bool has_prev() {
            return cur->pr != NULL;
        }

        /**
         * @brief returns the value of the node, the iterator points to
         * @return the node, the iterator points to
         */
        inline dll_node<T>* current() {
            return cur;
        }

        /**
         * @brief iterates forward, has_next() must return true
         * @return the node the iterator points to after iterating forward
         */
        inline dll_node<T>* next() {
            cur = cur->sc;
            return cur;
        }

        /**
         * @brief iterates forward, has_pred() must return true
         * @return the node the iterator points to after iterating backward
         */
        inline dll_node<T>* previous() {
            cur = cur->pr;
            return cur;
        }

        /**
         * @brief points the iterator to the node n
         * @param n a dll_node in l
         */
        inline void set(dll_node<T> *n) {
            cur = n;
        }
    };

    /**
     * @brief returns an iterator pointing to the node n
     * 
     * @param n a dll_node in the dl_list
     * @return an iterator
     */
    inline dl_list<T>::dll_it iterator(dll_node<T> *n) {
        return dl_list<T>::dll_it(this,n);
    }

    /**
     * @brief returns an iterator pointing to the head of the list
     * @return an iterator
     */
    inline dl_list<T>::dll_it iterator() {
        return dl_list<T>::dll_it(this,hd);
    }
};