#pragma once

#include <queue>
#include <omp.h>
#include <concurrentqueue.h>
#include <misc/avl_tree.hpp>
#include <misc/dl_list.hpp>
#include <misc/dg_io_nc.hpp>

// ############################# V2/3/4 #############################

template <typename uint_t> using pair_list_node = dll_node<std::pair<uint_t,uint_t>>;
template <typename uint_t> using pair_list = dl_list<std::pair<uint_t,uint_t>>;

template <typename uint_t> using pair_tree_node = avl_node<pair_list_node<uint_t>>;
template <typename uint_t> using pair_tree = avl_tree<pair_list_node<uint_t>>;

// ############################# V2 #############################

template <typename uint_t> using te_pair = std::pair<pair_list_node<uint_t>*,pair_tree_node<uint_t>*>;
template <typename uint_t> using te_node = avl_node<te_pair<uint_t>>;
template <typename uint_t> using te_tree = avl_tree<te_pair<uint_t>>;

// ############################# V3/V4 PARALLEL #############################

template <typename uint_t> using q_pair = std::pair<pair_list_node<uint_t>*,pair_list_node<uint_t>*>;

// ############################# V3 PARALLEL #############################

template <typename uint_t> using qt_v3 = std::vector<std::vector<std::queue<q_pair<uint_t>>>>;

// ############################# V4 PARALLEL #############################

template <typename uint_t> using qt_v4 = std::vector<std::vector<moodycamel::ConcurrentQueue<q_pair<uint_t>>>>;

/**
 * @brief balances a disjoint interval sequence
 * @tparam uint_t (integer) type of the interval starting positions
 */
template <typename uint_t>
class mds_builder {
    public:

    /**
     * @brief 
     * 
     * @param I 
     * @param n 
     * @param k 
     * @param D_p 
     * @param Q 
     * @param k_ 
     * @param l_max
     * @param a 
     * @param p 
     * @param log 
     * @param mf 
     */
    static void build_mds(
        std::pair<uint_t,uint_t>*& I,
        uint64_t n,
        uint64_t k,
        uint_t*& D_p,
        std::pair<uint_t,uint_t>*& D_offsidx,
        uint64_t &k_,
        uint64_t l_max,
        uint16_t a,
        uint16_t p,
        bool log,
        std::ofstream &mf
    ) {
        mds_builder<uint_t>(I,n,k,D_p,D_offsidx,k_,a,p,l_max,log,mf);
    };

    mds_builder() = delete;
    mds_builder(const mds_builder&) = delete;

    /**
     * @brief builds the move datastructure md
     * @param I disjoint interval sequence
     * @param D_p input interval starting positions
     * @param Q output interval starting positions
     * @param n n = p_{k-1} + d_{k-1}, k <= n
     * @param k k = |I|
     * @param l_max maximum input interval length
     * @param a balancing parameter, restricts size increase to the factor
     *          (1+1/(a-1)) and restricts move query runtime to 2a, 2 <= a
     * @param p number of threads to use
     * @param log enables log messages during build process
     * @param mf output stream to write runtime and space usage to if log is enabled
     */
    mds_builder(
        std::pair<uint_t,uint_t>*& I,
        uint64_t n,
        uint64_t k,
        uint_t*& D_p,
        std::pair<uint_t,uint_t>*& D_offsidx,
        uint64_t &k_,
        uint16_t a,
        uint16_t p,
        uint64_t l_max,
        bool log, std::ofstream &mf
    ) : I(I), D_p(D_p), D_offsidx(D_offsidx), k_(k_), mf(mf) {
        this->n = n;
        this->k = k;
        this->a = a;
        this->l_max = l_max;
        this->log = log;

        if (p > 1 && 1000*p > k) {
            p = std::max((uint64_t)1,k/1000);
            if (log) std::cout << std::endl << "warning: p > k/1000, setting p to k/1000 ~ " << std::to_string(p);
        }

        this->p = p;

        omp_set_num_threads(p);

        if constexpr (v == 1) {
            balance_v1();
        } else {
            balance_v2_v3_v4();
        }

        #ifndef NDEBUG
        verify_correctness();
        #endif

        free(pi);
        pi = NULL;
        free(D_q);
        D_q = NULL;
    }

    /**
     * @brief deletes the mds_builder
     */
    ~mds_builder() {}

    /**
     * @brief builds D_offs and D_idx
     */
    void build_doffsidx();

    /**
     * @brief verifies the correctness of the resulting interval sequence
     */
    void verify_correctness();

    // ############################# VARIABLES #############################

    protected:
    std::pair<uint_t,uint_t>*& I; // disjoint interval sequence
    uint_t*& D_p; // input interval starting positions
    std::pair<uint_t,uint_t>*& D_offsidx;
    uint_t* D_q; // output interval starting positions
    uint_t* pi = NULL; // permutation storing the order of the output intervals
    uint64_t n; // maximum value, n = p_{k-1} + d_{k-1}, k <= n
    uint64_t k; // number of intervals in the (possibly unbalanced) inteval sequence I, 0 < k
    uint64_t &k_; // number of intervals in the balanced inteval sequence B(I), 0 < k <= k'
    uint16_t a; // balancing parameter, restricts size increase to the factor (1+1/(a-1)), 2 <= a
    uint16_t p; // number of threads to use
    static constexpr uint8_t v = 3; // balancing method version
    uint64_t l_max; // maximum input interval length
    bool log;
    std::ofstream &mf; // measurement file
    size_t baseline; // baseline memory allocation
    std::chrono::steady_clock::time_point time; // time point of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time point of the start of the build process

    // ############################# V1 #############################

    /**
     * @brief builds the move datastructure md
     */
    void balance_v1();

    // ############################# V2/V3/V4 #############################

    /** 
     * @brief [0..p-1] doubly linked lists; L_in[i_p] stores the pairs (p_i,q_i) in ascending order of p_i,
     *        where s[i_p] <= p_i < s[i_p+1] and i_p in [0..p-1]. L_in[0]L_in[1]...L_in[p-1] = I.
     */
    std::vector<pair_list<uint_t>> L_in;
    /**
     * @brief [0..p-1] avl trees; T_out[i_p] stores nodes of lists in L_in in ascending order of q_i,
     *        for each pair (p_i,q_i), s[i_p] <= q_i < s[i_p+1] holds, with i_p in [0..p-1].
     */
    std::vector<pair_tree<uint_t>> T_out;
    /**
     * @brief [0..p-1] section start positions in the range [0..n], 0 = s[0] < s[1] < ... < s[p-1] = n.
     *        Before building T_out, s is chosen so that |L_in[0]| + |T_out[0]| ~ |L_in[1]| + |T_out[1]|
     *        ~ ... ~ |L_in[p-1]| + |T_out[p-1]|, that is, if
     *        s[i_p] = min {s' \in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * \lfloor 2k/p \rfloor, where 
                                x[i_p] = min {x' \in [0,k-1], , s.t. p_{x'} >= s'} and
                                u[i_p] = min {u' \in [0,k-1], s.t. q_{\pi[u']} >= s'}
	          } holds.
     */
    std::vector<uint_t> s;
    /**
     * @brief stores the nodes in L_in[0..p-1] and T_out[0..p-1] which they were initially created with
     */
    pair_tree_node<uint_t>* nodes = NULL;
    /**
     * @brief [0..p-1] new_nodes[i_p] stores the newly created nodes in L_in[0..p-1] and T_out[0..p-1],
     *        which were created by thread i_p.
     */
    std::vector<dg_io_nc<pair_tree_node<uint_t>>> new_nodes;

    /**
     * @brief builds the move datastructure md
     */
    void balance_v2_v3_v4() {
        if (log) {
            baseline = malloc_count_current() - sizeof(I[0])*k;
            time = now();
            time_start = time;
            std::cout << std::endl;
        }

        // Build L_in and T_out
        build_lin_tout();

        // Choose the correct balancing algorithm
        if (p > 1) {
            if constexpr (v == 3) {
                balance_v3_par();
            } else {
                balance_v4_par();
            }
        } else {
            if constexpr (v == 2) {
                balance_v2_seq();
            } else {
                balance_v3_seq();
            }
        }

        // Build D_p and D_q
        build_dp_dq();

        // Build D_offs and D_idx
        build_doffsidx();

        if (log) {
            log_message("move datastructure built");
            if (mf.is_open()) {
                mf << " time_build_mds=" << time_diff_ns(time_start);
            }
        }
    }

    // ############################# V2/V3/V4 SEQUENTIAL/PARALLEL #############################

    /**
     * @brief builds L_in[0..p-1] and T_out[0..p-1] out of the disjoint interval sequence I
     */
    void build_lin_tout();

    /**
     * @brief inserts the pairs in L_in[0..p-1] into D_pair and writes D_idx
     */
    void build_dp_dq();

    /**
     * @brief returns the length of the input/output interval starting at p_i/q_i
     * @param pln (p_i,q_i)
     * @return length of the input/output interval starting at p_i/q_i
     */
    inline uint_t interval_length_seq(pair_list_node<uint_t> *pln);

    /**
     * @brief checks if the output interval [q_j, q_j + d_j) is unbalanced and iterates
     *        pln_IpI_ tothe lastoutput interval connected to it in the permutation graph
     * @param pln_IpI_ (p_{i+i_},q_{i+i_}), [p_i, p_i + d_i) must be the first input
     *                 interval connected to [q_j, q_j + d_j) in the permutation graph
     * @param i_ 1 <= i_ <= 2a
     * @param ptn_J (p_j,q_j)
     * @param ptn_J_nxt (p_{j'},q_{j'}), with q_j + d_j = q_{j'}
     * @return (p_{i+a},q_{i+a}) if [q_j, q_j + d_j) is unbalanced, else NULL
     */
    inline pair_list_node<uint_t>* is_unbalanced(
        pair_list_node<uint_t> **pln_IpI_,
        uint_t* i_,
        pair_tree_node<uint_t> *ptn_J,
        pair_tree_node<uint_t> *ptn_J_nxt = NULL
    );

    // ############################# V2 #############################

    /**
     * @brief balances the disjoint interval sequence in L_in[0] and T_out[0] sequentially
     */
    void balance_v2_seq();

    // ############################# V3 SEQUENTIAL #############################

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all unbalanced output intervals
     *        starting before or at q_u that have become unbalanced in the process
     * @param pln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *                to [q_j, q_j + d_j) in the permutation graph
     * @param ptn_J (p_j,q_j), [q_j, q_j + d_j) must be the first unbalanced output interval
     * @param q_u starting position of an output interval
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline pair_tree_node<uint_t>* balance_upto_v3_seq(
        pair_list_node<uint_t> *pln_IpA,
        pair_tree_node<uint_t> *ptn_J,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in[0] and T_out[0] sequentially
     */
    void balance_v3_seq();

    // ############################# V3 PARALLEL #############################

    /**
     * @brief balances the output interval [q_j, q_j + d_j) by inserting the newly created pair into
     *        T_out[i_p] and Q[0..p-1][i_p]
     * @param Q reference to Q
     * @param pln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *                to [q_j, q_j + d_j) in the permutation graph
     * @param ptn_J (p_j,q_j), [q_j, q_j + d_j) must be the first unbalanced output interval starting 
     *              in [s[i_p]..q_u]
     * @param ptn_J_next (p_j',q_j'), [q_j', q_j' + d_j') must be the first output interval starting
     *                   after [q_j, q_j + d_j)
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]]
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline pair_tree_node<uint_t>* balance_upto_v3_par(
        qt_v3<uint_t> &Q,
        pair_list_node<uint_t>* pln_IpA,
        pair_tree_node<uint_t>* ptn_J,
        pair_tree_node<uint_t>* ptn_J_nxt,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in[0..p-1] and T_out[0..p-1] in parallel
     */
    void balance_v3_par();

    // ############################# V4 PARALLEL #############################

    /**
     * @brief balances the output interval [q_j, q_j + d_j) by inserting the newly created pair into
     * T_out[i_p] and Q[0..p-1][i_p]
     * @param Q reference to Q
     * @param pln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *                to [q_j, q_j + d_j) in the permutation graph
     * @param ptn_J (p_j,q_j), [q_j, q_j + d_j) must be the first unbalanced output interval starting
     *              in [s[i_p]..q_u]
     * @param ptn_J_next (p_j',q_j'), [q_j', q_j' + d_j') must be the first output interval starting
     *                   after [q_j, q_j + d_j)
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]]
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline pair_tree_node<uint_t>* balance_upto_v4_par(
        qt_v4<uint_t> &Q,
        pair_list_node<uint_t> *pln_IpA,
        pair_tree_node<uint_t> *ptn_J,
        pair_tree_node<uint_t>* ptn_J_nxt,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in[0..p-1] and T_out[0..p-1] in parallel
     */
    void balance_v4_par();
};

#include "mds_builder/v1.cpp"
#include "mds_builder/common.cpp"
#include "mds_builder/build/lin_tout.cpp"
#include "mds_builder/build/dp_dq.cpp"
#include "mds_builder/build/doffs_didx.cpp"
#include "mds_builder/balance/v2_seq.cpp"
#include "mds_builder/balance/v3_par.cpp"
#include "mds_builder/balance/v3_seq.cpp"
#include "mds_builder/balance/v4_par.cpp"