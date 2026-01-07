

#ifndef BP_HPP
#define BP_HPP





// ********************* inlcudes ************************

// standard lib
#include <limits.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <bitset>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include "sys/time.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <random>

// pthread
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>

#include "parallel_algorithms.hpp"

// gurobi
#include "gurobi_c++.h"



namespace pa = parallel_algorithms;
using std::size_t;


static std::default_random_engine e;

inline int rndint(int a, int b) {
	std::uniform_int_distribution<int> cc(a, b);
	return cc(e);
}

inline void shuffle(std::vector<int> & list) {
    std::shuffle(list.begin(), list.end(), e);
}



// ********************** versions **********************

#define CODE_VERSION "20251110, v6.6"




#define USE_VER (66)


#ifndef WORKERS
#define WORKERS (32)
#endif

#define SHOW_DETAILS (0)
#define CHECK_CORRECT (1)

// # type:
// # -1: NULL
// # 0: xor
// # 1: nxor
// # 2: and
// # 3: or
#define TYPE_NULL (-1)
#define TYPE_XOR (0)
#define TYPE_NXOR (1)
#define TYPE_AND (2)
#define TYPE_OR (3)


#define ENUM_BOUND (5)

#define GUROBI_THREADS (4)

#define TIE_BREAKING_STRA "NONE"
#if STRATEGY == 1
    #undef TIE_BREAKING_STRA
    #define TIE_BREAKING_STRA "MN"
#elif STRATEGY == 2
    #undef TIE_BREAKING_STRA
    #define TIE_BREAKING_STRA "ASM"
#elif STRATEGY == 3
    #undef TIE_BREAKING_STRA
    #define TIE_BREAKING_STRA "SOM"  
#endif

#define TAR_STRATEGY_STRA "NONE"
#if TARGET_COND == 1
    #undef TAR_STRATEGY_STRA
    #define TAR_STRATEGY_STRA "AT"
#elif TARGET_COND == 2
    #undef TAR_STRATEGY_STRA
    #define TAR_STRATEGY_STRA "ST"
#endif

#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef EXPAND_STRINGIFY
#define EXPAND_STRINGIFY(x) STRINGIFY(x)
#endif



inline void print_paramters() {
    std::cout << std::endl << "\t++++++++++++++++++++ paramters ++++++++++++++++++++" << std::endl;
    std::cout << "\tProgram version: " << CODE_VERSION << std::endl;
    std::cout << "\tExecute codes: " << USE_VER << std::endl << std::endl;

    std::cout << "\tOps file = " << EXPAND_STRINGIFY(FILENAME) << std::endl;
    std::cout << "\tDimension = " << DIMENSION << std::endl;
    std::cout << "\tTie breaking strategy = " << TIE_BREAKING_STRA << std::endl;
    std::cout << "\tTarget strategy = " << TAR_STRATEGY_STRA << std::endl;
    
    std::cout << "\tNL_NUM = " << NL_NUM << std::endl << std::endl;

    std::cout << "\tDEPTH_LIMIT = " << DEPTH_LIMIT << std::endl << std::endl;

    std::cout << "\tCHECK_CORRECT = " << CHECK_CORRECT << std::endl;
    std::cout << "\tSHOW_DETAILS = " << SHOW_DETAILS << std::endl << std::endl;

    std::cout << "\tWORKERS = " << WORKERS << std::endl;
    std::cout << "\tGUROBI_THREADS = " << GUROBI_THREADS << std::endl << std::endl;    
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;    
}








//////////////////////////////////////////////////////////////


#include "time.inl"


#define DISPLAY_GAP (4)

typedef uint64_t                     word64_t;
typedef int                          wind_t;


const wind_t DIM_LEN = std::to_string(DIMENSION + NL_NUM).size();


// b[i < DIM] refers to ei.
// b[DIM] means 1.
// b[i > DIM] means nonlinear g
typedef std::bitset<DIMENSION + NL_NUM> vec_t;
typedef std::pair<wind_t, wind_t> pair_t;
typedef std::vector<wind_t> idx_t;


struct VectorComparator {
    bool operator()(const vec_t & a, const vec_t & b) const {
        if (a == b) return false;
        if (a.count() != b.count()) return a.count() < b.count();

        for (wind_t i = 0; i < DIMENSION; ++i) {
            if (a[i] == 0 && b[i] == 1) return false;
            if (a[i] == 1 && b[i] == 0) return true;
        }
        return false;
    }
};




#define PRINT_ENTER() { std::cout << std::endl; }

inline void vector_print(vec_t const & v, bool enter = true) { 

    for (wind_t i = 0; i < DIMENSION + NL_NUM; ++i) {
        if (v.test(i)) std::cout << " 1";
        else std::cout << " _";

        // has one
        if (i < DIMENSION + NL_NUM && (i + 1) % DISPLAY_GAP == 0) std::cout << " ";
    }
    if (enter) std::cout << std::endl;
}


void pair_print(pair_t const & pa, bool const & enter = false) { 
    std::cout << "(" << pa.first << ", " << pa.second << ")";
    if (enter) std::cout << std::endl;
}



void linear_print(idx_t const & lin, bool enter) {
    for (auto & idx : lin) std::cout << " " << idx;
    if (enter) std::cout << std::endl;
}


inline std::string pad_to_len(wind_t const & num, wind_t const & tar_len) { 
    std::string st = std::to_string(num);
    while (st.size() < tar_len) st = " " + st;
    return st;
}

inline std::string pad_to_len(wind_t const & tar_len) { 
    std::string st = "";
    while (st.size() < tar_len) st = " " + st;
    return st;
}



inline word64_t dvalue(wind_t const & depth) {
    return (word64_t(1) << depth);
}



inline wind_t upper_log2(word64_t const & value) {
    wind_t d = 0;
    while (dvalue(d) < value) ++d;
    return d;
}


struct op_t {
    wind_t dst = -1; // dst = dst type src
    wind_t src_0 = -1;
    wind_t src_1 = -1;
    wind_t type = TYPE_NULL;
    wind_t tar_idx = -1;
};





struct bp_param_t {
    bool display = true;
    bool greedy = true;
};


void param_print(bp_param_t const & param) { 
    PRINTF_STAMP("paramters: display %d, greedy %d\n", param.display, param.greedy);
}


struct del_tail_idxs_t {
    wind_t current_idx = -1;
    std::vector<wind_t> del_idxs;
    std::vector<wind_t> tail_idxs;
};


inline uint64_t vector_to_uint64(vec_t const & v) { 
    uint64_t res = 0;
    for (wind_t j = 0; j < DIMENSION; ++j) {
        if (j >= 64) break;
        if (v.test(j))
            res ^= ((uint64_t)1 << (DIMENSION - 1 - j));
    }
    return res;
}

inline void vector_code_print(vec_t const & v, bool enter) {
    vector_print(v, false);
    std::cout << " code: " << vector_to_uint64(v);
    if (enter) std::cout << std::endl;
}





//////////////////////////////////////////////////////////////


#include "matrix.h"



#endif 
