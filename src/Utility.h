//
// Created by ssunah on 11/30/16.
//

#ifndef FSM_UTILITY_H
#define FSM_UTILITY_H

#include "Graph.h"
#include "Config.h"
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>
#include <x86intrin.h>

#define PACK32TO64(high, low) ((uint64_t)(high)) << 32 | (low)
#define UNPACK64HIGH(v) (uint32_t)(((v) & 0xFFFFFFFF00000000LL) >> 32)
#define UNPACK64LOW(v) (uint32_t)((v) & 0xFFFFFFFFLL)

class Utility {
private:
    Utility() {}

public:
    // Refer to paper: An O(m) Algorithm for Cores Decomposition of Networks.
    // The difference is that the vertices in our implementation are ranged from 0 to n-1.
    static void GetKCore(const Graph *graph, unsigned int *core_table);

    // Find all the automorphisms of a graph.
    static void FindAutomorphism(const Graph *graph, const unsigned int* enumerate_order, const unsigned int* pivot,
                                 const unsigned int enumerate_length, vector<vector<unsigned int> > &aut);

    // Time function.
    static timeval GetTime();
    static double TimeDiffInSeconds(const timeval start, const timeval end);

    // Common neighbors computing.
    static const unsigned int GallopingSearch(const unsigned int *src, const unsigned int begin, const unsigned int end,
                                              const unsigned int target);


    static const unsigned int FirstIndexGreaterThan(const unsigned int *src, const unsigned int begin, const unsigned int end, const unsigned int target);
    static bool Contain(const unsigned int* src, const unsigned int begin, const unsigned int end, const unsigned int target);
    static const unsigned int BinarySearch(const unsigned int *src, const unsigned int begin, const unsigned int end, const unsigned int target);

#if SI == 0
    static const unsigned int BinarySearchForGallopingSearchAVX2(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val);

    static const unsigned int GallopingSearchAVX2(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val);
#elif SI == 1
    static const unsigned int BinarySearchForGallopingSearchAVX512(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val);

    static const unsigned int GallopingSearchAVX512(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val);
#elif SI == 2
#endif

};

#endif //FSM_UTILITY_H
