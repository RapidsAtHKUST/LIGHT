//
// Created by ssunah on 11/30/17.
//

#ifndef FSE_COMPUTESETINTERSECTION_H
#define FSE_COMPUTESETINTERSECTION_H

#include "Utility.h"

/*
 * Because the set intersection is designed for computing common neighbors, the target is unsigned integer.
 */

class ComputeSetIntersection {
public:
#if HYBRID == 0
    static size_t galloping_cnt_;
    static size_t merge_cnt_;
#endif
    static void ComputeCandidates(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                  const unsigned int r_count, unsigned int *cn, unsigned int &cn_count);
    static void ComputeCandidates(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                  const unsigned int r_count, unsigned int &cn_count);

    static void ComputeCandidatesNoCount(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                  const unsigned int r_count, unsigned int *cn, unsigned int &cn_count);
    static void ComputeCandidatesNoCount(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                  const unsigned int r_count, unsigned int &cn_count);
#if SI == 0
    static void ComputeCNGallopingAVX2(const unsigned int *larray, const unsigned int l_count,
                                       const unsigned int *rarray, const unsigned int r_count, unsigned int *cn,
                                       unsigned int &cn_count);
    static void ComputeCNGallopingAVX2(const unsigned int *larray, const unsigned int l_count,
                                       const unsigned int *rarray, const unsigned int r_count, unsigned int &cn_count);

    static void ComputeCNMergeBasedAVX2(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                        const unsigned int r_count, unsigned int *cn, unsigned int &cn_count);
    static void ComputeCNMergeBasedAVX2(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                        const unsigned int r_count, unsigned int &cn_count);
#elif SI == 1

    static void ComputeCNGallopingAVX512(const unsigned int *larray, const unsigned int l_count,
                                         const unsigned int *rarray, const unsigned int r_count, unsigned int *cn,
                                         unsigned int &cn_count);
    static void ComputeCNGallopingAVX512(const unsigned int *larray, const unsigned int l_count,
                                         const unsigned int *rarray, const unsigned int r_count, unsigned int &cn_count);

    static void ComputeCNMergeBasedAVX512(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                          const unsigned int r_count, unsigned int *cn, unsigned int &cn_count);
    static void ComputeCNMergeBasedAVX512(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                          const unsigned int r_count, unsigned int &cn_count);

#elif SI == 2

    static void ComputeCNNaiveStdMerge(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                       const unsigned int r_count, unsigned int *cn, unsigned int &cn_count);
    static void ComputeCNNaiveStdMerge(const unsigned int *larray, const unsigned int l_count, const unsigned int *rarray,
                                       const unsigned int r_count, unsigned int &cn_count);

    static void ComputeCNGalloping(const unsigned int* larray, const unsigned int l_count, const unsigned int* rarray,
                                   const unsigned int r_count, unsigned int* cn, unsigned int& cn_count);
    static void ComputeCNGalloping(const unsigned int* larray, const unsigned int l_count, const unsigned int* rarray,
                                   const unsigned int r_count, unsigned int& cn_count);

#endif
};


#endif //FSE_COMPUTESETINTERSECTION_H
