//
// Created by ssunah on 11/30/17.
//

#include "ComputeSetIntersection.h"

void ComputeSetIntersection::ComputeCandidates(const unsigned int *larray, const unsigned int l_count,
                                               const unsigned int *rarray, const unsigned int r_count,
                                               unsigned int *cn, unsigned int &cn_count) {
#if HYBRID == 0
    #if SI == 0
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
    #elif SI == 1
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGallopingAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    }
    #elif SI == 2
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGalloping(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn, cn_count);
    }
    #endif
#elif HYBRID == 1
    #if SI == 0
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    #elif SI == 1
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    #elif SI == 2
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn, cn_count);
    #endif
#endif
}

void ComputeSetIntersection::ComputeCandidates(const unsigned int *larray, const unsigned int l_count,
                                               const unsigned int *rarray, const unsigned int r_count,
                                               unsigned int &cn_count) {
#if HYBRID == 0
    #if SI == 0
        if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn_count);
        }
    #elif SI == 1
        if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGallopingAVX512(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn_count);
        }
    #elif SI == 2
        if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGalloping(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn_count);
        }
    #endif
#elif HYBRID == 1
    #if SI == 0
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn_count);
    #elif SI == 1
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn_count);
    #elif SI == 2
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn_count);
    #endif
#endif
}

void ComputeSetIntersection::ComputeCandidatesNoCount(const unsigned int *larray, const unsigned int l_count,
                                               const unsigned int *rarray, const unsigned int r_count,
                                               unsigned int *cn, unsigned int &cn_count) {
#if HYBRID == 0
#if SI == 0
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    }
#elif SI == 1
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGallopingAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    }
#elif SI == 2
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        return ComputeCNGalloping(larray, l_count, rarray, r_count, cn, cn_count);
    }
    else {
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn, cn_count);
    }
#endif
#elif HYBRID == 1
    #if SI == 0
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn, cn_count);
    #elif SI == 1
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn, cn_count);
    #elif SI == 2
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn, cn_count);
    #endif
#endif
}

void ComputeSetIntersection::ComputeCandidatesNoCount(const unsigned int *larray, const unsigned int l_count,
                                               const unsigned int *rarray, const unsigned int r_count,
                                               unsigned int &cn_count) {
#if HYBRID == 0
#if SI == 0
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGallopingAVX2(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn_count);
        }
#elif SI == 1
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
            return ComputeCNGallopingAVX512(larray, l_count, rarray, r_count, cn_count);
        }
        else {
            return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn_count);
        }
#elif SI == 2
    if (l_count / 50 > r_count || r_count / 50 > l_count) {
        galloping_cnt_ += 1;
        return ComputeCNGalloping(larray, l_count, rarray, r_count, cn_count);
    }
    else {
        merge_cnt_ += 1;
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn_count);
    }
#endif
#elif HYBRID == 1
    #if SI == 0
        return ComputeCNMergeBasedAVX2(larray, l_count, rarray, r_count, cn_count);
    #elif SI == 1
        return ComputeCNMergeBasedAVX512(larray, l_count, rarray, r_count, cn_count);
    #elif SI == 2
        return ComputeCNNaiveStdMerge(larray, l_count, rarray, r_count, cn_count);
    #endif
#endif
}


#if SI == 0
void ComputeSetIntersection::ComputeCNGallopingAVX2(const unsigned int *larray, const unsigned int l_count,
                                                    const unsigned int *rarray, const unsigned int r_count,
                                                    unsigned int *cn, unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearchAVX2(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn[cn_count++] = larray[li];
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGallopingAVX2(const unsigned int *larray, const unsigned int l_count,
                                                    const unsigned int *rarray, const unsigned int r_count,
                                                    unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearchAVX2(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn_count += 1;
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX2(const unsigned int *larray, const unsigned int l_count,
                                                     const unsigned int *rarray, const unsigned int r_count,
                                                     unsigned int *cn, unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
    __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);
    unsigned int *cur_back_ptr = cn;

    auto size_ratio = (rc) / (lc);
    if (size_ratio > 2) {
        if (li < lc && ri + 7 < rc) {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask != 0) {
                    // at most 1 element
                    *cur_back_ptr = larray[li];
                    cur_back_ptr += 1;
                }
                if (larray[li] > rarray[ri + 7]) {
                    ri += 8;
                    if (ri + 7 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                } else {
                    li++;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        if (li + 1 < lc && ri + 3 < rc) {
            __m256i u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
            __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
            __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask << 16 != 0) {
                    *cur_back_ptr = larray[li];
                    cur_back_ptr += 1;
                }
                if (real_mask >> 16 != 0) {
                    *cur_back_ptr = larray[li + 1];
                    cur_back_ptr += 1;
                }


                if (larray[li + 1] == rarray[ri + 3]) {
                    li += 2;
                    ri += 4;
                    if (li + 1 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else if (larray[li + 1] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else {
                    li += 2;
                    if (li + 1 >= lc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                }
            }
        }
    }

    cn_count = (unsigned int)(cur_back_ptr - cn);
    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                ++li;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ++ri;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                // write back
                cn[cn_count++] = larray[li];

                ++li;
                ++ri;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX2(const unsigned int *larray, const unsigned int l_count,
                                                     const unsigned int *rarray, const unsigned int r_count,
                                                     unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    constexpr int parallelism = 8;

    int cn_countv[parallelism] = {0, 0, 0, 0, 0, 0, 0, 0};
    __m256i sse_cn_countv = _mm256_load_si256((__m256i *) (cn_countv));
    __m256i sse_countplus = _mm256_set1_epi32(1);
    __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
    __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);

    auto size_ratio = (rc) / (lc);
    if (size_ratio > 2) {
        if (li < lc && ri + 7 < rc) {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                mask = _mm256_and_si256(sse_countplus, mask);
                sse_cn_countv = _mm256_add_epi32(sse_cn_countv, mask);
                if (larray[li] > rarray[ri + 7]) {
                    ri += 8;
                    if (ri + 7 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                } else {
                    li++;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
            _mm256_store_si256((__m256i *) cn_countv, sse_cn_countv);
            for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
        }
    } else {
        if (li + 1 < lc && ri + 3 < rc) {
            __m256i u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
            __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
            __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                mask = _mm256_and_si256(sse_countplus, mask);
                sse_cn_countv = _mm256_add_epi32(sse_cn_countv, mask);

                if (larray[li + 1] == rarray[ri + 3]) {
                    li += 2;
                    ri += 4;
                    if (li + 1 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else if (larray[li + 1] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else {
                    li += 2;
                    if (li + 1 >= lc) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                }
            }
        }
        _mm256_store_si256((__m256i *) cn_countv, sse_cn_countv);
        for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
    }

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                ++li;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ++ri;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                cn_count++;
                ++li;
                ++ri;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}
#elif SI == 1
void ComputeSetIntersection::ComputeCNGallopingAVX512(const unsigned int *larray, const unsigned int l_count,
                                                          const unsigned int *rarray, const unsigned int r_count,
                                                          unsigned int *cn, unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearchAVX512(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn[cn_count++] = larray[li];
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGallopingAVX512(const unsigned int *larray, const unsigned int l_count,
                                                          const unsigned int *rarray, const unsigned int r_count,
                                                          unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearchAVX512(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn_count += 1;
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX512(const unsigned int *larray, const unsigned int l_count,
                                                       const unsigned int *rarray, const unsigned int r_count,
                                                       unsigned int *cn, unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);

    unsigned int *cur_back_ptr = cn;

    auto size1 = (rc) / (lc);
    if (size1 > 2) {
        if (li < lc && ri + 15 < rc) {
            __m512i u_elements = _mm512_set1_epi32(larray[li]);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements, v_elements);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li] > rarray[ri + 15]) {
                    ri += 16;
                    if (ri + 15 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                } else {
                    li += 1;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm512_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        if (li + 3 < lc && ri + 3 < rc) {
            __m512i u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
            __m512i u_elements_per = _mm512_permutevar_epi32(st, u_elements);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
            __m512i v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements_per, v_elements_per);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements_per);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li + 3] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                    v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);
                } else if (larray[li + 3] < rarray[ri + 3]) {
                    li += 4;
                    if (li + 3 >= lc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
                    u_elements_per = _mm512_permutevar_epi32(st, u_elements);
                } else {
                    li += 4;
                    ri += 4;
                    if (li + 3 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
                    u_elements_per = _mm512_permutevar_epi32(st, u_elements);
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                    v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);
                }
            }
        }
    }

    cn_count = (unsigned int)(cur_back_ptr - cn);

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                li += 1;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ri += 1;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                // write back
                cn[cn_count++] = larray[li];

                li += 1;
                ri += 1;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}

void ComputeSetIntersection::ComputeCNMergeBasedAVX512(const unsigned int *larray, const unsigned int l_count,
                                                       const unsigned int *rarray, const unsigned int r_count,
                                                       unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    constexpr int parallelism = 16;
    __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
    __m512i ssecountplus = _mm512_set1_epi32(1);
    int cn_countv[parallelism] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    __m512i ssecn_countv = _mm512_set1_epi32(0);
    auto size1 = (rc) / (lc);

    if (size1 > 2) {
        if (li < lc && ri + 15 < rc) {
            __m512i u_elements = _mm512_set1_epi32(larray[li]);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements, v_elements);
                ssecn_countv = _mm512_mask_add_epi32(ssecn_countv, mask, ssecn_countv, ssecountplus);

                if (larray[li] > rarray[ri + 15]) {
                    ri += 16;
                    if (ri + 15 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                } else {
                    li += 1;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm512_set1_epi32(larray[li]);
                }
            }
            _mm512_storeu_si512((__m512i *) cn_countv, ssecn_countv);
            for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
        }
    } else {
        if (li + 3 < lc && ri + 3 < rc) {
            __m512i u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
            __m512i u_elements_per = _mm512_permutevar_epi32(st, u_elements);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
            __m512i v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements_per, v_elements_per);
                ssecn_countv = _mm512_mask_add_epi32(ssecn_countv, mask, ssecn_countv, ssecountplus);

                if (larray[li + 3] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                    v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);
                } else if (larray[li + 3] < rarray[ri + 3]) {
                    li += 4;
                    if (li + 3 >= lc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
                    u_elements_per = _mm512_permutevar_epi32(st, u_elements);
                } else {
                    li += 4;
                    ri += 4;
                    if (li + 3 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray + li));
                    u_elements_per = _mm512_permutevar_epi32(st, u_elements);
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray + ri));
                    v_elements_per = _mm512_permute4f128_epi32(v_elements, 0b00000000);
                }
            }
            _mm512_storeu_si512((__m512i *) cn_countv, ssecn_countv);
            for (int cn_countvplus : cn_countv) { cn_count += cn_countvplus; }
        }
    }

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                li += 1;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ri += 1;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                cn_count += 1;
                li += 1;
                ri += 1;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
}

#elif SI == 2
void ComputeSetIntersection::ComputeCNNaiveStdMerge(const unsigned int *larray, const unsigned int l_count,
                                                    const unsigned int *rarray, const unsigned int r_count,
                                                    unsigned int *cn, unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        if (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }
        else if (larray[li] > rarray[ri]) {
            ri += 1;
            if (ri >= rc) {
                return;
            }
        }
        else {
            cn[cn_count++] = larray[li];

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNNaiveStdMerge(const unsigned int *larray, const unsigned int l_count,
                                                    const unsigned int *rarray, const unsigned int r_count,
                                                    unsigned int &cn_count) {
    cn_count = 0;

    if (l_count == 0 || r_count == 0)
        return;

    unsigned int lc = l_count;
    unsigned int rc = r_count;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        if (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }
        else if (larray[li] > rarray[ri]) {
            ri += 1;
            if (ri >= rc) {
                return;
            }
        }
        else {
            cn_count += 1;
            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGalloping(const unsigned int *larray, const unsigned int l_count,
                                                const unsigned int *rarray, const unsigned int r_count,
                                                unsigned int *cn, unsigned int &cn_count) {
    unsigned int lc = l_count;
    unsigned int rc = r_count;
    cn_count = 0;
    if (lc == 0 || rc == 0)
        return;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearch(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn[cn_count++] = larray[li];

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}

void ComputeSetIntersection::ComputeCNGalloping(const unsigned int *larray, const unsigned int l_count,
                                                const unsigned int *rarray, const unsigned int r_count,
                                                unsigned int &cn_count) {
    unsigned int lc = l_count;
    unsigned int rc = r_count;
    cn_count = 0;
    if (lc == 0 || rc == 0)
        return;

    if (lc > rc) {
        auto tmp = larray;
        larray = rarray;
        rarray = tmp;

        unsigned int tmp_count = lc;
        lc = rc;
        rc = tmp_count;
    }

    unsigned int li = 0;
    unsigned int ri = 0;

    while (true) {
        while (larray[li] < rarray[ri]) {
            li += 1;
            if (li >= lc) {
                return;
            }
        }

        ri = Utility::GallopingSearch(rarray, ri, rc, larray[li]);
        if (ri >= rc) {
            return;
        }

        if (larray[li] == rarray[ri]) {
            cn_count += 1;

            li += 1;
            ri += 1;
            if (li >= lc || ri >= rc) {
                return;
            }
        }
    }
}
#endif