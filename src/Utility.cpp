//
// Created by ssunah on 11/30/16.
//

#include "Utility.h"

timeval Utility::GetTime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp;
}

double Utility::TimeDiffInSeconds(const timeval start, const timeval end)
{
    double second_diff = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / (double)1000000.0;
    return second_diff;
}

void Utility::FindAutomorphism(const Graph *graph, const unsigned int* enumerate_order, const unsigned int* pivot,
                               const unsigned int enumerate_length, vector<vector<unsigned int> > &aut) {
    bool* visited = new bool[enumerate_length];
    fill(visited, visited + enumerate_length, false);
    unsigned int* embedding = new unsigned int[enumerate_length];
    unsigned int* mapping = new unsigned int[enumerate_length];
    unsigned int* candidate_idx = new unsigned int[enumerate_length];
    unsigned int* candidate_count = new unsigned int[enumerate_length];
    const unsigned int** candidate = new const unsigned int*[enumerate_length];

    for (unsigned int i = 0; i < enumerate_length; ++i) {
        unsigned int u = enumerate_order[0];
        unsigned int v = i;

        if (graph->Degree(u) <= graph->Degree(v)) {
            embedding[0] = v;
            mapping[u] = v;
            visited[v] = true;

            unsigned int cur_level = 1;
            candidate_idx[cur_level] = 0;
            candidate[cur_level] = graph->Neighbors(v, candidate_count[cur_level]);

            while (true) {
                while(candidate_idx[cur_level] < candidate_count[cur_level]) {
                    u = enumerate_order[cur_level];
                    v = candidate[cur_level][candidate_idx[cur_level]];
                    candidate_idx[cur_level] += 1;

                    // Check whether u and v is feasible.
                    bool is_feasible = false;
                    if (!visited[v] && graph->Degree(u) <= graph->Degree(v)) {
                        bool valid = true;
                        for (unsigned int j = 0; j < cur_level; ++j) {
                            unsigned int before_u = enumerate_order[j];

                            if (graph->IsEdge(before_u, u)) {
                                if (!graph->IsEdge(mapping[before_u], v)) {
                                    valid = false;
                                    break;
                                }
                            }
                        }

                        if (valid)
                            is_feasible = true;
                    }

                    if (is_feasible) {
                        if (cur_level == enumerate_length - 1) {
                            // Find an embedding
                            mapping[u] = v;
                            aut.push_back(vector<unsigned int> (mapping, mapping + enumerate_length));
                        }
                        else {
                            embedding[cur_level] = v;
                            mapping[u] = v;
                            visited[v] = true;

                            cur_level += 1;
                            candidate_idx[cur_level] = 0;
                            candidate[cur_level] = graph->Neighbors(mapping[pivot[cur_level]], candidate_count[cur_level]);
                        }
                    }
                }

                cur_level -= 1;
                visited[embedding[cur_level]] = false;

                if (cur_level == 0) {
                    break;
                }
            }
        }
    }

    delete[] visited;
    delete[] embedding;
    delete[] mapping;
    delete[] candidate;
    delete[] candidate_count;
    delete[] candidate_idx;
}

const unsigned int Utility::FirstIndexGreaterThan(const unsigned int *src, const unsigned int begin, const unsigned int end, const unsigned int target) {
    int offset_begin = begin;
    int offset_end = (int)end - 1;
    while (offset_begin <= offset_end) {
        int mid = offset_begin + ((offset_end - offset_begin) >> 1);
        _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
        if (src[mid] == target) {
            return mid + 1;
        } else if (src[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid - 1;
        }
    }

    return (unsigned int)offset_begin;
}

bool Utility::Contain(const unsigned int *src, const unsigned int begin, unsigned int end, const unsigned int target) {
    int offset_begin = begin;
    int offset_end = (int)end;
    while (offset_end - offset_begin >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_begin) + offset_end) / 2);
        _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
        if (src[mid] == target) {
            return true;
        } else if (src[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (src[offset] == target)
            return true;
        else if (src[offset] > target)
            return false;
    }

    return false;
}

const unsigned int Utility::GallopingSearch(const unsigned int *src, const unsigned int begin, const unsigned int end,
                                            const unsigned int target) {
    if (src[end - 1] < target) {
        return end;
    }
    // galloping
    if (src[begin] >= target) {
        return begin;
    }
    if (src[begin + 1] >= target) {
        return begin + 1;
    }
    if (src[begin + 2] >= target) {
        return begin + 2;
    }

    unsigned int jump_idx = 4;
    unsigned int offset_beg = begin;
    while (true) {
        unsigned int peek_idx = offset_beg + jump_idx;
        if (peek_idx >= end) {
            return BinarySearch(src, (jump_idx >> 1) + offset_beg + 1, end, target);
        }
        if (src[peek_idx] < target) {
            jump_idx <<= 1;
        } else {
            return src[peek_idx] == target ? peek_idx :
                   BinarySearch(src, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, target);
        }
    }
}

const unsigned int Utility::BinarySearch(const unsigned int *src, const unsigned int begin, const unsigned int end,
                                         const unsigned int target) {
    int offset_begin = begin;
    int offset_end = end;
    while (offset_end - offset_begin >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_begin) + offset_end) / 2);
        _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
        if (src[mid] == target) {
            return mid;
        } else if (src[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (src[offset] >= target) {
            return (unsigned int)offset;
        }
    }

    return (unsigned int)offset_end;
}

void Utility::GetKCore(const Graph *graph, unsigned int *core_table) {
    unsigned int vertices_num = graph->VerticesCount();
    unsigned int max_degree = graph->MaxDegree();

    unsigned int* vertices = new unsigned int[vertices_num];          // Vertices sorted by degree.
    unsigned int* position = new unsigned int[vertices_num];          // The position of vertices in vertices array.
    unsigned int* degree_bin = new unsigned int[max_degree + 1];      // Degree from 0 to max_degree.
    unsigned int* offset = new unsigned int[max_degree + 1];          // The offset in vertices array according to degree.

    memset(degree_bin, 0, sizeof(unsigned int) * (max_degree + 1));

    for (unsigned int i = 0; i < vertices_num; ++i) {
        unsigned int degree = graph->Degree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    unsigned int start = 0;
    for (unsigned int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (unsigned int i = 0; i < vertices_num; ++i) {
        unsigned int degree = graph->Degree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (unsigned int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (unsigned int i = 0; i < vertices_num; ++i) {
        unsigned int v = vertices[i];

        unsigned int count;
        const unsigned int* neighbors = graph->Neighbors(v, count);

        for(unsigned int j = 0; j < count; ++j) {
            unsigned int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                unsigned int cur_degree_u = core_table[u];
                unsigned int position_u = position[u];
                unsigned int position_w = offset[cur_degree_u];
                unsigned int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }
    }

    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}

#if SI == 0
const unsigned int Utility::BinarySearchForGallopingSearchAVX2(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val) {
    while (offset_end - offset_beg >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
        if (array[mid] == val) {
            return mid;
        } else if (array[mid] < val) {
            offset_beg = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback, be careful with operator>> and operation+ priority
    __m256i pivot_element = _mm256_set1_epi32(val);
    for (; offset_beg + 7 < offset_end; offset_beg += 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) {
            return offset_beg + (_popcnt32(mask) >> 2);
        }
    }
    if (offset_beg < offset_end) {
        auto left_size = offset_end - offset_beg;
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    }
    return offset_end;
}

const unsigned int Utility::GallopingSearchAVX2(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val) {
    if (array[offset_end - 1] < val) {
        return offset_end;
    }

    // linear search
    __m256i pivot_element = _mm256_set1_epi32(val);
    if (offset_end - offset_beg >= 8) {
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        if (mask != 0xffffffff) { return offset_beg + (_popcnt32(mask) >> 2); }
    } else {
        auto left_size = offset_end - offset_beg;
        __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
        __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
        int mask = _mm256_movemask_epi8(cmp_res);
        int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + (_popcnt32(mask) >> 2); }
    }

    // galloping, should add pre-fetch later
    auto jump_idx = 8u;
    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}
#elif SI == 1
const unsigned int Utility::BinarySearchForGallopingSearchAVX512(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val) {
    while (offset_end - offset_beg >= 32) {
        auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
        if (array[mid] == val) {
            return mid;
        } else if (array[mid] < val) {
            offset_beg = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    constexpr int parallelism = 16;
    __m512i pivot_element = _mm512_set1_epi32(val);
    for (; offset_beg + 15 < offset_end; offset_beg += parallelism) {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        if (mask != 0xffff) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    if (offset_beg < offset_end) {
        auto left_size = offset_end - offset_beg;
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        __mmask16 cmp_mask = ((__mmask16) 0xffff) >> (16 - left_size);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + _mm_popcnt_u32(mask); }
    }
    return offset_end;
}

const unsigned int Utility::GallopingSearchAVX512(const unsigned int *array, unsigned int offset_beg, unsigned int offset_end, unsigned int val) {
    if (array[offset_end - 1] < val) {
        return offset_end;
    }

    // front peeking
    __m512i pivot_element = _mm512_set1_epi32(val);
    auto left_size = offset_end - offset_beg;
    if (left_size >= 16) {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        if (mask != 0xffff) { return offset_beg + _mm_popcnt_u32(mask); }
    } else {
        __m512i elements = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(array + offset_beg));
        __mmask16 mask = _mm512_cmplt_epi32_mask(elements, pivot_element);
        __mmask16 cmp_mask = ((__mmask16) 0xffff) >> (16 - left_size);
        mask &= cmp_mask;
        if (mask != cmp_mask) { return offset_beg + _mm_popcnt_u32(mask); }
    }

    // galloping
    auto jump_idx = 16u;
    // pre-fetch
    auto jump_times = 32 - _lzcnt_u32((offset_end - offset_beg) >> 4);
    __m512i prefetch_idx = _mm512_set_epi64(16, 32, 64, 128, 256, 512, 1024, 2048);
    __mmask8 mask = jump_times >= 8 ? (__mmask8) 0xff : (__mmask8) 0xff << (8 - jump_times);
    _mm512_mask_prefetch_i64gather_ps(prefetch_idx, mask, (void *)(array + offset_beg), 1, _MM_HINT_T0);

    while (true) {
        auto peek_idx = offset_beg + jump_idx;
        if (peek_idx >= offset_end) {
            return BinarySearchForGallopingSearchAVX512(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
        }
        if (array[peek_idx] < val) {
            jump_idx <<= 1;
        } else {
            return array[peek_idx] == val ? peek_idx :
                   BinarySearchForGallopingSearchAVX512(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
        }
    }
}
#elif SI == 2
#endif