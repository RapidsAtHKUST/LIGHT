//
// Created by ssunah on 11/11/17.
//

#ifndef FSE_FSE_H
#define FSE_FSE_H

#include "Graph.h"
#include "TaskPool.h"
#include <pthread.h>

class FSE;

struct WorkerParam {
    FSE* fse;
    unsigned int logic_thread_id;
};

class FSE {
private:
    const Graph* data_graph_;

    size_t embedding_count_;

    // Parallel usage.
    pthread_mutex_t mtx_;
    pthread_cond_t cv_;

    unsigned int task_queue_capacity_;
    unsigned int idle_thread_num_;
    unsigned int thread_num_;

    TaskPool* task_queue_;

    // Buffer for each thread.
    size_t *p_embedding_count_;

    // Parallel execution statistics.
    double* p_thread_execution_time_;

    // Temporary for k-clique search.
    unsigned int k_clique_;

    // Is Plus.
    bool plus_;

    void (FSE::*process_)(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                         unsigned int *cache, TaskSlot &task, size_t &embedding_count);

public:
    double search_time_;
    double total_time_;

private:
    void Initialize();
    void Release();

    static void ThreadInitialize(FSE* fse, unsigned int *&idx, unsigned int *&cn_count, unsigned int **&cn, unsigned int *&offset, unsigned int *&cache);
    static void ThreadRelease(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset, unsigned int *cache);

private:
    size_t set_intersection_invoke_cnt_;

    void PrintStatistics();

private:
    // Parallel Execution.
    void ExecuteParallel();
    void FillTaskPool(const bool plus);
    void SplitTask(const bool plus, const unsigned int vertex_index, const unsigned int edge_index, TaskSlot &task, unsigned int &edge_end);
    static void* Execute(void* input);

private:
    void ParallelProcessSquare(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                                   unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessHouse(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                              unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessKClique(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                                unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessChordSquare(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                                    unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessNear5Clique(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                                    unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessTripleTriangles(unsigned int *idx, unsigned int *cn_count, unsigned int **cn, unsigned int *offset,
                                        unsigned int *cache, TaskSlot &task, size_t &embedding_count);

    void ParallelProcessQuadrupleTriangles(unsigned int *idx, unsigned int *cn_count, unsigned int **cn,
                                              unsigned int *offset,
                                              unsigned int *cache, TaskSlot &task, size_t &embedding_count);
public:
    size_t Enumerate(const Graph *data_graph, const string &pattern_name, const unsigned int thread_num);

    ~FSE();
};


#endif //FSE_FSE_H
