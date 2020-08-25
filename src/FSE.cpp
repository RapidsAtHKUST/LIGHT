//
// Created by ssunah on 11/11/17.
//

#include "FSE.h"
#include "Config.h"
#include "Utility.h"
#include "ComputeSetIntersection.h"
#include <set>
#include <iomanip>
#include <queue>

size_t FSE::Enumerate(const Graph *data_graph, const string &pattern_name, const ui thread_num) {
    data_graph_ = data_graph;
    thread_num_ = thread_num;

    Initialize();

    if (pattern_name == "p0") {
        plus_ = true;
        k_clique_ = 3;
        process_ = &FSE::ParallelProcessKClique;
    }
    else if (pattern_name == "p1") {
        plus_ = true;
        process_ = &FSE::ParallelProcessSquare;
    }
    else if (pattern_name == "p2") {
        plus_ = true;
        process_ = &FSE::ParallelProcessChordSquare;
    }
    else if (pattern_name == "p3") {
        plus_ = true;
        k_clique_ = 4;
        process_ = &FSE::ParallelProcessKClique;
    }
    else if (pattern_name == "p4") {
        plus_ = true;
        process_ = &FSE::ParallelProcessHouse;
    }
    else if (pattern_name == "p5") {
        plus_ = false;
        process_ = &FSE::ParallelProcessQuadrupleTriangles;
    }
    else if (pattern_name == "p6") {
        plus_ = true;
        process_ = &FSE::ParallelProcessNear5Clique;
    }
    else if (pattern_name == "p7") {
        plus_ = true;
        k_clique_ = 5;
        process_ = &FSE::ParallelProcessKClique;
    }
    else if (pattern_name == "p8") {
        plus_ = false;
        process_ = &FSE::ParallelProcessTripleTriangles;
    }
    else {
        cout << "This pattern is not supported yet." << endl;
        exit(-1);
    }

    FillTaskPool(plus_);

    timeval search_start = Utility::GetTime();

    ExecuteParallel();

    timeval search_end = Utility::GetTime();

    PrintStatistics();

    search_time_ = Utility::TimeDiffInSeconds(search_start, search_end);
    total_time_ = search_time_;

    Release();

    return embedding_count_;
}

void FSE::Initialize() {
    // Parallel usage.
    embedding_count_ = 0;
    idle_thread_num_ = 0;
    task_queue_capacity_ = thread_num_ * 4;
    task_queue_ = new TaskPool(task_queue_capacity_);

    p_embedding_count_ = new size_t[thread_num_];
    fill(p_embedding_count_, p_embedding_count_ + thread_num_, 0);

    p_thread_execution_time_ = new double[thread_num_];

    set_intersection_invoke_cnt_ = 0;
    pthread_mutex_init(&mtx_, NULL);
    pthread_cond_init(&cv_, NULL);
}

void FSE::Release() {
    delete task_queue_;
    delete[] p_embedding_count_;
    delete[] p_thread_execution_time_;
    pthread_mutex_destroy(&mtx_);
    pthread_cond_destroy(&cv_);
}

void FSE::ThreadInitialize(FSE* fse, ui *&idx, ui *&cn_count, ui **&cn, ui *&offset, ui *&cache) {
    idx = new ui[MAXIMUM_SIZE];
    cn_count = new ui[MAXIMUM_SIZE];
    cn = new ui*[MAXIMUM_SIZE];

    for (ui j = 0; j < MAXIMUM_SIZE; ++j) {
        cn[j] = new ui[fse->data_graph_->MaxDegree()];
    }

    offset = new ui[fse->data_graph_->MaxDegree() + 1];
    cache = new ui[fse->data_graph_->MaxDegreeSum()];
}

void FSE::ThreadRelease(ui *idx, ui *cn_count, ui **cn, ui *offset, ui *cache) {
    delete[] idx;
    delete[] cn_count;

    for (ui i = 0; i < MAXIMUM_SIZE; ++i) {
        delete[] cn[i];
    }

    delete[] cn;
    delete[] offset;
    delete[] cache;
}

FSE::~FSE() {
}

void FSE::ExecuteParallel() {
    WorkerParam* params = new WorkerParam[thread_num_];
    pthread_t* threads = new pthread_t[thread_num_];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (ui i = 0; i < thread_num_; ++i) {
        params[i].fse = this;
        params[i].logic_thread_id = i;
        pthread_create(&threads[i], &attr, Execute, (void*)&params[i]);
    }

    for (ui i = 0; i < thread_num_; ++i) {
        pthread_join(threads[i], NULL);
    }

    for (ui i = 0; i < thread_num_; ++i) {
        embedding_count_ += p_embedding_count_[i];
    }

    pthread_attr_destroy(&attr);

    delete[] threads;
    delete[] params;

}

void *FSE::Execute(void *input) {
    FSE* fse = ((WorkerParam*)input)->fse;
    ui logic_thread_id = ((WorkerParam*)input)->logic_thread_id;

    size_t local_embedding_count = 0;
    double local_time = 0;
    TaskSlot task;

    // Initialize the local buffer of each thread.
    ui *idx;
    ui *cn_count;
    ui **cn;
    ui *offset;
    ui *cache;

    ThreadInitialize(fse, idx, cn_count, cn, offset, cache);

    while (true) {
        if (fse->task_queue_->Pop(task)) {
            timeval start = Utility::GetTime();

            (fse->*fse->process_)(idx, cn_count, cn, offset, cache, task, local_embedding_count);

            timeval end = Utility::GetTime();
            local_time += Utility::TimeDiffInSeconds(start, end);
        }
        else {
            // Acquire mutex to update the number of idle thread
            pthread_mutex_lock(&fse->mtx_);
            fse->idle_thread_num_ += 1;

            if (fse->idle_thread_num_ == fse->thread_num_) {
                // Wake up all other threads
                pthread_cond_broadcast(&fse->cv_);

                // Unlock mutex and exit
                pthread_mutex_unlock(&fse->mtx_);

                fse->p_embedding_count_[logic_thread_id] = local_embedding_count;
                fse->p_thread_execution_time_[logic_thread_id] = local_time;

                // Release the local buffer of this thread.
                ThreadRelease(idx, cn_count, cn, offset, cache);
                pthread_exit(NULL);
            }
            else {
                // Wait for signal
                while (pthread_cond_wait(&fse->cv_, &fse->mtx_) != 0);

                if (fse->idle_thread_num_ == fse->thread_num_) {
                    // All threads have been finished, unlock mutex and exit
                    pthread_mutex_unlock(&fse->mtx_);

                    fse->p_embedding_count_[logic_thread_id] = local_embedding_count;
                    fse->p_thread_execution_time_[logic_thread_id] = local_time;

                    // Release the local buffer of this thread.
                    ThreadRelease(idx, cn_count, cn, offset, cache);
                    pthread_exit(NULL);
                }
                else {
                    fse->idle_thread_num_ -= 1;
                    pthread_mutex_unlock(&fse->mtx_);
                }
            }
        }
    }
}

void FSE::FillTaskPool(const bool plus) {
    ui candidate_edge_count = (plus ? data_graph_->MaxNeighborPlusSum() : data_graph_->EdgesCount());

    ui step_length = candidate_edge_count / task_queue_capacity_;
    ui step_remainder = candidate_edge_count % task_queue_capacity_;
    ui first_task_step_length = step_length + step_remainder;

    TaskSlot task;
    task.edge_offset_begin_ = 0;
    task.edge_begin_ = 0;

    ui temp_sum = 0;
    ui cur_vertex_id = 0;
    ui cur_edge_index = 0;
    ui cur_vertex_degree = plus ? data_graph_->DegreePlus(cur_vertex_id) : data_graph_->Degree(cur_vertex_id);

    for (ui i = 0; i < task_queue_capacity_; ++i) {
        // Add the remainder to the first task slot.
        ui cur_step_length = i == 0 ? first_task_step_length : step_length;

        while (true) {
            ui cur_edge_count = cur_vertex_degree - cur_edge_index;
            if (temp_sum + cur_edge_count >= cur_step_length) {
                break;
            }
            else {
                temp_sum += cur_edge_count;
                cur_vertex_id += 1;
                cur_edge_index = 0;
                cur_vertex_degree = plus ? data_graph_->DegreePlus(cur_vertex_id) : data_graph_->Degree(cur_vertex_id);
            }
        }

        cur_edge_index += cur_step_length - temp_sum;
        task.edge_offset_end_ = cur_vertex_id;
        task.edge_end_ = cur_edge_index;
        task.edge_size_ = cur_step_length;
        task_queue_->Push(task);

        temp_sum = 0;
        task.edge_offset_begin_ = cur_vertex_id;
        task.edge_begin_ = cur_edge_index;
    }
}

void FSE::SplitTask(const bool plus, const ui vertex_index, const ui edge_index, TaskSlot &task,
                    ui &edge_end) {
    if (idle_thread_num_ != 0 && task_queue_->Empty()) {
        if (task.edge_size_ >= MINIMUM_TASK_LENGTH) {
            ui step_length = task.edge_size_ / 2;

            ui temp_sum = 0;
            ui cur_vertex_id = vertex_index;
            ui cur_edge_index = edge_index;
            ui cur_vertex_degree = plus ? data_graph_->DegreePlus(cur_vertex_id) : data_graph_->Degree(cur_vertex_id);

            while (true) {
                ui cur_edge_count = cur_vertex_degree - cur_edge_index;
                if (temp_sum + cur_edge_count >= step_length) {
                    break;
                }
                else {
                    temp_sum += cur_edge_count;
                    cur_vertex_id += 1;
                    cur_edge_index = 0;
                    cur_vertex_degree = plus ? data_graph_->DegreePlus(cur_vertex_id) : data_graph_->Degree(cur_vertex_id);
                }
            }

            cur_edge_index += step_length - temp_sum;
            TaskSlot new_task (cur_vertex_id, task.edge_offset_end_, cur_edge_index, task.edge_end_, task.edge_size_ - step_length);

            task.edge_offset_end_ = cur_vertex_id;
            task.edge_end_ = cur_edge_index;
            task.edge_size_ = step_length;

            if (cur_vertex_id == vertex_index) {
                edge_end = cur_edge_index;
            }

            task_queue_->Push(new_task);
            pthread_cond_signal(&cv_);
        }
    }
}

void FSE::ParallelProcessSquare(ui *idx, ui *cn_count, ui **cn, ui *offset,
                                ui *cache, TaskSlot &task, size_t &embedding_count) {
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i){
        ui v0 = i;

        ui v0_neighbors_count_plus;
        const ui* v0_neighbors_plus = data_graph_->NeighborsPlus(v0, v0_neighbors_count_plus);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_neighbors_count_plus;

        for (ui j = edge_begin; j < edge_end; ++j) {
            ui v1 = v0_neighbors_plus[j];

            // Update the edge size.
            task.edge_size_ -= 1;
            SplitTask(plus_, v0, j + 1, task, edge_end);

            ui v1_neighbors_count;
            const ui* v1_neighbors = data_graph_->Neighbors(v1, v1_neighbors_count);
            const ui start_index_v1 = Utility::FirstIndexGreaterThan(v1_neighbors, 0, v1_neighbors_count, v0);

            for (ui k = j + 1; k < v0_neighbors_count_plus; ++k) {
                ui v2 = v0_neighbors_plus[k];

                ui v2_neighbors_count;
                const ui* v2_neighbors = data_graph_->Neighbors(v2, v2_neighbors_count);
                const ui start_index_v2 = Utility::FirstIndexGreaterThan(v2_neighbors, 0, v2_neighbors_count, v0);

                ui temp;

                if (v1_neighbors_count > start_index_v1 && v2_neighbors_count > start_index_v2) {
                    ComputeSetIntersection::ComputeCandidates(v1_neighbors + start_index_v1, v1_neighbors_count - start_index_v1,
                                             v2_neighbors + start_index_v2, v2_neighbors_count - start_index_v2, temp);

                    embedding_count += temp;
                }
            }
        }
    }
}

void FSE::ParallelProcessHouse(ui *idx, ui *cn_count, ui **cn, ui *offset,
                               ui *cache, TaskSlot &task, size_t &embedding_count) {
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {
        ui v0 = i;

        ui v0_neighbors_count_plus;
        const ui* v0_neighbors_plus = data_graph_->NeighborsPlus(v0, v0_neighbors_count_plus);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_neighbors_count_plus;

        if (data_graph_->Degree(v0) >= 3) {
            ui v0_neighbors_count;
            const ui* v0_neighbors = data_graph_->Neighbors(v0, v0_neighbors_count);

            for (ui j = edge_begin; j < edge_end; ++j) {
                ui v1 = v0_neighbors_plus[j];

                task.edge_size_ -= 1;
                SplitTask(plus_, v0, j + 1, task, edge_end);

                if (data_graph_->Degree(v1) >= 3) {
                    ui v1_neighbors_count;
                    const ui* v1_neighbors = data_graph_->Neighbors(v1, v1_neighbors_count);

#ifdef COLLECT_RESULTS
                    set_intersection_invoke_cnt_ += 1;
#endif

                    ComputeSetIntersection::ComputeCandidates(v0_neighbors, v0_neighbors_count, v1_neighbors, v1_neighbors_count, cn[0], cn_count[0]);

                    if (cn_count[0] < 1)
                        continue;

                    for (ui l = 0; l < v0_neighbors_count; ++l) {
                        ui v2 = v0_neighbors[l];

                        if (v2 != v1) {
                            ui v2_neighbors_count;
                            const ui* v2_neighbors = data_graph_->Neighbors(v2, v2_neighbors_count);

#ifdef COLLECT_RESULTS
                            set_intersection_invoke_cnt_ += 1;
#endif

                            ComputeSetIntersection::ComputeCandidates(v2_neighbors, v2_neighbors_count, v1_neighbors, v1_neighbors_count, cn[1], cn_count[1]);
                            if (cn_count[1] >= 2) {
#ifdef ENUMERATE_RESULTS
                                for (ui m = 0; m < cn_count[0]; ++m) {
                                    ui v3 = cn[0][m];
                                    if (v3 != v2) {
                                        for (ui n = 0; n < cn_count[1]; ++n) {
                                            ui v4 = cn[1][n];
                                            if (v4 != v0 && v4 != v3)
                                                embedding_count += 1;
                                        }
                                    }
                                }
#else
                                unsigned int temp_count;
                                ComputeSetIntersection::ComputeCandidates(cn[0], cn_count[0], cn[1], cn_count[1], temp_count);

                                unsigned int temp = cn_count[0];
                                if (Utility::Contain(cn[0], 0, cn_count[0], v2)) {
                                    temp -= 1;
                                }

                                embedding_count += (cn_count[1] - 2) * temp_count + (cn_count[1] - 1) * (temp - temp_count);
#endif
                            }
                        }
                    }
                }
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }
}

void FSE::ParallelProcessKClique(ui *idx, ui *cn_count, ui **cn, ui *offset,
                                 ui *cache, TaskSlot &task, size_t &embedding_count) {
    ui k_clique = k_clique_;
    ui minimum_neighbors_count = k_clique - 1;
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {

        ui begin = i;
        ui begin_neighbors_count;
        const ui *begin_neighbors = data_graph_->NeighborsPlus(begin, begin_neighbors_count);

        ui edge_begin = begin == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = begin == task.edge_offset_end_ ? task.edge_end_ : begin_neighbors_count;

        if (begin_neighbors_count >= minimum_neighbors_count) {
            for (ui j = edge_begin; j < edge_end; ++j) {
                ui end = begin_neighbors[j];

                task.edge_size_ -= 1;
                SplitTask(plus_, begin, j + 1, task, edge_end);

                ui end_neighbors_count;
                const ui *end_neighbors = data_graph_->NeighborsPlus(end, end_neighbors_count);

                ComputeSetIntersection::ComputeCandidates(begin_neighbors, begin_neighbors_count,
                                                        end_neighbors, end_neighbors_count,
                                                        cn[2], cn_count[2]);

                if (k_clique == 3) {
                    embedding_count += cn_count[2];
                }
                else {
                    if (cn_count[2] >= minimum_neighbors_count - 1) {
                        ui cur_level = 2;
                        ui next_level = 3;
                        idx[cur_level] = 0;

                        while (true) {
                            while (idx[cur_level] < cn_count[cur_level]) {
                                ui v = cn[cur_level][idx[cur_level]];
                                idx[cur_level] += 1;

                                ui v_neighbors_count;
                                const ui *v_neighbors = data_graph_->NeighborsPlus(v, v_neighbors_count);

                                ComputeSetIntersection::ComputeCandidates(cn[cur_level], cn_count[cur_level],
                                                                        v_neighbors, v_neighbors_count,
                                                                        cn[next_level], cn_count[next_level]);

                                if (cn_count[next_level] >= minimum_neighbors_count - cur_level) {
                                    if (k_clique == cur_level + 2) {
                                        embedding_count += cn_count[next_level];
                                    }
                                    else {
                                        cur_level += 1;
                                        next_level += 1;
                                        idx[cur_level] = 0;
                                    }
                                }
                            }

                            cur_level -= 1;
                            next_level -= 1;
                            if (cur_level <= 1) {
                                break;
                            }
                        }
                    }
                }
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }
}

void FSE::ParallelProcessChordSquare(ui *idx, ui *cn_count, ui **cn, ui *offset,
                                     ui *cache, TaskSlot &task, size_t &embedding_count) {
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {
        ui v0 = i;

        ui v0_neighbors_count_plus;
        const ui* v0_neighbors_plus = data_graph_->NeighborsPlus(v0, v0_neighbors_count_plus);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_neighbors_count_plus;

        if (data_graph_->Degree(v0) >= 3) {

            ui v0_neighbors_count;
            unsigned const int* v0_neighbors = data_graph_->Neighbors(v0, v0_neighbors_count);

            for (ui j = edge_begin; j < edge_end; ++j) {
                ui v1 = v0_neighbors_plus[j];

                task.edge_size_ -= 1;
                SplitTask(plus_, v0, j + 1, task, edge_end);

                ui v1_neighbors_count;
                const ui* v1_neighbors = data_graph_->Neighbors(v1, v1_neighbors_count);


#ifdef COLLECT_RESULTS
                set_intersection_invoke_cnt_ += 1;
#endif
                ComputeSetIntersection::ComputeCandidates(v0_neighbors, v0_neighbors_count, v1_neighbors, v1_neighbors_count, cn[0], cn_count[0]);

                if (cn_count[0] < 1)
                    continue;

#ifdef ENUMERATE_RESULTS
                for (ui k = 0; k < cn_count[0]; ++k) {
                    ui v2 = cn[0][k];
                    const ui start_index_v3 = Utility::FirstIndexGreaterThan(cn[0], 0, cn_count[0], v2);
                    if (start_index_v3 < cn_count[0])
                        embedding_count += cn_count[0] - start_index_v3;
                }
#else
                ui temp_count = cn_count[0];
                if (temp_count > 1) {
                    embedding_count += (temp_count - 1) * temp_count / 2;
                }
#endif
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }
}

void FSE::ParallelProcessNear5Clique(ui *idx, ui *cn_count, ui **cn, ui *offset,
                                     ui *cache, TaskSlot &task, size_t &embedding_count) {
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {
        ui v0 = i;

        ui v0_nbr_cnt_plus;
        const ui* v0_nbrs_plus = data_graph_->NeighborsPlus(v0, v0_nbr_cnt_plus);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_nbr_cnt_plus;

        if (data_graph_->Degree(v0) >= 4) {
            ui v0_nbr_cnt;
            const ui* v0_nbrs = data_graph_->Neighbors(v0, v0_nbr_cnt);

            for (ui j = edge_begin; j < edge_end; ++j) {
                ui v1 = v0_nbrs_plus[j];

                task.edge_size_ -= 1;
                SplitTask(plus_, v0, j + 1, task, edge_end);

                if (data_graph_->Degree(v1) >= 4) {
                    ui v1_nbr_cnt;
                    const ui* v1_nbrs = data_graph_->Neighbors(v1, v1_nbr_cnt);

#ifdef COLLECT_RESULTS
                    set_intersection_invoke_cnt_ += 1;
#endif

                    ComputeSetIntersection::ComputeCandidates(v0_nbrs, v0_nbr_cnt, v1_nbrs, v1_nbr_cnt,
                                                            cn[0], cn_count[0]);

                    if (cn_count[0] >= 3) {
                        for (ui k = 0; k < cn_count[0]; ++k) {
                            ui v2 = cn[0][k];

                            ui v2_nbrs_cnt_plus;
                            const ui* v2_nbrs_plus = data_graph_->NeighborsPlus(v2, v2_nbrs_cnt_plus);

                            if (v2_nbrs_cnt_plus >= 1) {

#ifdef COLLECT_RESULTS
                                set_intersection_invoke_cnt_ += 1;
#endif

                                ComputeSetIntersection::ComputeCandidates(v2_nbrs_plus, v2_nbrs_cnt_plus, cn[0], cn_count[0], cn[1], cn_count[1]);

#ifdef ENUMERATE_RESULTS
                                for (ui l = 0; l < cn_count[1]; ++l) {
                                    ui v3 = cn[1][l];

                                    for (ui m = 0; m < cn_count[0]; ++m) {
                                        ui v4 = cn[0][m];
                                        if (v4 != v3)
                                            embedding_count += 1;
                                    }
                                }
#else
                                ui temp_count = cn_count[1];
                                embedding_count += (cn_count[0] - 2) * temp_count;
#endif
                            }
                        }
                    }
                }
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }
}

void FSE::ParallelProcessTripleTriangles(ui *idx, ui *cn_count, ui **cn, ui *offset,
                                         ui *cache, TaskSlot &task, size_t &embedding_count) {

    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {
        ui v0 = i;
        ui v0_neighbors_count;
        const ui* v0_neighbors = data_graph_->Neighbors(v0, v0_neighbors_count);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_neighbors_count;

        if (data_graph_->Degree(v0) >= 4) {
            // Compute the common neighbors.
            for (ui j = 0; j <= edge_begin; ++j) {
                offset[j] = 0;
            }

            for (ui j = edge_begin; j < v0_neighbors_count; ++j) {
                ui v1 = v0_neighbors[j];
                offset[j + 1] = offset[j];

                if (data_graph_->Degree(v1) >= 3) {
                    ui v1_neighbors_count;
                    const ui* v1_neighbors = data_graph_->Neighbors(v1, v1_neighbors_count);

                    ui temp_count;
                    ComputeSetIntersection::ComputeCandidates(v0_neighbors, v0_neighbors_count, v1_neighbors, v1_neighbors_count, cache + offset[j], temp_count);
                    offset[j + 1] = offset[j] + temp_count;
                }
            }

            // Start enumeration.
            for (ui j = edge_begin; j < edge_end; ++j) {
                ui v1 = v0_neighbors[j];

                task.edge_size_ -= 1;
                SplitTask(plus_, v0, j + 1, task, edge_end);

                if (data_graph_->Degree(v1) >= 3) {
                    ui v0_v1_cn_count = offset[j + 1] - offset[j];
                    ui* v0_v1_cn = cache + offset[j];

                    if (v0_v1_cn_count >= 2) {
                        ui start_index = Utility::FirstIndexGreaterThan(v0_v1_cn, 0, v0_v1_cn_count, v1);

                        for (ui k = start_index; k < v0_v1_cn_count; ++k) {
                            ui v2 = v0_v1_cn[k];

                            if (data_graph_->Degree(v2) >= 3) {
                                ui v2_index = Utility::BinarySearch(v0_neighbors, 0, v0_neighbors_count, v2);

                                ui v0_v2_cn_count = offset[v2_index + 1] - offset[v2_index];
                                ui *v0_v2_cn = cache + offset[v2_index];

                                if (v0_v2_cn_count >= 2) {
#ifdef ENUMERATE_RESULTS
                                    for (ui l = 0; l < v0_v1_cn_count; ++l) {
                                        ui v3 = v0_v1_cn[l];
                                        if (v3 != v2) {
                                            for (ui m = 0; m < v0_v2_cn_count; ++m) {
                                                ui v4 = v0_v2_cn[m];
                                                if (v4 != v1 && v4 != v3)
                                                    embedding_count += 1;
                                            }
                                        }
                                    }
#else
                                    unsigned int temp_count;
                                    ComputeSetIntersection::ComputeCandidates(v0_v1_cn, v0_v1_cn_count, v0_v2_cn, v0_v2_cn_count, temp_count);

                                    embedding_count += (v0_v1_cn_count - temp_count - 1) * (v0_v2_cn_count - 1) + temp_count * (v0_v2_cn_count - 2);
#endif
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }

}

void FSE::ParallelProcessQuadrupleTriangles(ui *idx, ui *cn_count, ui **cn,
                                               ui *offset, ui *cache, TaskSlot &task,
                                               size_t &embedding_count) {
    for (ui i = task.edge_offset_begin_; i <= task.edge_offset_end_; ++i) {
        ui v0 = i;
        ui v0_nbrs_cnt;
        const ui* v0_nbrs = data_graph_->Neighbors(v0, v0_nbrs_cnt);

        ui edge_begin = v0 == task.edge_offset_begin_ ? task.edge_begin_ : 0;
        ui edge_end = v0 == task.edge_offset_end_ ? task.edge_end_ : v0_nbrs_cnt;

        if (v0_nbrs_cnt >= 5) {
            // Compute the common neighbors.

            offset[0] = 0;
            for (ui j = 0; j < v0_nbrs_cnt; ++j) {
                ui v1 = v0_nbrs[j];
                offset[j + 1] = offset[j];

                if (data_graph_->Degree(v1) >= 3) {
                    ui v1_nbr_cnt;
                    const ui* v1_nbrs = data_graph_->Neighbors(v1, v1_nbr_cnt);

                    ui temp_count;
                    ComputeSetIntersection::ComputeCandidates(v0_nbrs, v0_nbrs_cnt, v1_nbrs, v1_nbr_cnt, cache + offset[j], temp_count);
                    offset[j + 1] = offset[j] + temp_count;
                }
            }

            for (ui j = edge_begin; j < edge_end; ++j) {
                task.edge_size_ -= 1;
                SplitTask(plus_, v0, j + 1, task, edge_end);
                ui v1 = v0_nbrs[j];

                ui v0_v1_cn_count = offset[j + 1] - offset[j];
                ui* v0_v1_cn = cache + offset[j];

                if (v0_v1_cn_count >= 2) {
                    for (ui k = 0; k < v0_v1_cn_count - 1; ++k) {
                        ui v2 = v0_v1_cn[k];

                        if (data_graph_->Degree(v2) < 3)
                            continue;

                        ui v2_index = Utility::BinarySearch(v0_nbrs, 0, v0_nbrs_cnt, v2);

                        ui v0_v2_cn_count = offset[v2_index + 1] - offset[v2_index];
                        ui* v0_v2_cn = cache + offset[v2_index];

                        if (v0_v2_cn_count >= 2) {
                            for (ui l = k + 1; l < v0_v1_cn_count; ++l) {
                                ui v3 = v0_v1_cn[l];

                                if (data_graph_->Degree(v3) < 3)
                                    continue;

                                ui v3_index = Utility::BinarySearch(v0_nbrs, 0, v0_nbrs_cnt, v3);

                                ui v0_v3_cn_count = offset[v3_index + 1] - offset[v3_index];
                                ui* v0_v3_cn = cache + offset[v3_index];

#ifdef ENUMERATE_RESULTS
                                for (ui m = 0; m < v0_v2_cn_count; ++m) {
                                    ui v4 = v0_v2_cn[m];
                                    if (v4 != v1 && v4 != v3) {
                                        for (ui n = 0; n < v0_v3_cn_count; ++n) {
                                            ui v5 = v0_v3_cn[n];
                                            if (v5 != v1 && v5 != v4 && v5 != v2)
                                                embedding_count += 1;
                                        }
                                    }
                                }

#else
                                if (v0_v3_cn_count >= 2) {
                                    unsigned int v4_temp_count = v0_v2_cn_count - 1;
                                    unsigned int v5_temp_count = v0_v3_cn_count - 1;
                                    if (Utility::Contain(v0_v2_cn, 0, v0_v2_cn_count, v3)) {
                                        v4_temp_count -= 1;
                                        v5_temp_count -= 1;
                                    }

                                    unsigned int temp_count;
                                    ComputeSetIntersection::ComputeCandidates(v0_v2_cn, v0_v2_cn_count, v0_v3_cn, v0_v3_cn_count, temp_count);
                                    temp_count -= 1;

                                    if (v4_temp_count >= temp_count && v5_temp_count > 0) {
                                        unsigned int temp_sum = (v4_temp_count - temp_count) * v5_temp_count + temp_count * (v5_temp_count - 1);

                                        embedding_count += temp_sum;
                                    }
                                }
#endif
                            }
                        }
                    }
                }
            }
        }
        else {
            task.edge_size_ -= (edge_end - edge_begin);
        }
    }
}

void FSE::PrintStatistics() {
    cout << PRINT_SEPARATOR << endl;
    double p_thread_execution_time_sum = 0;
    for (ui i = 0; i < thread_num_; ++i) {
        p_thread_execution_time_sum += p_thread_execution_time_[i];
        cout << "Thread " << i << " Execution Time: " << setprecision(4) << p_thread_execution_time_[i] << " seconds." << endl;
    }

    cout << "Thread Execution Time Sum: " << setprecision(4) << p_thread_execution_time_sum << " seconds." << endl;
    cout << PRINT_SEPARATOR << endl;
    task_queue_->PrintStatistics();
    cout << PRINT_SEPARATOR << endl;
    cout << "Set Intersection Invoke Count: " << set_intersection_invoke_cnt_ << endl;
    cout << PRINT_SEPARATOR << endl;
}