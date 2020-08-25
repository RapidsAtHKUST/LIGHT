//
// Created by ssunah on 12/4/17.
//

#ifndef FSE_TASKPOOL_H
#define FSE_TASKPOOL_H


#include <pthread.h>
#include "TaskSlot.h"

class TaskPool {
public:
    TaskPool(const unsigned int queue_size);
    ~TaskPool();
    bool Empty() { return queue_count_ == 0; }

    bool Push(const TaskSlot& task);
    bool Pop(TaskSlot& task);

    void PrintStatistics();
private:
    // Concurrent queue data member
    TaskSlot* concurrent_queue_;
    unsigned int queue_front_;
    unsigned int queue_count_;
    unsigned int queue_capacity_;
    pthread_mutex_t queue_mutex_;

    // Queue statistics
    size_t task_count_;
};


#endif //FSE_TASKPOOL_H
