//
// Created by ssunah on 12/4/17.
//

#include "TaskPool.h"
#include <iostream>
using namespace std;

TaskPool::TaskPool(const unsigned int queue_size) {
    task_count_ = 0;
    queue_front_ = 0;
    queue_count_ = 0;
    queue_capacity_ = queue_size;

    concurrent_queue_ = new TaskSlot[queue_capacity_];
    pthread_mutex_init(&queue_mutex_, NULL);
}

TaskPool::~TaskPool() {
    delete[] concurrent_queue_;

    pthread_mutex_destroy(&queue_mutex_);
}

bool TaskPool::Push(const TaskSlot &task) {
    bool res = false;

    pthread_mutex_lock(&queue_mutex_);

    if (queue_count_ < queue_capacity_) {
        concurrent_queue_[(queue_front_ + queue_count_) % queue_capacity_] = task;

        queue_count_ += 1;
        res = true;

        task_count_ += 1;
    }

    pthread_mutex_unlock(&queue_mutex_);

    return res;
}

bool TaskPool::Pop(TaskSlot &task) {
    bool res = false;

    pthread_mutex_lock(&queue_mutex_);

    if (queue_count_ > 0) {
        task = concurrent_queue_[queue_front_];
        queue_front_ = (queue_front_ + 1) % queue_capacity_;

        queue_count_ -= 1;
        res = true;
    }

    pthread_mutex_unlock(&queue_mutex_);

    return res;
}

void TaskPool::PrintStatistics() {
    cout << "Total Task Number: " << task_count_ << endl;
}
