//
// Created by ssunah on 12/4/17.
//

#ifndef FSE_TASK_H
#define FSE_TASK_H


class TaskSlot {
public:
    TaskSlot() {}
    TaskSlot(const unsigned int edge_offset_begin, const unsigned int edge_offset_end,
             const unsigned int edge_begin, const unsigned int edge_end, const unsigned int edge_size) {
        edge_offset_begin_ = edge_offset_begin;
        edge_offset_end_ = edge_offset_end;
        edge_begin_ = edge_begin;
        edge_end_ = edge_end;
        edge_size_ = edge_size;
    }

    unsigned int edge_offset_begin_;
    unsigned int edge_offset_end_;
    unsigned int edge_begin_;
    unsigned int edge_end_;
    unsigned int edge_size_;
};


#endif //FSE_TASK_H
