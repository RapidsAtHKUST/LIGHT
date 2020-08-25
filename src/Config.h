//
// Created by ssunah on 11/30/17.
//

#ifndef FSE_CONFIG_H
#define FSE_CONFIG_H

#define PRINT_SEPARATOR "------------------------------"

typedef unsigned int ui;

#define MINIMUM_TASK_LENGTH 4
#define MAXIMUM_SIZE 8
// Simple nested for loop to compute Cartesian products without any optimizations.
// #define ENUMERATE_RESULTS
// #define COLLECT_RESULTS

/*
 *  0: AVX2; 1: AVX512; 2: Basic;
 */
#define SI 0
/*
 * 0: Hybrid; 1: Merge;
 */
#define HYBRID 0

#endif //FSE_CONFIG_H
