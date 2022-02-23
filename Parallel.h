//
// Created by Yuanyuan Qi on 2/23/22.
//

#ifndef UNIPPM_PARALLEL_H
#define UNIPPM_PARALLEL_H

#include "Solver.h"


class Parallel_Solver {
    static void sampling_parallel(Parallel_Solver &);
    int n_threads,n_sample;

    int job_index, phase, sum_sols, c_sols;
    std::mutex job_index_mutex, c_mutex, res_lock;
    std::vector<CNF> jobs;
    std::vector<int> n_sols;

    Solver & sol;

    std::map<std::vector<std::pair<int, int> >, int> & res;

public:
    Parallel_Solver(Solver &, int N_threads, int N_sample, std::map<std::vector<std::pair<int, int> >, int> &);

    void sample();
};


#endif //UNIPPM_PARALLEL_H
