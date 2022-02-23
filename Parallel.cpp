//
// Created by Yuanyuan Qi on 2/23/22.
//

#include "Parallel.h"
#include <thread>

Parallel_Solver::Parallel_Solver(Solver & solver, int N_threads, int N_sample
                                 ,std::map<std::vector<std::pair<int, int> >, int> &Res):
                                 n_threads(N_threads), sol(solver), n_sample(N_sample), res(Res){
    int n_d1 = sol.In.outd(sol.In.In.r).size();
    int n_cases = 1<<n_d1;
    jobs.resize(n_cases,sol.F);
    std::vector<CMSat::Lit> enumerate(n_d1);
    for (int i = 0; i < n_d1; i++){
        enumerate[i] = sol.edge2var[std::pair<int,int>(sol.In.In.r,sol.In.outd(sol.In.In.r)[i])];
    }
    for (int i = 0; i < n_d1; i++){
        jobs[0].add_clause({enumerate[i]});
    }
    for (int ca = 1; ca < n_cases; ca++){
        int cas = ca^(ca-1);
        for(int i = 0; i < n_d1; i++){
            if (cas & 1){
                enumerate[i] = ~enumerate[i];
            }
            cas >>= 1;
            jobs[ca].add_clause({enumerate[i]});
        }
    }
}

void Parallel_Solver::sampling_parallel(Parallel_Solver & p_s ) {
    int job_id;
    while (true) {
        p_s.job_index_mutex.lock();
        job_id = p_s.job_index++;
        p_s.job_index_mutex.unlock();
        if (job_id >= 2 * p_s.jobs.size()) {
            break;
        }
        if (job_id < p_s.jobs.size()) {
            p_s.n_sols[job_id] = p_s.jobs[job_id].Counting();
            p_s.sum_sols += p_s.n_sols[job_id];
            p_s.c_mutex.lock();
            p_s.c_sols++;
            if (p_s.c_sols >= p_s.jobs.size()) {
                p_s.phase = 1;
            }
            p_s.c_mutex.unlock();
            continue;
        }
        while (p_s.phase < 1);
        job_id -= p_s.jobs.size();
        p_s.sol.sampling(p_s.jobs[job_id], (int) (1.0 * p_s.n_sols[job_id] / p_s.sum_sols * p_s.n_sample),
                         p_s.res, &p_s.res_lock);
    }
}

void Parallel_Solver::sample() {
    std::vector<std::thread> pool;
    for(int i = 0; i< n_threads; i++){
        std::thread a(Parallel_Solver::sampling_parallel,std::ref(*this));
    }
}
