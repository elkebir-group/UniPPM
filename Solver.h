//
// Created by Yuanyuan Qi on 2/14/22.
//

#ifndef UNIPPM_SOLVER_H
#define UNIPPM_SOLVER_H

#include "AncestryGraph.h"
#include "CNF.h"
#include <approxmc/approxmc.h>
#include <unigen/unigen.h>
#include <map>
#include <mutex>

class Solver {
    friend class Parallel_Solver;
public:
    Solver(const AncestryGraph & In, int n_threads);

    ~Solver();

    bool attempt(CNF &, std::vector<std::pair<int,int> > * = NULL);

    void sampling(int n_sample, std::map<std::vector<std::pair<int,int> >,int > & res);


    CNF& self_solver();

    const AncestryGraph &In;
    
    void interpret(const std::vector<int> & vars, std::vector<std::pair<int,int> > & ans);

private:

    std::map<int,std::pair<int,int> > var2edge;
    std::map<std::pair<int,int>, CMSat::Lit> edge2var;

    std::vector<CMSat::Lit> to_bin_list(int n, const int &N);

    CNF F;

};


#endif //UNIPPM_SOLVER_H
