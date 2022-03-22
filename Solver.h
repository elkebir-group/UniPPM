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
    friend class CNF;
public:
    explicit Solver(const AncestryGraph & In, int timeout = 10000, int rec_t=1000);

    ~Solver();

    bool attempt(CNF &, std::vector<std::pair<int,int> > * = NULL);

    void sampling(int n_sample, std::map<std::vector<std::pair<int,int> >,int > & res);


    CNF& self_solver();

    const AncestryGraph &In;
    
    void interpret(const std::vector<int> & vars, std::vector<std::pair<int,int> > & ans);

    void set_up_recursive();

private:

    std::vector<std::vector<uint32_t> > CNF_recursive_sets;

    std::map<int,std::pair<int,int> > var2edge;
    std::map<std::pair<int,int>, CMSat::Lit> edge2var;

    std::vector<CMSat::Lit> to_bin_list(int n, const int &N);

    CNF F;

    int timeout,rec_t;
};


#endif //UNIPPM_SOLVER_H
