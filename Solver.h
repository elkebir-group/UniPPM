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

class Solver {
public:
    Solver(const AncestryGraph & In, int n_threads);

    ~Solver();

    bool attempt(CNF &, std::vector<int> * = NULL);

    int counting(CNF &);

    std::vector<std::vector<int> > sampling(CNF &,int n_sample);

    CNF& self_solver();

    const AncestryGraph &In;

private:

    std::map<int,std::pair<int,int> > var2edge;
    std::map<std::pair<int,int>, CMSat::Lit> edge2var;

    std::vector<CMSat::Lit> to_bin_list(int n, const int &N);

    CNF F;

    bool divide_n_conquer;
};


#endif //UNIPPM_SOLVER_H
