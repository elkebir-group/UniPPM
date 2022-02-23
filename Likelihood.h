//
// Created by Yuanyuan Qi on 2/14/22.
//

#ifndef UNIPPM_LIKELIHOOD_H
#define UNIPPM_LIKELIHOOD_H

#include "AncestryGraph.h"
#include "Input.h"
#include "ortools/linear_solver/linear_solver.h"


class Likelihood {
public:
    Likelihood(const Input_int &,const Input_Reads &,const int &);

    ~Likelihood();

    const Input_int &In;
    const Input_Reads &Reads;
    const int n_split;

    double LLH(const std::vector<std::pair<int,int> > & edge_set);

private:

//    std::unique_ptr<operations_research::MPSolver> solver;
//    std::vector<std::vector<operations_research::MPVariable*> > f;
//    std::vector<std::vector<std::vector<operations_research::MPVariable*> > > lambda;

    double **F_lower, **F_upper, *ptr1, ***split, **ptpt2, *ptr2;
};

double log_(double x, int n_bits);

double lower_bound_llh(const Input_Reads &, const Input &, const int &n_bits);

#endif //UNIPPM_LIKELIHOOD_H