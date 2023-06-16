//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_ILP_H
#define UNIPPM_ILP_H

#include "gurobi_c++.h"

class Input;
class AncestryGraph;
struct ILP_base{
    const Input &data;
    const AncestryGraph & GF;
    GRBModel model;
    std::vector<std::vector<GRBVar> > f;
    std::vector<GRBVar> root;
    std::vector<GRBVar> arc;
    std::vector<std::vector<GRBVar> > arc_f;
    ILP_base(const Input &data, const AncestryGraph &GF, GRBEnv & env);
};

class Hashing;
struct ILP {
    GRBModel model;
    Hashing &hash;
    std::vector<GRBVar> arc;
    ILP(const ILP_base& base, Hashing & hash, int n_constr);
};


#endif //UNIPPM_ILP_H
