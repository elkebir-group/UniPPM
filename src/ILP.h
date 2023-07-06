//
// Created by Yuanyuan Qi on 7/6/23.
//

#ifndef UNIPPM_ILP_H
#define UNIPPM_ILP_H


//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_ILP_RANGE_H
#define UNIPPM_ILP_RANGE_H

#include <list>
#include <stack>
#include "gurobi_c++.h"

class Input;
class AncestryGraph;
struct ILP_base{
    const Input &data;
    const AncestryGraph & GF;
    GRBModel model;
    std::vector<GRBVar> arc;
    ILP_base(const AncestryGraph &GF, GRBEnv & env);
};

class Hashing;
struct ILP {
    GRBModel model;
    Hashing &hash;
    std::vector<GRBVar> arc;
    ILP(const ILP_base& base, Hashing & hash, int n_constr);
};

#endif //UNIPPM_ILP_RANGE_H


#endif //UNIPPM_ILP_H
