//
// Created by Yuanyuan Qi on 7/6/23.
//

#include "ILP.h"
#include "Input.h"
#include "AncestryGraph.h"
#include "Hashing.h"

ILP_base::ILP_base(const AncestryGraph &GF, GRBEnv & env):
        data(GF.data),GF(GF), model(env),
        arc(GF.arc_set.size())
{
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    model.set(GRB_IntParam_LazyConstraints, 1);
    for (int i = 0; i < GF.arc_set.size(); i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(GF.arc_set[i].first) + "," +
                              std::to_string(GF.arc_set[i].second)+"]");
    }

    GRBLinExpr sum;
    for (int i = 0; i < data.m; i++) {
        for (int p = 0; p < data.n; p++) {
            for (auto q: GF.possible_children[p]) {
                sum += arc[i]*data.data[i][q];
            }
            model.addConstr(data.data[i][p] >= sum );
            sum.clear();
        }
    }

    //////one incoming edge
    for (int i = 0; i < data.n; i++) {
        if(i==data.r) continue;
        for (auto j: GF.possible_parent[i]) {
            sum += arc[GF.arc_set_index[j][i]];
        }
        model.addConstr(sum == 1);
        sum.clear();
    }

    model.update();
}

ILP::ILP(const ILP_base &base, Hashing &hash, int n_constr): hash(hash), model(base.model), arc(base.arc.size()) {
    GRBLinExpr sum;
    GRBVar tmp;
    for (int i = 0; i < arc.size(); i++) {
        arc[i] = model.getVarByName(base.arc[i].get(GRB_StringAttr_VarName));
    }
    for (int i = 0; i < n_constr; i++) {
        hash.generate();
        for (int j = 1; j < hash.hash_f_coef.size(); j++) {
            if (hash.hash_f_coef[j]) {
                sum += arc[j-1];
            }
        }
//        //TODO: anyother ways or hashing functions?
//        tmp = model.addVar(0,(base.data.n>>1),0,GRB_INTEGER, "hash_"+std::to_string(i));
//        model.addConstr(sum == hash.hash_f_coef[0]+2*tmp);
        model.addConstr(sum == ((base.data.n-1)>>1));

        sum.clear();
    }
    model.optimize();
}
