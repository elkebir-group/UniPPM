//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "ILP.h"
#include "Input.h"
#include "AncestryGraph.h"
#include "Hashing.h"

ILP_base::ILP_base(const Input &data, const AncestryGraph &GF, GRBEnv & env):
data(data),GF(GF), model(env),
f(data.m, std::vector<GRBVar>(data.n)),
root(data.n),
arc(GF.arc_set.size()),
arc_f(data.m,std::vector<GRBVar>(GF.arc_set.size()))
{
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    for (int i = 0; i < data.m; i++) {
        for (int j = 0; j < data.n; j++) {
            f[i][j] = model.addVar(data.data[i][j].first, data.data[i][j].second, 0, GRB_CONTINUOUS,
                                   std::string("f") + char(i + 1) + char(j + 1));
        }
    }
    for (int i = 0; i < GF.arc_set.size(); i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              std::string("a") + char(GF.arc_set[i].first + 1) + char(GF.arc_set[i].second + 1));

        for (int j = 0; j < data.m; j++) {
            arc_f[j][i] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                       std::string("fa") + char(GF.arc_set[i].first + 1) +
                                       char(GF.arc_set[i].second + 1));
            model.addConstr(arc_f[j][i] <= arc[i]);
            model.addConstr(arc_f[j][i] <= f[j][GF.arc_set[i].second]);
            model.addConstr(arc_f[j][i] >= arc[i] + f[j][GF.arc_set[i].second] - 1);
        }
    }

    GRBLinExpr sum;
    for (int i = 0; i < data.m; i++) {
        for (int p = 0; p < data.n; p++) {
            for (auto q: GF.possible_children[p]) {
                sum += arc_f[i][GF.arc_set_index[p][q]];
            }
            model.addConstr(f[i][p] >= sum);
            sum.clear();
        }
    }

    ////tree constraints
    //////one root
    for (int i = 0; i < data.n; i++) {
        root[i] = model.addVar(0,1,0,GRB_INTEGER,
                               std::string("r")+char(i+1));
        sum+=root[i];
    }
    model.addConstr(sum==1);
    sum.clear();
    //////one incoming edge
    for (int i = 0; i < data.n; i++) {
        sum+=root[i];
        for(auto j:GF.possible_parent[i]){
            sum+=arc[GF.arc_set_index[j][i]];
        }
        model.addConstr(sum==1);
        sum.clear();
    }
    //////TODO:prevent cycle
    std::vector<GRBVar> cycle(data.n);
    for (int i = 0; i < data.n; i++) {
        cycle[i]=model.addVar(0,data.n-1,0,GRB_CONTINUOUS,
                              std::string("d")+char(i));
//        model.addConstr(cycle[i]<=data.n-root[i]*data.n); //Theoretically not necessary
    }
    for (int i = 0; i < GF.arc_set.size(); i++){
        model.addConstr(cycle[GF.arc_set[i].second] >= cycle[GF.arc_set[i].first]+arc[i]*data.n-data.n+1);
    }

    model.update();
}

ILP::ILP(const ILP_base &base, Hashing &hash, int n_constr):hash(hash), model(base.model), arc(base.arc.size()) {
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
        //TODO: anyother ways or hashing functions?
        tmp = model.addVar(0,(base.data.n>>1),0,GRB_INTEGER, std::string("h")+char(i+1));
        model.addConstr(sum == hash.hash_f_coef[0]+2*tmp);
    }
    model.optimize();
}
