//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "ILP.h"
#include "Input.h"
#include "AncestryGraph.h"
#include "Hashing.h"

ILP_base::ILP_base(const AncestryGraph &GF, GRBEnv & env):
data(GF.data),GF(GF), model(env), f(data.m, std::vector<GRBVar>(data.n)),
//root(data.n),
arc(GF.arc_set.size()), arc_f(data.m,std::vector<GRBVar>(GF.arc_set.size()))
{
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    model.set(GRB_IntParam_LazyConstraints, 1);
    for (int i = 0; i < data.m; i++) {
        for (int j = 0; j < data.n; j++) {
            f[i][j] = model.addVar(data.data[i][j].first, data.data[i][j].second, 0, GRB_CONTINUOUS,
                                   "f[" + std::to_string(i ) + "," + std::to_string(j)+"]");
        }
    }
    for (int i = 0; i < GF.arc_set.size(); i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(GF.arc_set[i].first) + "," +
                              std::to_string(GF.arc_set[i].second)+"]");

        for (int j = 0; j < data.m; j++) {
            arc_f[j][i] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                       "f_a[" + std::to_string(j) + ";" + std::to_string(GF.arc_set[i].first) +
                                       "," + std::to_string(GF.arc_set[i].second)+"]");
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
            model.addConstr(f[i][p] >= sum );
            sum.clear();
        }
    }

//    ////tree constraints
//    //////one root
//    for (int i = 0; i < data.n; i++) {
//        root[i] = model.addVar(0, 1, 0, GRB_INTEGER,
//                               "r[" + std::to_string(i)+"]");
//        sum += root[i];
//    }
//    model.addConstr(sum == 1);
//    sum.clear();
    //////one incoming edge
    for (int i = 0; i < data.n; i++) {
//        sum += root[i];
        if(i==data.r) continue;
        for (auto j: GF.possible_parent[i]) {
            sum += arc[GF.arc_set_index[j][i]];
        }
        model.addConstr(sum == 1);
        sum.clear();
    }
//    //////prevent cycle
//    std::vector<GRBVar> cycle(data.n);
//    for (int i = 0; i < data.n; i++) {
//        cycle[i] = model.addVar(0, data.n - 1, 0, GRB_CONTINUOUS,
//                                "d_help[" + std::to_string(i)+"]");
//        model.addConstr(cycle[i]<=data.n-root[i]*data.n); //Theoretically not necessary
//    }
//    for (int i = 0; i < GF.arc_set.size(); i++) {
//        model.addConstr(cycle[GF.arc_set[i].second] >= cycle[GF.arc_set[i].first] + arc[i] * data.n - data.n + 1);
//    }
//
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
//        tmp = model.addVar(0,(base.data.n>>1),0,GRB_INTEGER, "hash_"+std::to_string(i));
//        model.addConstr(sum == hash.hash_f_coef[0]+2*tmp);

        model.addConstr(sum == ((base.data.n-1)>>1));
//        if (hash.hash_f_coef[0]){
//            model.addConstr(sum <= ((base.data.n-1)>>1) );
//        }
//        else {
//            model.addConstr(sum >= ((base.data.n-1)>>1)+1 );
//        }
        sum.clear();
    }
    auto pt = new ILP_Callback(*this, base.GF);
    model.setCallback(pt);
    model.optimize();
    delete pt;
}

ILP_Callback::ILP_Callback(const ILP &LP, const AncestryGraph &GF):LP(LP), GF(GF), data(GF.data){
}

void ILP_Callback::callback() {
    if (where != GRB_CB_MIPSOL) return;

    std::vector<std::list<int> > ee (data.n);
    std::vector<bool>visited(data.n, false);
    std::list<std::pair<int,std::list<int>::iterator > > dfs;

    for (int i = 0; i < LP.arc.size(); i++) {
        if (getSolution(LP.arc[i]) >= 0.99) {
            ee[GF.arc_set[i].first].push_back(GF.arc_set[i].second);
        }
    }
    bool flag = false;
    int loop_at;
    for (int i = -1; i < data.n; i++) {
        if (i < 0) {
            visited[data.r] = true;
            dfs.push_front({data.r, ee[data.r].begin()});
        } else {
            if (visited[i]) continue;
            visited[i] = true;
            dfs.push_front({i, ee[i].begin()});
        }
        while (!dfs.empty()) {
            if (dfs.front().second == ee[dfs.front().first].end()) {
                dfs.pop_front();
                continue;
            }
            auto next_node = *dfs.front().second;
            dfs.front().second++;
            if (visited[next_node]) {
                flag = true;
                loop_at = next_node;
                break;
            }
            visited[next_node] = true;
            dfs.push_front({next_node, ee[next_node].begin()});

        }
        if (flag) break;
    }

//    printf("%d %d\n",loop_at, flag);
    if (flag) {
        GRBLinExpr sum;
        int loop_size = 0, fr = -1, to = loop_at;
        while (fr != loop_at) {
            fr = dfs.front().first;
            sum += LP.arc[GF.arc_set_index[fr][to]];
            printf("%d->%d;", fr, to);
            to = fr;
            dfs.pop_front();
            loop_size++;
        }
        addLazy(sum <= loop_size - 1);
        printf("\n");
    }

}

