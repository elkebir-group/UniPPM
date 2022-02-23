//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "AncestryGraph.h"

#include <vector>

AncestryGraph::AncestryGraph(const Input_int &in): In(in), out_d(in.n), in_d(in.n), E(edge_set),
Is_single(true),is_single(Is_single){
    for (int i = 0; i < in.n; i++) {
        for (int j = 0; j < in.n; j++){
            if (i == j || j == In.r) continue;
            bool ij_flag = 1;
            for (int p = 0; p < In.m; p++) {
                if (In.F_l[p][j] > In.F_u[p][i]) {
                    ij_flag = 0;
                    break;
                }
            }

            if (ij_flag) {
                out_d[i].push_back(j);
                in_d[j].push_back(i);
                edge_set.emplace_back(i,j);
            }
        }
    }

    for (int p = 0; p < In.m; p++) {
        for (int i = 0; i < In.n; i++) {
            if(In.F_u[p][i] - In.F_l[p][i] > 0){
                Is_single = false;
                break;
            }
        }
        if(!Is_single) break;
    }
}

const std::vector<int> &AncestryGraph::ind(int i) const {
    return in_d[i];
}

const std::vector<int> &AncestryGraph::outd(int i) const {
    return out_d[i];
}


