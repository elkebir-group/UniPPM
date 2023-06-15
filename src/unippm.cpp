//
// Created by Yuanyuan Qi on 6/14/23.
//
#include <cmath>
#include "input.h"
#include "ILP.h"
#include "Hashing.h"
#include "AncestryGraph.h"

int hash_size;
int n_rep_has=10;
int n_samples;

int main(int argc, char ** argv) {
    Input input(argv[1]);
    AncestryGraph GF(input);
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    ILP_base ilp_base(input, GF, env);
    Hashing hash(GF.arc_set.size() + 1);
    int ub = int(log2(input.n) * (input.n - 2)) + 1;
    int lb = 1;
    int mid = (ub + 10 * lb) / 11;
    printf("%d %d %d\n",lb,ub,mid);
    int n_fail;
    while (ub - lb > 0) {
        n_fail = 0;
        for (int i = 0; i < n_rep_has; i++) {
            ILP test(ilp_base, hash, mid);
            if (test.model.get(GRB_IntAttr_Status)==GRB_INFEASIBLE) {
                n_fail += 1;
            }
        }

        printf("%d %d\n",mid,n_fail);
        if (n_fail > 5) { // Fails a lot, maybe too many constraints
            ub = mid;
        } else { // Succeeds a lot, maybe too few constraints
            lb = mid + 1;
        }
        mid = (lb + ub) / 2;
    }
//    return 0;
//    sampling
    if (argc >= 3) n_samples = std::stoi(argv[2]);
    else n_samples = 100;
    int n_sampled = 0;
    while (n_sampled < n_samples) {
        ILP test(ilp_base, hash, lb);
        if (test.model.get(GRB_IntAttr_Status)==GRB_OPTIMAL) {
            printf("------------------[sample %d]------------------\n", ++n_sampled);
            for (int i = 0; i < test.arc.size(); i++) {
                if (test.arc[i].get(GRB_DoubleAttr_X) > 0.5) {
                    printf("%d %d\n", GF.arc_set[i].first, GF.arc_set[i].second);
                }
            }
        }
    }
    return 0;
}