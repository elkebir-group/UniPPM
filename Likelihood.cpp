//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "Likelihood.h"
#include <cmath>
#include "gurobi_c++.h"
#include <cassert>

#define logcomb(N,K) (lgamma((N)+1)-lgamma((K)+1)-lgamma((N)-(K)+1))
#define logbinom(N,K,P,NN) (logcomb(N,K)+(K)*log_(P,NN)+((N)-(K))*log_(1-(P),NN))

Likelihood::Likelihood(const Input_int &Gf, const Input_Reads &reads, const int &N_intervals, bool mul):
In(Gf),Reads(reads),n_split(N_intervals+1),mul(mul)
{
    F_lower = new double *[In.m];
    F_upper = new double *[In.m];
    ptr1 = new double[In.m*In.n*2];
    split = new double ** [In.m];
    ptpt2 = new double * [In.m*In.n];
    ptr2 = new double[In.m*In.n*n_split];

    for (int i = 0; i < In.m; i++) {
        F_lower[i] = ptr1 + i * In.n;
        F_upper[i] = ptr1 + (i + In.m) * In.n;
        split[i] = ptpt2 + i * In.n;
        for (int j = 0; j < In.n; j++) {
            F_lower[i][j] = (double)In.F_l[i][j] / (1<<In.n_bits);
            F_upper[i][j] = (double)In.F_u[i][j] / (1<<In.n_bits);
            split[i][j] = ptr2 + (i*In.n+j) * n_split;
            for (int k = 0; k < n_split; k++){
                split[i][j][k] = F_lower[i][j]+(F_upper[i][j]-F_lower[i][j])/(n_split-1)*k;
            }
        }
    }

//    cl_ptr = new int[reads.n];
//    cluster_info = new int*[reads.n_cluster];
//    for(int i = 0, *tmp=cl_ptr ; i < reads.n_cluster ;i++){
//        cluster_info[i] = tmp;
//        tmp += reads.cl_cnt[i];
//    }

}

double Likelihood::LLH(const std::vector<std::pair<int,int> > & edge_set, std::vector<std::vector<double> > *f_r){

    std::vector<std::vector<int> > ind(In.n), outd(In.n);
    for (auto it = edge_set.begin(); it!=edge_set.end();it++){
        ind[it->second].push_back(it->first);
        outd[it->first].push_back(it->second);
    }
    int r;
    for (int i = 0; i < In.n; i++){
        if (ind[i].empty()){
            r = i;
            break;
        }
    }

//    auto solver(operations_research::MPSolver::CreateSolver("SCIP"));
    GRBEnv env;
    GRBModel model(env);
    std::vector<std::vector<GRBVar> >
            f(In.m,std::vector<GRBVar>(In.n));
    std::vector<std::vector<std::vector<GRBVar> > >
            lambda(In.m,std::vector<std::vector<GRBVar> >
            (In.n,std::vector<GRBVar> (n_split) ));
    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < In.n; j++) {
            f[i][j] = model.addVar(F_lower[i][j], F_upper[i][j], 0, GRB_CONTINUOUS,
                                          "f["+std::to_string(i)+"]["+std::to_string(j)+"]");
            for(int k = 0; k < n_split; k++){
                lambda[i][j][k] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                                      "lambda["+std::to_string(i)+"]["+
                                                      std::to_string(j)+"]["+
                                                      std::to_string(k)+"]");
            }
        }
    }

    GRBLinExpr sum;
    for (int i = 0; i < In.m; i++){
        for (int j = 0 ; j < In.n; j++){
            for (auto it = outd[j].begin();it!=outd[j].end();it++)
                sum += f[i][*it];
            model.addConstr(f[i][j] >= sum);
            sum.clear();
        }
    }

    GRBLinExpr sum2;
    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < In.n; j++) {
            for (int k = 0; k < n_split; k++){
                sum += lambda[i][j][k];
                sum2 += lambda[i][j][k] * split[i][j][k];
            }
            model.addConstr(sum == 1);
            model.addConstr(sum2 == f[i][j]);
            sum.clear();
            sum2.clear();
        }
    }

    for (int i = 0; i < In.m; i++){
        for (int j = 0; j < Reads.n; j++){
            for (int k = 0; k < n_split; k++){
                // for the case where we allow multiple children from the normal clone (mul == true),
                // there will be an artificial mutation (indexed by n-1), which must be skipped as it
                // won't have read counts associated with it
//                if(mul && j==(In.n-1)) continue;

                sum += lambda[i][Reads.cl[j]][k] * (Reads.var[i][j]*log_(split[i][Reads.cl[j]][k],In.n_bits)+
                                          Reads.ref[i][j]*log_(1-split[i][Reads.cl[j]][k],In.n_bits));
            }
        }
    }

    model.setObjective(sum, GRB_MAXIMIZE);
//    model.update();
    model.optimize();

    auto obj_val_inferred = model.get(GRB_DoubleAttr_ObjVal);
    double log_binom_coeffs = 0;
    double obj_val = 0;
    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < Reads.n; j++) {
//            if(mul && j==(In.n-1)) continue;
            log_binom_coeffs += logcomb(Reads.var[i][j] + Reads.ref[i][j],
                                        Reads.var[i][j]);
            obj_val += logbinom(Reads.var[i][j] + Reads.ref[i][j],
                            Reads.var[i][j],
                            f[i][Reads.cl[j]].get(GRB_DoubleAttr_X), In.n_bits);
        }
    }

    std::cout << "GUROBI: " << obj_val_inferred + log_binom_coeffs << " -- recomputed: " << obj_val << std::endl;
    if (f_r!= nullptr){
        for(int i = 0; i < In.m; i++){
            for (int j = 0; j < In.n; j++){
                (*f_r)[i][j] = f[i][j].get(GRB_DoubleAttr_X);
            }
        }
    }

    return obj_val;
}


Likelihood::~Likelihood() {
    delete [] F_upper;
    delete [] F_lower;
    delete [] ptr1;
    delete [] split;
    delete [] ptpt2;
    delete [] ptr2;
}

double lower_bound_llh(const Input_Reads & reads, const Input & intervals, const int & n_bits, bool mul) {
    assert(reads.m==intervals.m && reads.n_cluster==(intervals.n-mul));
    double ans = 0;
    for(int i = 0; i < reads.m; i++)
        for(int j = 0; j < reads.n; j++) {
            double a = logbinom(reads.var[i][j] + reads.ref[i][j], reads.var[i][j], intervals.F_l[i][reads.cl[j]], n_bits),
                    b = logbinom(reads.var[i][j] + reads.ref[i][j], reads.var[i][j], intervals.F_l[i][reads.cl[j]], n_bits);
            ans += std::max(a,b);
        }
    return ans;
}

double log_(double x, int n_bits) {
    if(x<(1.0/(1<<n_bits))){
        return -log(n_bits)-(1.0/(1<<n_bits)-x)*(1<<n_bits)*(1<<n_bits);
    }
    return log(x);
}

#undef logbinom

