//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "Likelihood.h"
#include <cmath>

#define logbinom(N,K,P,NN) (lgamma((N)+1)-lgamma((K)+1)-lgamma((N)-(K)+1)+(K)*log_(P,NN)+((N)-(K))*log_(1-(P),NN))
#define logcomb(N,K) (lgamma((N)+1)-lgamma((K)+1)-lgamma((N)-(K)+1))

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

}

double Likelihood::LLH(const std::vector<std::pair<int,int> > & edge_set){

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

    std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("SCIP"));
    std::vector<std::vector<operations_research::MPVariable*> >
            f(In.m,std::vector<operations_research::MPVariable *>(In.n,NULL));
    std::vector<std::vector<std::vector<operations_research::MPVariable*> > >
            lambda(In.m,std::vector<std::vector<operations_research::MPVariable *> >
            (In.n,std::vector<operations_research::MPVariable *> (n_split,NULL) ));
    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < In.n; j++) {
            f[i][j] = solver ->MakeNumVar(F_lower[i][j],F_upper[i][j],
                                          "f["+std::to_string(i)+"]["+std::to_string(j)+"]");
            for(int k = 0; k < n_split; k++){
                lambda[i][j][k] = solver ->MakeNumVar(0,1,
                                                      "lambda["+std::to_string(i)+"]["+
                                                      std::to_string(j)+"]["+
                                                      std::to_string(k)+"]");
            }
        }
    }

    for (int i = 0; i < In.m; i++){
        for (int j = 0 ; j < In.n; j++){
            auto sc = solver ->MakeRowConstraint(0,solver->infinity());
            sc ->SetCoefficient(f[i][j],1);
            for (auto it = outd[j].begin();it!=outd[j].end();it++)
                sc ->SetCoefficient(f[i][*it],-1);
        }
    }

    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < In.n; j++) {
            auto fc = solver -> MakeRowConstraint(0,0);
            fc->SetCoefficient(f[i][j],1);
            auto lc = solver ->MakeRowConstraint(1,1);
            for (int k = 0; k < n_split; k++){
                fc ->SetCoefficient(lambda[i][j][k],-split[i][j][k]);
                lc ->SetCoefficient(lambda[i][j][k],1);
            }
        }
    }

    auto objective = solver ->MutableObjective();
    for (int i = 0; i < In.m; i++){
        for (int j = 0; j < In.n; j++){
            for (int k = 0; k < n_split; k++){
                if(mul && j==(In.n-1)) continue;
                objective ->SetCoefficient(lambda[i][j][k],
                               Reads.var[i][j]*log_(split[i][j][k],In.n_bits)+
                               Reads.ref[i][j]*log_(1-split[i][j][k],In.n_bits));
            }
        }
    }
    objective->SetMaximization();
    solver->Solve();
    auto tmp_ans = objective->Value();
    double ans = 0;
    for (int i = 0; i < In.m; i++) {
        for (int j = 0; j < In.n; j++) {
            if(mul && j==(In.n-1)) continue;
            ans += logbinom(Reads.var[i][j]+Reads.ref[i][j],Reads.var[i][j],f[i][j]->solution_value(),In.n_bits);
        }
    }

    return ans;
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
    assert(reads.m==intervals.m && reads.n==(intervals.n-mul));
    double ans = 0;
    for(int i = 0; i < reads.m; i++)
        for(int j = 0; j < reads.n; j++) {
            double a = logbinom(reads.var[i][j] + reads.ref[i][j], reads.var[i][j], intervals.F_l[i][j], n_bits),
                    b = logbinom(reads.var[i][j] + reads.ref[i][j], reads.var[i][j], intervals.F_l[i][j], n_bits);
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

