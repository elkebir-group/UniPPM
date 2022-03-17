//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "Solver.h"
#include <map>
#include <chrono>
#include <utility>

Solver::Solver(const AncestryGraph &in, int rec_size, int rec_T, int rec_min):F(),In(in),rec_size(rec_size),
rec_para(rec_T,rec_min){

    std::vector<CMSat::Lit> r_vars;
    if (In.In.r < 0){
        //if no specific root, exact one root
        r_vars = F.new_vars(In.In.n);
        F.Exact_One(r_vars);
    }

    for (int i = 0; i < In.E.size(); i++){
        int var = F.new_var(true);
        var2edge[var] = In.E[i];
        edge2var[In.E[i]] = CMSat::Lit(var,false);
    }

    // exact one parent, except root
    for (int j = 0; j < In.In.n; j++) {
        int _e_o_size;
        if (In.In.r < 0) _e_o_size = In.ind(j).size() + 1;
        else if (j == In.In.r) continue;
        else _e_o_size = In.ind(j).size();

        std::vector<CMSat::Lit> j_p(_e_o_size);
        std::vector<CMSat::Lit>::iterator it = j_p.begin();
        if (In.In.r < 0) {
            *it = r_vars[j];
            it++;
        }
        for (int i = 0; i < In.ind(j).size(); i++) {
            int i_ = In.ind(j)[i];
            *it = edge2var[std::pair<int, int>(i_, j)];

            it++;
        }
        F.Exact_One(j_p);
    }

    //generating Freqs.
    std::vector< std::vector< std::vector<CMSat::Lit> > > F_var
        (In.In.m,std::vector< std::vector<CMSat::Lit> >(In.In.n));
    if (In.is_single){
        for (int p = 0; p < In.In.m; p++)
            for (int i = 0; i < In.In.n; i++)
                F_var[p][i] = to_bin_list(In.In.F_l[p][i],In.In.n_bits);
    }
    else{
        for (int p = 0; p < In.In.m; p++)
            for (int i = 0; i < In.In.n; i++)
                F_var[p][i] = F.new_vars(In.In.n_bits);
    }

    //Sum Condition
    for (int i = 0; i < In.In.n; i++) {
        for (int p = 0; p < In.In.m; p++){
            std::vector<CMSat::Lit> Sum(In.In.n_bits,F.False), tmp(In.In.n_bits);
            for (auto j = In.outd(i).begin(); j != In.outd(i).end(); j++){
                for (int k = 0; k < In.In.n_bits; k++){
                    tmp[k] = F.AND(edge2var[std::pair<int,int>(i,*j)], F_var[p][*j][k]);
                }
                Sum = F.add(Sum,tmp);
            }
            if (In.is_single){
                F.leq(Sum,F_var[p][i]);
            }
            else{
                std::vector<CMSat::Lit> F_u_var = to_bin_list(In.In.F_u[p][i], In.In.n_bits),
                                        F_l_var = to_bin_list(In.In.F_l[p][i], In.In.n_bits);
                F.leq(Sum,F_u_var);
                F.eq(F_var[p][i],F.max(Sum,F_l_var));
            }
        }
    }

    //Cycle prevention
    std::vector<std::vector<CMSat::Lit> > relation(In.In.n, std::vector<CMSat::Lit>(In.In.n));
    for (int i = 0; i < In.In.n; i++){
        for (int j = 0; j < In.In.n; j++) {
            if (i==j) {
                relation[i][i] = F.True;
                continue;
            }
            relation[i][j] = CMSat::Lit(F.new_var(), false);
        }
    }
    for (int i = 0; i < In.In.n; i++){
        for (int j = 0; j < In.In.n; j++){
            if (i==j) continue;
            std::vector<CMSat::Lit> inj(In.ind(j).size());
            for (int ij = 0;ij<In.ind(j).size();ij++){
                inj[ij] = F.AND(relation[i][In.ind(j)[ij]],edge2var[std::pair<int,int>(In.ind(j)[ij],j)]);
            }
            F.OR_Lits(inj,relation[i][j]);
        }
    }
    for (int i = 0; i < In.In.n; i++)
        for (int j = 0; j < i; j++) {
            F.add_clause({~relation[i][j], ~relation[j][i]});
        }

}

Solver::~Solver() {
}

std::vector<CMSat::Lit> Solver::to_bin_list(int n, const int &N) {
    std::vector<CMSat::Lit> a(N);
    for(int i=0;i<N;i++){
        a[i]=(n&1)?F.True:F.False;
        n>>=1;
    }
    return a;
}

CNF& Solver::self_solver() {
    return F;
}

bool Solver::attempt(CNF &Fptr, std::vector<std::pair<int,int> > *result) {
    if (result==NULL) {
        std::vector<CMSat::Lit> res = Fptr.Solve(false);
        return (res[0]==F.True);
    }
    std::vector<CMSat::Lit> res = Fptr.Solve(true);
    if (res.empty())
        return false;
    std::vector<int> var_res(res.size());
    for (int i = 0; i < res.size(); i++){
        var_res[i] = (res[i].sign()?-(int)(res[i].var()):(res[i].var()));
    }
    interpret(var_res,*result);
    return true;
}

void Solver::interpret(const std::vector<int> & vars, std::vector<std::pair<int, int> > &ans) {
    ans.clear();
    for(auto it = vars.begin(); it!=vars.end(); it++){
        if (*it > 0) {
            ans.push_back(var2edge[(*it)-1]);
        }
    }
}

void Solver::sampling(int n_sample, std::map<std::vector<std::pair<int, int> >, int> &res) {
    int n_d1 = In.outd(In.In.r).size();
    std::vector<uint32_t> enumerate(n_d1);

    for (int i = 0; i < n_d1; i++){
        enumerate[i] = edge2var[std::pair<int,int>(In.In.r,In.outd(In.In.r)[i])].var();
    }

    std::list<std::pair<std::vector<int>, int> >  unigen_res;
    std::vector<std::pair<int,int> > tmp;
    std::list<std::vector<int> > data;
//    F.Enum_Sampling(enumerate,n_sample, data, rec_T, rec_size);
//    F.Enum_Sampling({},n_sample, data, rec_T, rec_size);
    set_up_recursive();

    auto appmc = new ApproxMC::AppMC;
    if (In.In.n >=15 && In.In.n/In.In.m>=10){
        F.UniPPM_Sampling(nullptr,n_sample,std::pair(0,rec_size),this,data,rec_para);
    }
    else{
        auto Counting = CNF::Counting(F,{},appmc);
        F.UniPPM_Sampling(&Counting,n_sample,std::pair(0,rec_size),this,data,rec_para);
    }




    std::cout<<"[UniPPM] Sampling finished, sorting solutions. "<<std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    data.sort();

    for (auto it = data.begin(), it_pre = data.begin();it!=data.end(); it++) {
        if (it==data.begin() ) {
            unigen_res.emplace_back(*it,1);
            continue;
        }
        if ((*it)!=(*it_pre)){
            unigen_res.emplace_back(*it,1);
        }
        else {
            unigen_res.back().second++;
        }
        it_pre++;
    }

    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "[UniPPM] Sorting takes "
              << (std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())/1e6
              << " seconds.\n[UniPPM] Verifying solutions."<< std::endl;

    for(auto it=unigen_res.begin();it!=unigen_res.end();it++){
        interpret(it->first, tmp);
        res[tmp] = it->second;
    }

    start = std::chrono::high_resolution_clock::now();
    std::cout << "[UniPPM] Verifying takes "
              << (std::chrono::duration_cast<std::chrono::microseconds>(start - stop).count())/1e6
              << " seconds."<< std::endl;
}

void Solver::set_up_recursive() {
    std::vector<int> step_order(In.In.n);
    CNF_recursive_sets.resize(In.In.n);
    for (int i = 0; i < In.In.n; i++){
        step_order[i] = i;
    }

    struct comp{
        bool operator()(const int & a, const int &b) const{
            return ptr->In.ind(a).size() > ptr->In.ind(b).size();
        }
        Solver * ptr;
        comp(Solver * ptr):ptr(ptr){}
    };

    std::sort(step_order.begin(),step_order.end(),comp(this));

    for (int i=0;i<step_order.size();i++){
        if (In.ind(step_order[i]).empty()) {
            CNF_recursive_sets.resize(i);
            break;
        }
        for(auto j: In.ind(step_order[i])){//j,i
            CNF_recursive_sets[i].push_back(edge2var[std::pair(j,step_order[i])].var());
        }
    }
}










