//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "CNF.h"
#include <random>
#include <cassert>
#include <chrono>
#include "Callback.h"

CNF::CNF():True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), //appmc(NULL), unigen(NULL),
n_variables(1), //Results(callbackdata),
clauses({std::vector<CMSat::Lit>{ CMSat::Lit(0, false)} }),
ind_vs()
{
}

CNF::CNF(const CNF & cp):True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), //appmc(NULL), unigen(NULL),
clauses(cp.clauses), n_variables(cp.n_variables), ind_vs(cp.ind_vs)//, Results(callbackdata)
{
}

int CNF::new_var(bool is_ind) {
    if (is_ind){
        ind_vs.push_back(n_variables);
    }
    return n_variables++;
}

std::vector<CMSat::Lit> CNF::new_vars(int n,bool is_ind) {
    std::vector<CMSat::Lit> vars(n);
    for (int i=0; i<n; i++){
        vars[i] = CMSat::Lit(new_var(is_ind), false);
    }
    return vars;
}

CMSat::Lit CNF::AND(const CMSat::Lit &a, const CMSat::Lit &b) {
    CMSat::Lit r(new_var(), false);
    AND(a,b,r);
    return r;
}

void CNF::AND(const CMSat::Lit &a, const CMSat::Lit &b, const CMSat::Lit &r) {

    clauses.push_back(std::vector<CMSat::Lit>{r, ~a, ~b});
    clauses.push_back(std::vector<CMSat::Lit>{~r, a});
    clauses.push_back(std::vector<CMSat::Lit>{~r, b});

}

CMSat::Lit CNF::OR(const CMSat::Lit &a, const CMSat::Lit &b) {
    CMSat::Lit r(new_var(), false);
    OR(a,b,r);
    return r;
}

void CNF::OR(const CMSat::Lit &a, const CMSat::Lit &b, const CMSat::Lit &r) {

    clauses.push_back(std::vector<CMSat::Lit>{~r, a, b});
    clauses.push_back(std::vector<CMSat::Lit>{r, ~a});
    clauses.push_back(std::vector<CMSat::Lit>{r, ~b});

}

CMSat::Lit CNF::XOR(const CMSat::Lit &a, const CMSat::Lit &b) {
    CMSat::Lit r(new_var(), false);
    XOR(a,b,r);
    return r;
}

void CNF::XOR(const CMSat::Lit &a, const CMSat::Lit &b, const CMSat::Lit &r) {

    clauses.push_back(std::vector<CMSat::Lit>{~r, a, b});
    clauses.push_back(std::vector<CMSat::Lit>{~r, ~a, ~b});
    clauses.push_back(std::vector<CMSat::Lit>{r, a, ~b});
    clauses.push_back(std::vector<CMSat::Lit>{r, ~a, b});
}

void CNF::Exact_One(const std::vector<CMSat::Lit> &Lits) {
    clauses.push_back(Lits);
    for (auto it=Lits.begin();it!=Lits.end();it++){
        for (auto it2 = Lits.begin();it2!=it;it2++){
            clauses.push_back(std::vector<CMSat::Lit>{~(*it),~(*it2)});
        }
    }
}

void CNF::half_adder(const CMSat::Lit &a, const CMSat::Lit &b, const CMSat::Lit &result, const CMSat::Lit &carry) {
    XOR(a,b,result);
    AND(a,b,carry);
}

void CNF::full_adder(const CMSat::Lit &a, const CMSat::Lit &b, const CMSat::Lit &c, const CMSat::Lit &result,
                     const CMSat::Lit &carry) {
    CMSat::Lit r1(new_var(), false),
    c1 (new_var(), false),
    c2 (new_var(), false);
    half_adder(a,b,r1,c1);
    half_adder(r1,c,result,c2);
    OR(c1,c2,carry);
}

std::vector<CMSat::Lit> CNF::add(const std::vector<CMSat::Lit> &a, const std::vector<CMSat::Lit> &b) {
    assert(a.size()==b.size());
    std::vector<CMSat::Lit> r = new_vars(a.size());
    add(a,b,r);
    return r;
}

void CNF::add(const std::vector<CMSat::Lit> &a, const std::vector<CMSat::Lit> &b, const std::vector<CMSat::Lit> &r) {
    assert(a.size()==b.size() && a.size()==r.size());
    std::vector<CMSat::Lit> c = new_vars(a.size());
    for (int i=0; i<a.size(); i++){
        if (i) full_adder(a[i], b[i],c[i-1], r[i], c[i]);
        else half_adder(a[i], b[i], r[i], c[i]);
    }
    clauses.emplace_back(1,~c[a.size()-1]);
}

CMSat::Lit CNF::leq(std::vector<CMSat::Lit> a, const std::vector<CMSat::Lit> &b, bool free_result) {
    assert(a.size()==b.size());
    std::vector<CMSat::Lit> r,c;
    r = new_vars(a.size());
    c = new_vars(a.size());

    for (int i=0; i<a.size(); i++){
        a[i] = ~a[i];
    }
    for (int i=0; i<a.size(); i++){
        if (i) full_adder(a[i], b[i],c[i-1], r[i], c[i]);
        else full_adder(a[i], b[i], True, r[i], c[i]);
    }
    if(!free_result){
        clauses.emplace_back(1,c[a.size()-1]);
    }
    return c[a.size()-1];
}

void CNF::eq(const std::vector<CMSat::Lit> &a, const std::vector<CMSat::Lit> &b) {
    for (int i=0;i<a.size();i++){
        clauses.push_back(std::vector<CMSat::Lit>{~a[i],b[i]});
        clauses.push_back(std::vector<CMSat::Lit>{a[i],~b[i]});
    }
}

std::vector<CMSat::Lit> CNF::max(const std::vector<CMSat::Lit> &a, const std::vector<CMSat::Lit> &b) {
    assert(a.size()==b.size());
    std::vector<CMSat::Lit> r = new_vars(a.size());
    max(a,b,r);
    return r;
}

void CNF::max(const std::vector<CMSat::Lit> &a, const std::vector<CMSat::Lit> &b, const std::vector<CMSat::Lit> &r) {
    assert(a.size()==b.size() && a.size()==r.size());
    CMSat::Lit alb = leq(a, b, true);
    for (int i=0; i<a.size(); i++){
        OR(AND(~alb,a[i]), AND(alb,b[i]),r[i]);
    }
}

//std::vector<CMSat::Lit> CNF::increment(const std::vector<CMSat::Lit> &a, const CMSat::Lit &b = True) {
//    std::vector<CMSat::Lit> r(0);
//    for (int i=0; i<a.size(); i++){
//        r.push_back(CMSat::Lit(new_var(), false));
//    }
//    increment(a,b,r);
//    return r;
//}
//
//void CNF::increment(const std::vector<CMSat::Lit> &a, const CMSat::Lit &b, const std::vector<CMSat::Lit> &r) {
//    std::vector<CMSat::Lit> c(0);
//    for (int i=0; i<a.size(); i++){
//        c.push_back(CMSat::Lit(new_var(), false));
//    }
//    for (int i=0; i<a.size(); i++){
//        if (i) half_adder(a[i], c[i-1], r[i], c[i]);
//        else half_adder(a[i], b, r[i], c[i]);
//    }
// //    clauses.push_back(std::vector<CMSat::Lit>(1,~c[a.size()-1]));
//}

CMSat::Lit CNF::OR_Lits(const std::vector<CMSat::Lit> &a) {
    CMSat::Lit r(new_var(), false);
    OR_Lits(a,r);
    return r;
}

void CNF::OR_Lits(std::vector<CMSat::Lit> a, const CMSat::Lit &r) {
    for (int i=0;i<a.size();i++)
        clauses.push_back(std::vector<CMSat::Lit>{r,~a[i]});
    a.push_back(~r);
    clauses.push_back(a);
}

void CNF::add_clause(const std::vector<CMSat::Lit> &cl) {
    clauses.push_back(cl);
}


void CNF::to_file(const char *filename) const {
    FILE * out = fopen(filename,"w");
    fprintf(out,"p cnf %d %d\n",n_variables, clauses.size());
    fprintf(out,"c ind");
    for(int i = 0; i < ind_vs.size();i++){
        fprintf(out, " %d", ind_vs[i]+1);
    }
    fprintf(out," 0\n");
    for(auto it = clauses.begin();it != clauses.end(); it++){
        for(int i = 0; i < it->size(); i++){
            fprintf(out,"%d", ((*it)[i].sign()?-1:1)*((*it)[i].var()+1));
            fprintf(out, (i==it->size()-1?" 0\n":" "));
        }
    }
    fclose(out);
}

std::vector<CMSat::Lit> CNF::Solve(bool ind) {
    minisat = new CMSat::SATSolver;
    minisat -> new_vars(n_variables);
    minisat -> set_sampling_vars(&ind_vs);
    for (auto it = clauses.begin(); it!=clauses.end();it++){
        minisat->add_clause(*it);
    }
    CMSat::lbool ret = minisat ->solve();
    if (ret==CMSat::l_False){
        delete minisat;
        if (ind){
            return {};
        }
        return {False};
    }
    if (ind){
        std::vector<CMSat::Lit> result;
        int tmp=0;
        for (auto it = ind_vs.begin(); it != ind_vs.end(); it++){
            while(tmp<*it ){
                tmp++;
            }
            if (minisat->get_model()[tmp] == CMSat::l_True) result.emplace_back(tmp+1, false);
            else result.emplace_back(tmp+1, true);
        }
        delete minisat;
        return result;
    }
    delete minisat;
    return {True};
}

CNF::~CNF() {
}

void CNF::Enum_Sampling(const std::vector<uint32_t> & enum_set, int n_samples, std::map<std::vector<int>,int> & result,
                        std::vector<Callback> & data) {
    std::vector<ApproxMC::AppMC*> appmcs (1<<enum_set.size());
    std::vector<ApproxMC::SolCount> appmc_res (1<<enum_set.size());
    std::vector<CMSat::Lit>  additional_clauses (enum_set.size());

    long long tot_sol = 0;
    for (int i = 0; i < enum_set.size(); i++) {
        additional_clauses[i] = CMSat::Lit(enum_set[i],false);
    }
    for(int i = 0, _tmp; i < (1<<enum_set.size()); i++ ) {
        auto start = std::chrono::high_resolution_clock::now();
        appmcs[i] = new ApproxMC::AppMC;
        if (i > 0) {
            _tmp = i ^ (i - 1);
            for (auto it = additional_clauses.begin(); it != additional_clauses.end(); it++, _tmp >>= 1) {
                if (_tmp & 1)
                    (*it) = ~(*it);
            }
        }
        appmc_res[i] = Counting(*this, additional_clauses, appmcs[i]);

        if (appmc_res[i].cellSolCount == 0) {
            delete appmcs[i];
        }

        tot_sol += (1LL<<appmc_res[i].hashCount)*appmc_res[i].cellSolCount;

        auto stop = std::chrono::high_resolution_clock::now();

        std::cout << "ApproxMC: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << std::endl;
    }

    for(int i = 0, _tmp; i < (1<<enum_set.size()); i++ ) {
        if (appmc_res[i].cellSolCount > 0) {
            auto start = std::chrono::high_resolution_clock::now();
            _tmp = (1LL<<appmc_res[i].hashCount)*appmc_res[i].cellSolCount*n_samples/tot_sol+1;
            Sampling(_tmp, appmcs[i], appmc_res[i], &data[i]);
            auto stop = std::chrono::high_resolution_clock::now();

            std::cout << "UniGen: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl << std::endl;
        }
//        delete appmc;
    }
    for (int i = 0, _tmp, _c; i < (1<<enum_set.size()); i++ ) {
        if (appmc_res[i].cellSolCount <= 0) {
            continue;
        }

        _tmp = (1LL<<appmc_res[i].hashCount)*appmc_res[i].cellSolCount*n_samples/tot_sol+1;
        printf("appmc res: %d es..\n",_tmp);
        _c = 0;
        for (auto it = data[i].data.begin(); _c < _tmp; _c++, it++) {
            if(result.find(*it)!=result.end()){
                result[*it]++;
            } else {
                result[*it]=1;
            }
        }
    }
}

ApproxMC::SolCount CNF::Counting(const CNF &origin, const  std::vector<CMSat::Lit>  & additional_clauses,
                   ApproxMC::AppMC *appmc, int verbosity) {
    appmc -> set_verbosity(verbosity);
    appmc -> set_seed(rand());
    appmc -> new_vars(origin.n_variables);
    for (auto it = origin.clauses.begin(); it!=origin.clauses.end(); it++){
        appmc -> add_clause(*it);
    }
    for (auto it = additional_clauses.begin(); it!=additional_clauses.end(); it++){
        appmc -> add_clause({*it});
    }
    appmc ->set_projection_set(origin.ind_vs);
    appmc ->setup_vars();
    return appmc -> count();
}

void CNF::Sampling(int n_samples, ApproxMC::AppMC *appmc, const ApproxMC::SolCount & sol_count,
                   Callback * ptr_) {
    UniGen::UniG * unigen = new UniGen::UniG(appmc);
    unigen -> set_callback(Callback::callback, ptr_);
    unigen -> sample(&sol_count,  n_samples);

    delete unigen;
}










