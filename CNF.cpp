//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "CNF.h"
#include <random>
#include <cassert>

CNF::CNF():True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), appmc(NULL), unigen(NULL), n_variables(1), Results(callbackdata),
clauses({std::vector<CMSat::Lit>{ CMSat::Lit(0, false)} }),
ind_vs()
{
}

CNF::CNF(const CNF & cp):True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), appmc(NULL), unigen(NULL),
clauses(cp.clauses), n_variables(cp.n_variables), ind_vs(cp.ind_vs), Results(callbackdata)
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
    clauses.push_back(std::vector<CMSat::Lit>(1,~c[a.size()-1]));
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
        clauses.push_back(std::vector<CMSat::Lit>(1,c[a.size()-1]));
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
            if (minisat->get_model()[tmp] == CMSat::l_True) result.push_back(CMSat::Lit(tmp+1, false));
            else result.push_back(CMSat::Lit(tmp+1, true));
        }
        delete minisat;
        return result;
    }
    delete minisat;
    return {True};
}

int CNF::Counting(bool sampling , int verbosity) {
    if (appmc != NULL){
        delete appmc;
        delete unigen;
        unigen = NULL;
    }
    appmc = new ApproxMC::AppMC;
    appmc ->set_verbosity(verbosity);

    appmc -> new_vars(n_variables);
    for (auto it = clauses.begin(); it!=clauses.end();it++){
        appmc -> add_clause(*it);
    }

    appmc -> set_projection_set (ind_vs);
    appmc -> setup_vars();

    if (sampling) {
        callbackdata.clear();
        unigen = new UniGen::UniG(appmc);
        unigen -> set_callback(callback, &callbackdata);
    }

    sol_count = appmc -> count();

    return (1<<sol_count.hashCount)*sol_count.cellSolCount;
}

void CNF::Sampling(int n_samples, int verbosity, bool keep) {
    if (unigen == NULL){
        Counting(true, verbosity);
    }

    unigen -> sample(&sol_count, n_samples);

    if (!keep) {
        delete unigen;
        unigen = NULL;
    }
}

CNF::~CNF() {
    delete appmc;
    delete unigen;
}

void callback(const std::vector<int> & solution, void* ptr_data) {
    std::vector< std::vector<int> > *callbackdata = (std::vector< std::vector<int> > *)ptr_data;
    callbackdata->push_back(solution);
}





