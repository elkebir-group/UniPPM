//
// Created by Yuanyuan Qi on 2/14/22.
//

#include "CNF.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "Solver.h"
//#include <pthread.h>
#include <thread>

CNF::CNF():True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), //appmc(NULL), unigen(NULL),
n_variables(1), //Results(callbackdata),
clauses({std::vector<CMSat::Lit>{ CMSat::Lit(0, false)} }),
ind_vs()
{
}

CNF::CNF(const CNF & cp):True(CMSat::Lit(0, false)),False(CMSat::Lit(0, true)),
minisat(NULL), //appmc(NULL), unigen(NULL),
clauses(cp.clauses), n_variables(cp.n_variables), ind_vs(cp.ind_vs)//Results(callbackdata)
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
    fprintf(out,"p cnf %d %ld\n",n_variables, clauses.size());
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

void CNF::Counting(const CNF &origin, const  std::list<CMSat::Lit>  & additional_clauses,
                   ApproxMC::AppMC * appmc, ApproxMC::SolCount &res, int verbosity) {
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
    res = appmc -> count();
}

void callback(const std::vector<int> & solution, void* ptr_data) {
    auto * callbackdata = (std::list<std::vector<int> > *)ptr_data;
    auto it1 = solution.begin();
//    auto it2 = callbackdata->it->begin();
    std::vector<int> v2;
//    auto it2 = v2.begin();
    while (it1!=solution.end()){
        if (*it1 > 0){
            v2.push_back(*it1);
        }
        it1++;
    }
//    callbackdata->it ++;
    callbackdata->emplace_back(v2);
}


void CNF::Sampling(int n_samples, ApproxMC::AppMC *appmc, const ApproxMC::SolCount & sol_count,
                   std::list<std::vector<int> > * ptr_) {
    UniGen::UniG * unigen = new UniGen::UniG(appmc);
    unigen -> set_callback(callback, ptr_);
    unigen -> sample(&sol_count,  n_samples);
    delete unigen;
}

void mytimeout(ApproxMC::AppMC * & appmc, int ms){
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    if (appmc) appmc->signal_stop();
}

void CNF::UniPPM_Preparing(int timeout, /*int rec_t,*/ int rec_step, int force_layer, Solver *ptr,
                           std::list<CMSat::Lit> &additional_clauses,
                           CNF::rec_node *root, const std::string &info_tag) {
    if(root == nullptr){
        root = &this->root;
    }

    if (rec_step >= force_layer) {

        root->appmc = new ApproxMC::AppMC;
        std::cout << "[UniPPM][" << info_tag << "] Estimating solutions.." << std::endl;

        std::thread count_t(Counting, std::ref(*this), std::ref(additional_clauses), root->appmc,
                            std::ref(root->res), 1);

        std::thread timer(mytimeout, std::ref(root->appmc), timeout);
        timer.detach();

        count_t.join();
    }
    else {
        std::cout << "[UniPPM][" << info_tag << "] forced layer.." << std::endl;
    }

    if (root->res.hashCount < 0x7fffffff){
        root->count = (1LL<<root->res.hashCount)*root->res.cellSolCount;
//        if (root->count > rec_t){
//            std::cout << "[UniPPM][" << info_tag << "] Estimated " << root->count << " solutions,";
//            split = true;
//        }
//        else {
            std::cout << "[UniPPM][" << info_tag << "] Estimated " << root->count << " solutions" << std::endl;
//        }
        root->split = false;
    }
    else
        root->split = true;

    delete root->appmc;
    root->appmc = nullptr;

    if(root->split)
    {
//        root->appmc->signal_stop();
        root->n=ptr->CNF_recursive_sets[rec_step].size();
        std::cout << "[UniPPM][" << info_tag << "] Too hard. Split into "<<root->n<<" branches."<<std::endl;
        root->splits = new rec_node[root->n];
        root->count = 0;
        for(int i = 0;i < root->n; i++){
            additional_clauses.emplace_back(ptr->CNF_recursive_sets[rec_step][i],false);
            UniPPM_Preparing(timeout,rec_step+1,force_layer,ptr,additional_clauses,& root->splits[i],
                             info_tag+"."+ std::to_string(i));
            additional_clauses.pop_back();
            root->count += root->splits[i].count;
        }
    }
}

void CNF::UniPPM_Sampling(int n_samples,int rec_step,Solver *ptr, std::list<CMSat::Lit> &additional_clauses,
                          std::list<std::vector<int> > &data, CNF::rec_node *root,
                          const std::string &info_tag) {
    if (root == nullptr){
        root = &this->root;
    }
    if (!root->count){
        return;
    }
    if (!root->split){
        root->appmc = new ApproxMC::AppMC;
        std::cout << "[UniPPM][" << info_tag << "] sampling with unigen: ("
                  << n_samples << " trees from " << root->count << " solutions)." << std::endl;
        Counting(*this, additional_clauses,root->appmc,root->res);
        Sampling(n_samples,root->appmc,root->res,&data);
        delete root->appmc;
    }
    else{
        std::cout<< "[UniPPM][" << info_tag << "] sampling (branching)." <<std::endl;
        for(int i = 0; i < root->n; i++){
            additional_clauses.emplace_back(ptr->CNF_recursive_sets[rec_step][i],false);
            UniPPM_Sampling(n_samples*root->splits[i].count/root->count+1,rec_step+1,ptr,additional_clauses,
                            data,&root->splits[i],
                            info_tag+"."+ std::to_string(i));
            additional_clauses.pop_back();
        }
    }
}


