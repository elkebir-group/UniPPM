//
// Created by Yuanyuan Qi on 2/14/22.
//

#ifndef UNIPPM_CNF_H
#define UNIPPM_CNF_H

#include <vector>
#include <cryptominisat5/cryptominisat.h>
#include <approxmc/approxmc.h>
#include <unigen/unigen.h>

class CNF {
public:
    CNF();
    CNF(const CNF &);

    ~CNF();

    const CMSat::Lit True,False;

    int new_var(bool is_ind = false);
    std::vector<CMSat::Lit> new_vars(int, bool is_ind = false);

    CMSat::Lit AND(const CMSat::Lit &, const CMSat::Lit &);
    void AND(const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &);

    CMSat::Lit OR(const CMSat::Lit&, const CMSat::Lit &);
    void OR(const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &);

    CMSat::Lit XOR(const CMSat::Lit&, const CMSat::Lit &);
    void XOR(const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &);

    void Exact_One(const std::vector<CMSat::Lit> &);

    void half_adder(const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &);

    void full_adder(const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &, const CMSat::Lit &);

    std::vector<CMSat::Lit> add(const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &);
    void add(const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &);

    CMSat::Lit leq(std::vector<CMSat::Lit>, const std::vector<CMSat::Lit> &, bool free_result = false);

    void eq(const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &);

    std::vector<CMSat::Lit> max(const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &);
    void max(const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &, const std::vector<CMSat::Lit> &);

//    std::vector<CMSat::Lit> increment(const std::vector<CMSat::Lit> &, const CMSat::Lit &);
//    void increment(const std::vector<CMSat::Lit> &, const CMSat::Lit &, const std::vector<CMSat::Lit> &);

    CMSat::Lit OR_Lits(const std::vector<CMSat::Lit> &);
    void OR_Lits(std::vector<CMSat::Lit>, const CMSat::Lit &);

    void add_clause(const std::vector<CMSat::Lit> &);

    void to_file(const char *)const;

    std::vector<CMSat::Lit> Solve(bool =true);

    int Counting(bool = true, int = 1);
    void Sampling(int, int = 1, bool = false);

    const std::vector< std::vector<int> > &Results;


private:
    std::vector< std::vector<CMSat::Lit> > clauses;
    int n_variables;
    std::vector<uint32_t> ind_vs;
    CMSat::SATSolver *minisat;

    ApproxMC::AppMC *appmc;
    ApproxMC::SolCount sol_count;
    UniGen::UniG *unigen;

    std::vector< std::vector<int> > callbackdata;
};

void callback(const std::vector<int> &, void *);

#endif //UNIPPM_CNF_H
