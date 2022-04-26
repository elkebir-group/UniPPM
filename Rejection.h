//
// Created by Yuanyuan Qi on 4/20/22.
//

#ifndef UNIPPM_REJECTION_H
#define UNIPPM_REJECTION_H

#include "AncestryGraph.h"

#include <random>
#include <list>
#include <map>

class Rejection {
public:
    Rejection(const AncestryGraph & Gf);

    bool generate_tree(double up_sig=0.5, double lower_sig=0.00, double decay_rate=0.5);

    const std::vector<std::pair<int,int> > & get_tree();

    bool SC();

    void try_sample(int n_sample, std::map<std::vector<std::pair<int,int> >,int> & res);

    void sample(int n_sample, std::map<std::vector<std::pair<int,int> >,int> & res);
private:
    bool verify_SC(int v);
    bool random_tree(double sig);
    int random_walk(int start_v, double sig);
    const AncestryGraph & in;

    std::vector<bool> cover, in_walk;
    std::vector<int> r_seq;

    std::mt19937 mt;

    std::vector<int> rw; int rw_size;
    std::vector<int> rw_pointer;

    std::vector<std::pair<int,int> > tree; int tree_size,r;

    std::vector<int> F_verify;
    std::vector<std::list<int> > tree_out_d;
};


#endif //UNIPPM_REJECTION_H
