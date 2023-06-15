//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_HASHING_H
#define UNIPPM_HASHING_H


#include <vector>

struct Hashing {
    std::vector<int> hash_f_coef;
    int n;
    Hashing(int n);
    void generate();
};


#endif //UNIPPM_HASHING_H
