//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "Hashing.h"

#include <random>

Hashing::Hashing(int n):n(n),hash_f_coef(n) {
}

void Hashing::generate() {
    for (auto &a:hash_f_coef)
        a=rand()%2;
}
