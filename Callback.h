//
// Created by Yuanyuan Qi on 3/8/22.
//

#ifndef UNIPPM_CALLBACK_H
#define UNIPPM_CALLBACK_H

#include <vector>

class Callback {
public:
    int n;//, n_sol;
    std::vector<std::vector<int> >::iterator it;

    std::vector<std::vector<int> > data;
    Callback(int n_, int n_sol_);
    static void callback(const std::vector<int> & solution, void * _ptr);
};


#endif //UNIPPM_CALLBACK_H
