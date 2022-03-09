//
// Created by Yuanyuan Qi on 3/8/22.
//

#include "Callback.h"


void Callback::callback(const std::vector<int> & solution, void* ptr_data) {
    Callback * callbackdata = (Callback *)ptr_data;
    auto it1 = solution.begin();
//    auto it2 = callbackdata->it->begin();
    auto it2 = callbackdata->data[callbackdata->index].begin();
    while (it1!=solution.end()){
        if (*it1 > 0){
            *(it2++) = *it1;
        }
        it1++;
    }
//    callbackdata->it ++;
    callbackdata->index++;
}

Callback::Callback(int n_, int n_sol_):n(n_), data(n_sol_, std::vector<int>(n_)) {//, n_sol(n_sol_)
//    it = data.begin();
    index = 0;
}
