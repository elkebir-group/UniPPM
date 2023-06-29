//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_INPUT_H
#define UNIPPM_INPUT_H

#include <vector>
#include <string>

struct Input {
    int n,m; //n mutations, m samples
    int r;
    std::vector<std::vector<std::pair<double,double> > >  data;
    explicit Input(const std::string& filename, int r = 0);
};


#endif //UNIPPM_INPUT_H
