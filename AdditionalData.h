//
// Created by Yuanyuan Qi on 4/4/22.
//

#ifndef UNIPPM_ADDITIONALDATA_H
#define UNIPPM_ADDITIONALDATA_H

#include <string>
#include "CNF.h"
#include "AncestryGraph.h"
#include <vector>
#include <map>

class AdditionalData {
    int n_c,n;
    short ** data;
    short * ptr;
public:
    AdditionalData(const int &n, const std::string & file_name);

    void Enforce_constraints(CNF * F, const std::vector<std::vector<CMSat::Lit> > & relation) const;

    ~AdditionalData();
};


#endif //UNIPPM_ADDITIONALDATA_H
