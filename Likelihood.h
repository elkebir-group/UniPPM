//
// Created by Yuanyuan Qi on 2/14/22.
//

#ifndef UNIPPM_LIKELIHOOD_H
#define UNIPPM_LIKELIHOOD_H

#include "AncestryGraph.h"


class Likelihood {
public:
    Likelihood(const AncestryGraph &Gf);

    ~Likelihood();

    const AncestryGraph &Gf;
private:

};


#endif //UNIPPM_LIKELIHOOD_H
