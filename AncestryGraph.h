//
// Created by Yuanyuan Qi on 2/14/22.
//

#ifndef UNIPPM_ANCESTRYGRAPH_H
#define UNIPPM_ANCESTRYGRAPH_H

#include "Input.h"
#include <vector>


class AncestryGraph {
public:
    AncestryGraph(const Input_int &in);
    const Input_int &In;

    const std::vector<std::pair<int,int> > & E;

    const std::vector<int> & ind(int i) const;
    const std::vector<int> & outd(int i) const;

    const bool &is_single;

private:
    std::vector<std::vector<int> > out_d;
    std::vector<std::vector<int> > in_d;
    std::vector< std::pair<int,int> > edge_set;

    bool Is_single;
};


#endif //UNIPPM_ANCESTRYGRAPH_H
