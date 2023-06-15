//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_ANCESTRYGRAPH_H
#define UNIPPM_ANCESTRYGRAPH_H

#include <vector>
#include <list>

class Input;

struct AncestryGraph {
    const Input & data;
    explicit AncestryGraph(const Input &data);
    std::vector<std::list<int> > possible_children;
    std::vector<std::list<int> > possible_parent;
    std::vector<std::vector<int> > arc_set_index;
    std::vector<std::pair<int,int> > arc_set;
};


#endif //UNIPPM_ANCESTRYGRAPH_H
