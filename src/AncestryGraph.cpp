//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "AncestryGraph.h"
#include "Input.h"

AncestryGraph::AncestryGraph(const Input &data):
data(data),possible_children(data.n), possible_parent(data.n),
arc_set_index(data.n,std::vector<int>(data.n,-1))
{
    int count=0;
    bool tmp_flag;
    for (auto i=0; i< data.n;++i){
        for(auto j=0;j< data.n;++j){
            if(i==j or j==data.r) continue;
            tmp_flag= true;
            for (auto _sample=0;_sample<data.m;_sample++){
                if(data.data[_sample][i] < data.data[_sample][j]){ //i can't be ancestor of j
                    tmp_flag= false;
                    break;
                }
            }
            if (tmp_flag){
                count+=1;
                possible_children[i].push_back(j);
                possible_parent[j].push_back(i);
            }
        }
    }
    arc_set.resize(count);
    for (auto i=0,_count=0; i< data.n;++i){
        for(auto j:possible_children[i]){
            arc_set_index[i][j]=_count;
            arc_set[_count]={i,j};
            ++_count;
        }
    }
}
