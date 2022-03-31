//
// Created by Yuanyuan Qi on 3/30/22.
//

#include <iostream>
#include <queue>

#include "Likelihood.h"

using namespace std;

map<string,string> vm;
string input_file,output_file;
int n_bits,n_intervals;
bool mul;

void parse_argument(int argc,char * argv[]){
    string options,value;
    set<string> p_options({"-i","-o","-N","-I","-M"});

    for (int i = 1; i < argc; i++){
        options = argv[i];
        if(options == "-h" || options == "--help"){
            exit(0);
        }
        value = argv[i+1];
        vm [options] = value;
        i++;
    }
    for (auto it = vm.begin();it!=vm.end();it++){
        p_options.erase(it->first);
        switch (it->first[1]) {
            case 'i':
                input_file = it -> second;
                break;
            case 'o':
                output_file = it -> second;
                break;
            case 'N':
                n_bits = stoi(it->second);
                break;
            case 'I':
                n_intervals = stoi(it->second);
                break;
            case 'M':
                mul = stoi(it->second)!=0;
        }
    }

    for (auto it = p_options.begin();it!=p_options.end();it++) {
        switch ((*it)[1]) {
            case 'i':
                cerr << "Input file" << endl;
//                input_file = "input.txt";
                exit(1);
            case 'o':
                cerr << "Output file" <<endl;
//                output_file = "tmp.txt";
                exit(1);
            case 'N':
                n_bits = 20;
                break;
            case 'I':
                n_intervals = 100;
                break;
            case 'M':
                mul = false;
        }
    }
}



void get_edge_set(vector<int> prufer_seq, int root, list<vector<pair<int,int> > > & store){
    int n = prufer_seq.size()+2;
    vector<vector<int> > link_list(n);
    set<int> remaining;
    map<int,int> prufer_rem;
    for (auto i:prufer_seq){
        if (prufer_rem.find(i)!=prufer_rem.end()) prufer_rem[i]++;
        else prufer_rem[i]=1;
    }
    for (int i = 0; i < n; i++){
        if(prufer_rem.find(i)==prufer_rem.end())
            remaining.emplace_hint(remaining.end(),i);
    }

    for (int i = 0; i < n-2; i++){
        auto v1 = *remaining.begin();
//        edge_set.emplace_back(v1,prufer_seq[i]);
        link_list[v1].push_back(prufer_seq[i]);
        link_list[prufer_seq[i]].push_back(v1);
        remaining.erase(v1);
        prufer_rem[prufer_seq[i]]--;
        if(!prufer_rem[prufer_seq[i]]){
            prufer_rem.erase(prufer_seq[i]);
            remaining.insert(prufer_seq[i]);
        }
    }
    link_list[*remaining.begin()].push_back(*remaining.rbegin());
    link_list[*remaining.rbegin()].push_back(*remaining.begin());

    vector<pair<int,int> >  arc_set(n-1);
    queue<int> bfs;
    vector<bool> visited(n,false);

    for (int i = 0; i < n; i++){
        if (root < 0 || i == root) {
            fill(visited.begin(),visited.end(), false);
            auto it = arc_set.begin();
            bfs.push(i);
            visited[i] = true;
            while (!bfs.empty()) {
                auto tmp = bfs.front();
                bfs.pop();
                for (auto k: link_list[tmp]) {
                    if (!visited[k]) {
                        bfs.push(k);
                        visited[k] = true;
                        it->first = tmp;
                        it->second = k;
                        it++;
                    }
                }
            }
            store.push_back(arc_set);
        }
    }
}

bool prufer_seq_increment(vector<int> & prufer_seq){
    int n = prufer_seq.size()+2;
    prufer_seq[0]++;
    for(int i = 1; i < n-2; i++){
        prufer_seq[i]+=prufer_seq[i-1]/n;
        prufer_seq[i-1]%=n;
    }
    return prufer_seq[n-3] < n;
}

int main(int argc,char * argv[]){
    parse_argument(argc,argv);
    Input_Reads raw_in(input_file.c_str());
    Input llhrange(raw_in,1,mul);
    Input_int llhrange_int(llhrange,n_bits);
    Likelihood llh(llhrange_int,raw_in,n_intervals,mul);

    list<vector<pair<int,int> > > all_trees;
    vector<int> prufer_seq(llhrange_int.n-2,0);
    prufer_seq[0]=-1;
    int i;
    for(i=0;prufer_seq_increment(prufer_seq);i++){
        get_edge_set(prufer_seq,llhrange_int.r,all_trees);
    }

    i=0;
    auto fout = fopen(output_file.c_str(),"w");
    for(auto t:all_trees){
        fprintf(fout,"# %d edges, tree %d, log-likelihood: %lf\n",(int)(t.size()),i,llh.LLH(t));
        for(auto edge:t){
            if (mul && edge.first == raw_in.n){
                fprintf(fout,"%d %d\n",-1,edge.second);
            }
            else {
                fprintf(fout, "%d %d\n", edge.first, edge.second);
            }
        }
        i++;
    }

    return 0;
}