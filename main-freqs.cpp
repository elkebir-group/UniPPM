//
// Created by Yuanyuan Qi on 3/21/22.
//

#include <iostream>
#include <map>
#include <string>
#include <set>
#include <random>

#include "Solver.h"

const double EPS(1e-9);

using namespace std;
map<string,string> vm;

string input_file,output_file;
int n_samples,n_bits;//,n_intervals;
//double approx_coef = -1, help_approx_coef;
long long seed;

int timeout;

void parse_argument(int argc,char * argv[]){
    string options,value;
    set<string> p_options({"-i","-o","-n","-N","-s","-T"});

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
            case 'n':
                n_samples = stoi(it->second);
                break;
            case 'N':
                n_bits = stoi(it->second);
                break;
            case 's':
                seed = stoi(it->second);
                break;
            case 'T':
                timeout = stoi(it->second);
                break;
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
            case 'n':
                n_samples = 20000;
                break;
            case 'N':
                n_bits = 20;
                break;
            case 's':
                seed = 1;//time(0);
                break;
            case 'T':
                timeout = 2000;
                break;
        }
    }
    srand(seed);
}

int main(int argc, char * argv[]) {
    parse_argument(argc,argv);
    Input raw_in(input_file.c_str());

    Input_int in(raw_in,n_bits);

    double l_a = EPS,r_a = 1-EPS, t_alpha;

    AncestryGraph Gf(in);
    Solver solver(Gf,timeout);

    map<vector<pair<int,int> >,int> res;
    solver.sampling(n_samples * 2, res);

//    for(auto i:res){
//        cout <<"----------- sample:"<<i.second<<endl;
//        for (auto j:i.first){
//            cout<<j.first<<' '<<j.second<<endl;
//        }
//    }

    int unique = 0, total = 0;

//    stringstream sout;

//    ofstream fout(output_file);
    cout<<"[UniPPM] calculating likelihood."<<endl;
    FILE * fout = fopen(output_file.c_str(),"w");
    for(auto it = res.begin();it!=res.end();it++) {

        fprintf(fout,"# %d edges, tree %d, %d sample\n",raw_in.n-1,unique,it->second);
        unique ++ ;
        total += it->second;
        for(auto edge:it->first){
            fprintf(fout,"%d %d\n",edge.first,edge.second);
        }
    }

    fprintf(fout,"# %d samples, # %d unique trees\n",total,unique);

    fclose(fout);
    return 0;
}
