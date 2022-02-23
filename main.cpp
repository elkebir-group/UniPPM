#include <iostream>
#include <map>
#include <string>
#include <set>
#include <ctime>
#include <random>
//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

//po::options_description main_options = po::options_description("Main options");

#include "Solver.h"
#include "Likelihood.h"

using namespace std;
map<string,string> vm;

string input_file,output_file;
int n_samples,n_bits,n_intervals;
double approx_coef = -1, help_approx_coef;
long long seed;

void parse_argument(int argc,char * argv[]){
    string options,value;
    set<string> p_options({"-i","-o","-n","-N","-a","-A","-I","-s"});

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
            case 'a':
                approx_coef = stod(it->second);
                break;
            case 'A':
                help_approx_coef = stod(it->second);
                break;
            case 'I':
                n_intervals = stoi(it->second);
                break;
            case 's':
                seed = stoi(it->second);
        }
    }
    for (auto it = p_options.begin();it!=p_options.end();it++) {
        switch ((*it)[1]) {
            case 'i':
                cerr << "Input file" << endl;
                exit(1);
            case 'o':
                cerr << "Output file" <<endl;
                exit(1);
            case 'n':
                n_samples = 1000;
                break;
            case 'N':
                n_bits = 10;
                break;
            case 'a':
                approx_coef = 0.95;
                break;
            case 'A':
                if (approx_coef < 0)
                    help_approx_coef = 0.9;
                else
                    help_approx_coef = approx_coef*0.95;
                break;
            case 'I':
                n_intervals = 10;
                break;
            case 's':
                seed = time(0);
        }
    }
    srand(seed);
}

int main(int argc, char * argv[]) {
    parse_argument(argc,argv);
    Input_Reads raw_in(input_file.c_str());
    double l_a = 1e-9,r_a = 1-1e-9, t_alpha;
    vector<pair<int,int> > edges;
#define m_a ((l_a+r_a)/2)
    while (r_a-l_a>1e-7){
        Input transform_in(raw_in,m_a);
        Input_int in(transform_in,n_bits);
        AncestryGraph Gf(in);
        Solver tester(Gf,0);
        if(tester.attempt(tester.self_solver(),&edges))
            r_a = m_a;
        else
            l_a = m_a;
    }
    t_alpha = r_a;
    double llh;
    {
        Input transform_in(raw_in,t_alpha);
        Input_int in(transform_in, n_bits);
        Likelihood LLH(in,raw_in,n_intervals);
        llh = LLH.LLH(edges);
    }
    cout<<"estimate lowerbound of log likelihood:"<<llh<<endl;

    double target = llh/help_approx_coef, filtering = llh/approx_coef;
    l_a = 1e-9,r_a = 1-1e-9;
    double lower_bound_ll = 1e300;
    while (r_a-l_a>1e-7 && abs(lower_bound_ll-target)>1e-7){
        Input transform_in(raw_in,m_a);
        lower_bound_ll = lower_bound_llh(raw_in,transform_in,n_bits);
        if(lower_bound_ll>target)
            l_a = m_a;
        else
            r_a = m_a;
    }
    if (m_a < t_alpha){
        cout<<"approx alpha overrided by min alpha"<<endl;
    }
    else{
        t_alpha = m_a;
    }
#undef m_a
    cout<<"using "<< t_alpha << " as the final alpha."<<endl;

    Input transform_in(raw_in,pow(t_alpha,1.0/(raw_in.n*raw_in.m)));
    Input_int in(transform_in,n_bits);
    AncestryGraph Gf(in);
    Solver solver(Gf,0);
    Likelihood LLH(in,raw_in,n_bits);

    map<vector<pair<int,int> >,int> res;
    for (int i = 0; i < in.m; i++){
        for(int j = 0; j < in.n; j++){
            cout<<in.F_l[i][j]<<' ' <<in.F_u[i][j]<<endl;
        }
    }
    solver.sampling(solver.self_solver(),n_samples*2,res);

    for(auto i:res){
        cout <<"----------- sample:"<<i.second<<endl;
        for (auto j:i.first){
            cout<<j.first<<' '<<j.second<<endl;
        }
    }

    ostringstream sout;
    int unique = 0, total = 0;

    for(auto it = res.begin();it!=res.end();it++) {
        double ll = LLH.LLH(it->first);
        if (ll < filtering){
            continue;
        }
        sout<< "#" <<(it->first.size())<<" edges, tree "<<unique<<", "<<(it->second)<<" samples, log-likelihood:"<<ll<<'\n';
        unique ++ ;
        total += it->second;
        for(auto edge:it->first){
            sout<<edge.first<<' '<<edge.second<<'\n';
        }
    }
    ofstream fout(output_file);
    fout<<"#"<<total<<"tree, #"<<unique<<" unique trees"<<endl;
    fout<<sout.str();
    fout.close();
    return 0;
}
