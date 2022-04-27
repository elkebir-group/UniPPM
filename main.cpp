#include <iostream>
#include <map>
#include <string>
#include <set>
#include <random>
#include <chrono>
#include <iomanip>

//po::options_description main_options = po::options_description("Main options");

#include "Solver.h"
#include "Likelihood.h"
#include "Rejection.h"

const double EPS(1e-9);

using namespace std;
map<string,string> vm;

string input_file,output_file;
int n_samples,n_bits,n_intervals;
double approx_coef = -1, help_approx_coef;
long long seed;
bool mul;
string add_data;

int force_layer;//,timeout;
string cnf_file;
int rej_att,rej_thr;

void parse_argument(int argc,char * argv[]){
    string options,value;
    set<string> p_options({"-i","-o","-n","-N","-a","-A","-I","-s","-R","-M","-C","-F","-J","-j"});

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
                break;
            case 'R':
                force_layer = stoi(it->second);
                break;
            case 'M':
                mul = stoi(it->second)!=0;
                break;
            case 'C':
                add_data = it->second;
                break;
            case 'F':
                cnf_file = it->second;
                break;
            case 'J':
                rej_att = stoi(it->second);
                break;
            case 'j':
                rej_thr = stoi(it->second);
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
                seed = 1;//time(0);
                break;
            case 'R':
                force_layer = 2;
                break;
            case 'M':
                mul = false;
                break;
            case 'C':
                add_data = "";
                break;
            case 'F':
                cnf_file = "";
                break;
            case 'J':
                rej_att = -1;
                break;
            case 'j':
                rej_thr = -1;
        }
    }
    srand(seed);
}

int main(int argc, char * argv[]) {
    auto start = chrono::high_resolution_clock::now();
    parse_argument(argc,argv);
    Input_Reads raw_in(input_file.c_str());
    AdditionalData additionalData (raw_in.n, add_data);

    double l_a = EPS,r_a = 1-EPS, t_alpha;
    vector<pair<int,int> > edges;
#define m_a ((l_a+r_a)/2)
    while (r_a-l_a>EPS){
        cout<<"[UniPPM] attempt solving with alpha = " << setprecision(11) << m_a << endl;
        Input transform_in(raw_in,m_a,mul);
        Input_int in(transform_in,n_bits);
        AncestryGraph Gf(in);
        Solver tester(Gf);
//        tester.extract_CNF("tmp.cnf");
        tester.add_additional_constraints(additionalData);
        bool att = tester.attempt(tester.self_solver(),&edges);

        if(att)
            r_a = m_a;
        else
            l_a = m_a;
    }
    t_alpha = r_a;
    double llh;
    {
        Input transform_in(raw_in,t_alpha,mul);
        Input_int in(transform_in, n_bits);
        Likelihood LLH(in,raw_in,n_intervals,mul);
        llh = LLH.LLH(edges);
    }


    double target = llh/help_approx_coef, filtering = llh/approx_coef;
    cout<<"[UniPPM] estimate lowerbound of log likelihood: "<<llh<<", target llh: "<<target<<endl;

    l_a = EPS,r_a = 1-EPS;
    double lower_bound_ll = 1e300;
    while (r_a-l_a>EPS && abs(lower_bound_ll-target)>EPS){
        Input transform_in(raw_in,m_a,mul);
        lower_bound_ll = lower_bound_llh(raw_in,transform_in,n_bits,mul);
        if(lower_bound_ll>target)
            l_a = m_a;
        else
            r_a = m_a;
    }
    if (m_a < t_alpha){
        cout<<"[UniPPM] warning: given help_approx_coef is too large, overrided by min_alpha required."<<endl;
    }
    else{
        t_alpha = m_a;
    }
#undef m_a
    cout<<"[UniPPM] using "<< t_alpha << " as the final alpha."<<endl;

    Input transform_in(raw_in,pow(t_alpha,1.0/(raw_in.n*raw_in.m)),mul);
    transform_in.Show(cout);
    Input_int in(transform_in,n_bits);
    AncestryGraph Gf(in);
    Solver solver(Gf, force_layer);
    solver.add_additional_constraints(additionalData);
    Likelihood LLH(in,raw_in,n_bits,mul);

    map<vector<pair<int,int> >,int> res;

    Rejection rej(Gf);
    if (rej_att < 0){
        rej_att = int(pow(in.n,(in.n-2)/4.0))*20;
    }
    if (rej_thr < 0 ){
        rej_thr = 7;
    }
    rej.try_sample(rej_att,res);
    int tmp=0;
    for(auto it = res.begin();it!=res.end();it++){
        tmp+=it->second;
    }
    cout<<"[UniPPM] get "<< tmp << " sols from "<< rej_att << " rej attempt."<<endl;
    res.clear();
    if(tmp>=rej_thr){
        rej.sample(n_samples*2, res);
    }
    else {
        solver.sampling(n_samples * 2, res);
    }

//    for(auto i:res){
//        cout <<"----------- sample:"<<i.second<<endl;
//        for (auto j:i.first){
//            cout<<j.first<<' '<<j.second<<endl;
//        }
//    }

    int unique = 0, total = 0;

//    stringstream sout;

//    ofstream fout(output_file);
    std::vector<string> out_tag(in.n,"");
    for(int i=0;i<raw_in.n;i++){
        if(out_tag[raw_in.cl[i]].empty()) {
            out_tag[raw_in.cl[i]] += to_string(i);
        }
        else{
            out_tag[raw_in.cl[i]] += "_";
            out_tag[raw_in.cl[i]] += to_string(i);
        }
    }
    if(mul){
        out_tag[in.n-1] = "-1";
    }

    auto e1 = chrono::high_resolution_clock::now();
    chrono::duration<double,milli> dur = (e1-start);
    cout<<"[UniPPM] sampling finished, takes total of "<< dur.count()/1000.0 <<" seconds."<<endl;

    cout<<"[UniPPM] calculating likelihood."<<endl;
    FILE * fout = fopen(output_file.c_str(),"w");
    for(auto it = res.begin();it!=res.end();it++) {
        double ll = LLH.LLH(it->first);
        if (ll < filtering){
            continue;
        }
        fprintf(fout,"# %d edges, tree %d, %d sample, log-likelihood: %lf\n",(int)(it->first.size()),unique,it->second,ll);
        unique ++ ;
        total += it->second;
//        for(auto edge:it->first){
//            if (mul && edge.first == raw_in.n){
//                fprintf(fout,"%d %d\n",-1,edge.second);
//            }
//            else {
//                fprintf(fout, "%d %d\n", edge.first, edge.second);
//            }
//        }
        for(auto edge:it->first){
            fprintf(fout,"%s %s\n",out_tag[edge.first].c_str(),out_tag[edge.second].c_str());
        }
    }

    fprintf(fout,"# %d samples, # %d unique trees\n",total,unique);
    auto e2 = chrono::high_resolution_clock::now();
    chrono::duration<double,milli> dur2 = (e2-e1);
    cout<<"[UniPPM] likelihood calculations take "<< dur2.count()/1000.0 <<" seconds."<<endl;
    chrono::duration<double,milli> dur_t = (e2-start);
    cout<<"[UniPPM] total: "<< dur_t.count()/1000.0 << " seconds."<<endl;

    fclose(fout);
    return 0;
}
