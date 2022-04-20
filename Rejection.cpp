//
// Created by Yuanyuan Qi on 4/20/22.
//

#include "Rejection.h"


bool Rejection::random_tree(double sig) {
//    int init_v = in.In.r;
    int is_f=0,c_rw_v;
    tree_size=0;
    std::fill(cover.begin(),cover.end(), false);
    std::shuffle(r_seq.begin(),r_seq.end(),mt);
    for (auto it=r_seq.begin();it!=r_seq.end();it++) {
        if (cover[*it]) continue;
        c_rw_v = random_walk(*it, sig);
        if (c_rw_v < 0) {
            is_f++;
            if (is_f > 1) return false;
        }
        for (auto i=rw_size-1; i>=0; i--) {
            if(c_rw_v >= 0) {
                tree[tree_size++] = {c_rw_v, rw[i]};
            }
            else {
                r = rw[i];
            }
            cover[rw[i]]=true;
            c_rw_v = rw[i];
        }
    }
    return true;
}

int Rejection::random_walk(int start_v, double sig) {
    std::fill(in_walk.begin(),in_walk.end(), false);
//    rw.clear();
    rw_size = 1;
    rw[0] = (start_v);
    rw_pointer[start_v] = 0;
    in_walk[start_v]= true;
    int c_v,n_v;
    while (true){
        c_v = rw[rw_size-1];
        auto ctmp = rand()/(double)RAND_MAX;
        if (in.ind(c_v).size()==0 || (ctmp < sig) )
            return -1;
        n_v = in.ind(c_v)[rand()%in.ind(c_v).size()];
        if (in_walk[n_v]){
//            std::fill(rw.begin()+rw_pointer[n_v]+1,rw.begin()+rw_size, false);
            for (auto it = rw.begin()+rw_pointer[n_v]+1;it!=rw.begin()+rw_size;it++){
                in_walk[*it]=false;
            }
            rw_size = rw_pointer[n_v]+1;
            continue;
        }
        if (cover[n_v]){
            return n_v;
        }
        rw_pointer[n_v] = rw_size;
        rw[rw_size++]=n_v;
        in_walk[n_v] = true;
    }
}

const std::vector<std::pair<int, int> > &Rejection::get_tree() {
    return tree;
}

bool Rejection::generate_tree(double up_sig, double lower_sig, double decay_rate) {
    double sig = up_sig;
    while(!random_tree(sig)){
        sig = (sig - lower_sig)*decay_rate+lower_sig;
    }
    return false;
}

bool Rejection::SC() {
    tree_out_d.clear();
    tree_out_d.resize(in.In.n);
//    std::vector<std::list<int> > tree_in_d(in.In.n);
    for(auto e:tree){
        tree_out_d[e.first].push_back(e.second);
//        tree_in_d[e.second].push_back(e.first);
    }
    std::fill(F_verify.begin(),F_verify.end(),0);
    return verify_SC(r);
}

bool Rejection::verify_SC(int v) {
#define F_mc(A,B) F_verify[(A)*in.In.m+(B)]
    for (auto n_v:tree_out_d[v]){
        if(!verify_SC(n_v))
            return false;
    }
    for (int i = 0; i < in.In.m; i++){
        for (auto n_v:tree_out_d[v]){
            F_mc(i,v) += F_mc(i,n_v);
        }
        if (F_mc(i,v) > in.In.F_u[i][v])
            return false;
        F_mc(i,v) = std::max(F_mc(i,v),in.In.F_l[i][v]);
    }
    return true;
#undef F_mc
}

void Rejection::try_sample(int n_sample, std::map<std::vector<std::pair<int, int>>, int> &res) {
    for (int i = 0; i < n_sample; i++){
//        do{
            generate_tree();
//        }while(!SC());
        if(!SC()) continue;
        auto it = res.find(tree);
        if (it !=res.end()){
            it->second ++;
        }
        else {
            res[tree]=1;
        }
    }
}

void Rejection::sample(int n_sample, std::map<std::vector<std::pair<int, int> >, int> &res) {
    std::cout<<"[UniPPM] sampling with rejection sampling.." << std::endl;
    for (int i = 0; i < n_sample; i++){
        std::cout<<"[UniPPM] sampling " <<(i+1) <<"th tree."<< std::endl;
        do{
            generate_tree();
        }while(!SC());
        if(!SC()) continue;
        auto it = res.find(tree);
        if (it !=res.end()){
            it->second ++;
        }
        else {
            res[tree]=1;
        }
    }
}

Rejection::Rejection(const AncestryGraph &Gf):in(Gf),mt(rand()),
cover(Gf.In.n),in_walk(Gf.In.n),r_seq(Gf.In.n),
rw(Gf.In.n),rw_pointer(Gf.In.n),
tree(Gf.In.n-1),
F_verify(Gf.In.m*Gf.In.n){
    for (auto it=r_seq.begin();it!=r_seq.end();it++){
        *it = (int)(it-r_seq.begin());
    }
}


