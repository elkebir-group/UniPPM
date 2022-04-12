//
// Created by Yuanyuan Qi on 2/11/22.
//

#include "Input.h"

#include <boost/math/special_functions/beta.hpp>

#include <cmath>
#include <fstream>



Input::Input(const char * filename, bool multi) : m(M), n(N), r(R), F_l(F_lower), F_u(F_upper) {
    std::ifstream fin(filename);
    fin >> M >> N;
    if (multi) N++;

    F_lower = new double *[m];
    F_upper = new double *[m];
    ptr = new double[m*n*2];

    for (int i = 0; i < m; i++) {
        F_lower[i] = ptr + i * n;
        F_upper[i] = ptr + (i + m) * n;
        for (int j = 0; j < n-(multi?1:0); j++){
            fin >> F_lower[i][j] >> F_upper[i][j];
        }
        if (multi) {
            F_lower[i][n - 1] = 1-1e9;
            F_upper[i][n - 1] = 1-1e9;
        }
    }
    fin >> R;
    if (multi) R = n-1;
}

Input::~Input() {
    delete [] ptr;
    delete [] F_upper;
    delete [] F_lower;
}

Input::Input(const Input_Reads &In, const double &alpha, bool multi): m(M), n(N), r(R), F_l(F_lower), F_u(F_upper) {
    M = In.m;
//    N = In.n;
    N = In.n_cluster;
    if (multi) N++;

    F_lower = new double *[m];
    F_upper = new double *[m];
    ptr = new double[m*n*2];

    for (int i = 0; i < m; i++) {
        F_lower[i] = ptr + i * n;
        F_upper[i] = ptr + (i + m) * n;
        for (int j = 0; j < n-(multi?1:0); j++){
//            beta_distribution<> dis(In.var[i][j]+0.5, In.ref[i][j]+0.5)
            if (alpha>=1){
                F_lower[i][j] = 0;
                F_upper[i][j] = 1-1e-9;
            }
            else {
                F_lower[i][j] = In.ave_var[i][j] ? boost::math::ibeta_inv(In.ave_var[i][j] + 0.5, In.ave_ref[i][j] + 0.5,
                                                                      (1 - alpha) / 2) : 0;
                F_upper[i][j] = In.ave_ref[i][j] ? boost::math::ibeta_inv(In.ave_var[i][j] + 0.5, In.ave_ref[i][j] + 0.5,
                                                                      (1 + alpha) / 2) : (1-1e-9);
            }
        }
        if (multi) {
            F_lower[i][n - 1] = 0.5;
            F_upper[i][n - 1] = 0.5;
        }
    }

    R = In.r;
    if (multi) R = n-1;
}

void Input::Show(std::ostream & out) {
#define tout (out<<"[UniPPM] ")
    out << "---------------[UniPPM]----------------" <<std::endl;
    tout << "m:"<< m <<" n:"<< n << std::endl;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++){
            tout<< "["<<i<<"]["<<j<<"]:F_lower="<<F_lower[i][j]<<", F_upper="<<F_upper[i][j]<<std::endl;
        }
    out << "---------------[UniPPM]----------------" <<std::endl;
#undef tout
}


Input_int::Input_int(const Input & In, const int &N_bits):m(In.m),n(In.n),r(In.r),Nn(N_bits),n_bits(Nn),F_l(F_lower),F_u(F_upper) {
    F_lower = new int *[m];
    F_upper = new int *[m];
    ptr = new int[m * n * 2];

    for (int i = 0; i < m; i++) {
        F_lower[i] = ptr + i * n;
        F_upper[i] = ptr + (i + m) * n;
        for (int j = 0; j < n; j++) {
            F_lower[i][j] = int(In.F_l[i][j] * (1 << n_bits));
            F_upper[i][j] = int(In.F_u[i][j] * (1 << n_bits));
        }
    }
}

Input_int::~Input_int() {
    delete [] ptr;
    delete [] F_upper;
    delete [] F_lower;
}



Input_Reads::Input_Reads(const char *filename): m(M), n(N), r(R), var(VAR), ref(REF),n_cluster(n_c),n_c(0),
ave_ref(AVE_REF),ave_var(AVE_VAR),cl(cluster),cl_cnt(cluster_cnt){
    std::ifstream fin(filename);
    fin >> M >> N;

    VAR = new int *[m];
    REF = new int *[m];
    ptr = new int[m*n*2];
    cluster = new int [n];

    for (int i = 0; i < m; i++) {
        VAR[i] = ptr + i * n;
        REF[i] = ptr + (i + m) * n;
        for (int j = 0; j < n; j++){
            fin >> VAR[i][j] >> REF[i][j];
        }
    }

    for(int i = 0; i < n; i++){
        fin >> cluster[i];
        n_c = std::max(n_c,cluster[i]);
    }
    n_c+=1;
    cluster_cnt=new int[n_c];
    std::fill(cluster_cnt, cluster_cnt + n_c, 0);
    for(int i = 0; i < n; i++) {
        cluster_cnt[cluster[i]]++;
    }

    AVE_VAR = new double *[m];
    AVE_REF = new double *[m];
    ptr_ = new double[m*n_c*2];
    std::fill(ptr_,ptr_+n_c*m*2,0);
    for (int i = 0; i < m; i++) {
        AVE_VAR[i] = ptr_ + i * n_c;
        AVE_REF[i] = ptr_ + (i + m) * n_c;
    }
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++) {
            AVE_REF[i][cluster[j]]+=(VAR[i][j]+REF[i][j]); //depth
            AVE_VAR[i][cluster[j]]+=(double)VAR[i][j]/(VAR[i][j]+REF[i][j]); //vaf
        }
    }
    for(int i = 0; i < m; i++){
        double tmp_depth,tmp_vaf;
        for(int j = 0; j < n_c; j++){
            tmp_depth = AVE_REF[i][j] / cluster_cnt[j];
            tmp_vaf = AVE_VAR[i][j] / cluster_cnt[j];
            AVE_REF[i][j] = tmp_depth*(1-tmp_vaf);
            AVE_VAR[i][j] = tmp_depth*tmp_vaf;
        }
    }


    fin >> R;
}

Input_Reads::~Input_Reads() {
    delete [] ptr;
    delete [] VAR;
    delete [] REF;
    delete [] cluster_cnt;
}


