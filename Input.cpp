//
// Created by Yuanyuan Qi on 2/11/22.
//

#include "Input.h"

#include <boost/math/special_functions/beta.hpp>

#include <cmath>
#include <fstream>



Input::Input(const char * filename) : m(M), n(N), r(R), F_l(F_lower), F_u(F_upper) {
    std::ifstream fin(filename);
    fin >> M >> N;

    F_lower = new double *[m];
    F_upper = new double *[m];
    ptr = new double[m*n*2];

    for (int i = 0; i < m; i++) {
        F_lower[i] = ptr + i * n;
        F_upper[i] = ptr + (i + m) * n;
        for (int j = 0; j < n; j++){
            fin >> F_lower[i][j] >> F_upper[i][j];
        }
    }

    fin >> R;
}

Input::~Input() {
    delete [] ptr;
    delete [] F_upper;
    delete [] F_lower;
}

Input::Input(const Input_Reads &In, const double &alpha): m(M), n(N), r(R), F_l(F_lower), F_u(F_upper) {
    M = In.m;
    N = In.n;

    F_lower = new double *[m];
    F_upper = new double *[m];
    ptr = new double[m*n*2];

    for (int i = 0; i < m; i++) {
        F_lower[i] = ptr + i * n;
        F_upper[i] = ptr + (i + m) * n;
        for (int j = 0; j < n; j++){
//            beta_distribution<> dis(In.var[i][j]+0.5, In.ref[i][j]+0.5)
            F_lower[i][j] = In.var[i][j] ? boost::math::ibeta_inv(In.var[i][j]+0.5, In.ref[i][j]+0.5,(1-alpha)/2) : 0;
            F_upper[i][j] = In.ref[i][j] ? boost::math::ibeta_inv(In.var[i][j]+0.5, In.ref[i][j]+0.5,(1+alpha)/2) : 1;
        }
    }

    R = In.r;
}



Input_int::Input_int(const Input & In, const int &N_bits):m(In.m),n(In.n),r(In.r),Nn(N_bits),n_bits(Nn),F_l(F_lower),F_u(F_upper) {
        F_lower = new int *[m];
        F_upper = new int *[m];
        ptr = new int[m*n*2];

        for (int i = 0; i < m; i++) {
            F_lower[i] = ptr + i * n;
            F_upper[i] = ptr + (i + m) * n;
            for (int j = 0; j < n; j++){
                F_lower[i][j] = int(In.F_l[i][j]*(1<<n_bits));
                F_upper[i][j] = int(In.F_u[i][j]*(1<<n_bits));
            }
        }
}

Input_int::~Input_int() {
    delete [] ptr;
    delete [] F_upper;
    delete [] F_lower;
}



Input_Reads::Input_Reads(const char *filename): m(M), n(N), r(R), var(VAR), ref(REF){
    std::ifstream fin(filename);
    fin >> M >> N;

    VAR = new int *[m];
    REF = new int *[m];
    ptr = new int[m*n*2];

    for (int i = 0; i < m; i++) {
        VAR[i] = ptr + i * n;
        REF[i] = ptr + (i + m) * n;
        for (int j = 0; j < n; j++){
            fin >> VAR[i][j] >> REF[i][j];
        }
    }

    fin >> R;
}

Input_Reads::~Input_Reads() {
    delete [] ptr;
    delete [] VAR;
    delete [] REF;
}


