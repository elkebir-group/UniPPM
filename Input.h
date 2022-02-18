//
// Created by Yuanyuan Qi on 2/11/22.
//

#ifndef UNIPPM_INPUT_H
#define UNIPPM_INPUT_H



class Input_Reads {
public:
    Input_Reads(const char * filename);

    ~Input_Reads();

    const int &m,&n,&r;

    const int *const *const &var, *const *const &ref;
private:
    int M,N;
    int **VAR, **REF, *ptr;
    int R;
};


class Input {
public:
    Input(const char * filename);

    Input(const Input_Reads& In, const double &alpha);

    ~Input();

    const int &m,&n,&r;

    const double *const *const &F_l, *const *const &F_u;

private:
    int M,N;
    double **F_lower, **F_upper, *ptr;
    int R;

};


class Input_int {
public:
    Input_int(const Input &, const int & N_bits);

    ~Input_int();

    const int &m, &n, &r, &n_bits;

    const int *const *const &F_l, *const *const &F_u;
private:
    int **F_lower, **F_upper, *ptr;
    int Nn;
};


#endif //UNIPPM_INPUT_H
