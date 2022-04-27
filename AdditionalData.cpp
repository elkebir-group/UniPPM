//
// Created by Yuanyuan Qi on 4/4/22.
//

#include "AdditionalData.h"

#include <fstream>

AdditionalData::AdditionalData(const int &n, const std::string &file_name): n(n) {
    if (file_name.empty()) {
        n_c = 0;
        ptr = nullptr;
        data = nullptr;
    } else {
        std::ifstream fin(file_name);
        fin >> n_c;
        ptr = new short[n_c * n];
        data = new short *[n_c];
        for (int i = 0; i < n_c; i++) {
            data[i] = ptr + i * n;
            for (int j = 0; j < n; j++)
                fin >> data[i][j];
        }
    }
    std::cout << "[UniPPM] " << n_c << " additional lines of data." << std::endl;
}

AdditionalData::~AdditionalData() {
    delete [] ptr;
    delete [] data;
}

void AdditionalData::Enforce_constraints(CNF * F, const std::vector<std::vector<CMSat::Lit> > & relation) const{
    for (int i = 0; i < n_c; i++) {
        for (int p = 0; p < n; p++) {
            if (!data[i][p]) continue;
            for (int q = p+1; q < n; q++){
                if(!data[i][q]) continue;
                switch (data[i][p]) {
                    case 1:
                        switch (data[i][q]) {
                            case 1:
                                F->add_clause({relation[p][q],relation[q][p]});
                                break;
                            case -1:
                                F->add_clause({~relation[q][p]});
                                break;
                        }
                        break;
                    case -1:
                        switch (data[i][q]) {
                            case 1:
                                F->add_clause({~relation[p][q]});
                                break;
                            case -1:
                                break;
                        }
                        break;
                }
            }
        }
    }
}
