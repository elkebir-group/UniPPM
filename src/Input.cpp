//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "Input.h"

#include <fstream>

Input::Input(const std::string &filename, int r): r(r) {
    std::ifstream fin(filename);
    fin >> m >> n;
    data = std::vector<std::vector<std::pair<double, double> > >(m,std::vector<std::pair<double, double> >(n));
    for (auto i = 0; i < m; ++i) {
        for (auto j = 0; j < n; j++) {
            fin >> data[i][j].first >> data[i][j].second;
        }
    }
}
