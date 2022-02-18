#include <iostream>

#include "Input.h"

using namespace std;

int main() {
    int t = 2;
    Input ifs = ("test.txt");
    Input_int ifi(ifs, 10);
    for (int i = 0; i < ifi.m; i++){
        for (int j = 0 ; j < ifi.n; j++)
            cout<<ifi.F_l[i][j]<<' ';
        cout<<endl;
    }
    return 0;
}
