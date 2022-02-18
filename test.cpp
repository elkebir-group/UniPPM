#include "Input.h"
#include "AncestryGraph.h"
#include "Solver.h"

using namespace std;
int main () {
    Input_Reads raw_in("test.txt");
    Input transform_in(raw_in,0.1);
    Input_int In(transform_in, 4);
    cout<<"input created"<<endl;
    AncestryGraph Gf(In);
    cout<<"GF created"<<endl;
    Solver solver(Gf,0);
    cout<<"class created"<<endl;

    solver.self_solver().to_file("test.cnf");
    vector<int> res;
    auto b = solver.attempt(solver.self_solver(),&res);
    cout<<res.size()<<endl;
    for(auto i:res){
        cout<<i<<' ';
    }
    cout<<endl;
    auto a = solver.sampling(solver.self_solver(),4);
    for (auto vec:a) {
        for (auto i:vec)
            cout<<i<<' ';
        cout<<endl;
    }
    return 0;
}