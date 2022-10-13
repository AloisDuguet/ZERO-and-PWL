#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;

arma::sp_mat readMatrix_CSV(string filename, int size_row, int size_col) {
    vector<long long unsigned int> location_u;
    vector<long long unsigned int> location_m;
    vector<double> values;

    arma::sp_mat M(size_row, size_col);

    ifstream file(filename);
    int a, b;
    double c;
    while(file >> a >> b >> c) {
        M(a,b) = c;
        location_u.push_back(a);
        location_m.push_back(b);
        values.push_back(c);
    }

    umat lu(location_u);
    umat lm(location_m);
    umat location(join_rows(lu, lm).t());

    //return arma::sp_mat(location, vec(values));
    return M;
}

arma::vec readVector_CSV(string filename, int size_row) {
    vector<double> values;

    arma::vec v(size_row);

    ifstream file(filename);
    double val;
    for(int i = 0; file >> val; i++)
       v.at(i) = val;

    return v;
}

arma::vec readSparseVector_CSV(string filename, int size_row) {
    // parse linear terms c which are written in "row value" format
    arma::vec c(size_row);

    ifstream file(filename);
    int a;
    double val;
    while(file >> a >> val)
        c.at(a) = val;

    return c;
}

arma::vec readIntegerIndexes_CSV(string filename) {
    vector<int> values;

    ifstream file(filename);
    int v;
    while(file >> v)
       values.push_back(v);

    return conv_to<arma::vec>::from(values);
}



// Some of your code here
// MatrixA.csv should have on each line "row column value"
//A = readMatrix_CSV("MatrixA.csv")
// Vectorb.csv should have on each line "value"
//b = readVector_CSV("Vectorb.csv")


// Some of your code here
