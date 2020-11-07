//
// Created by mikhail on 02.02.19.
//

#ifndef VM1_MATRIX_H
#define VM1_MATRIX_H

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class Matrix {
public:
    vector<vector<float>> matr;
    vector<vector<float>> lu;

    vector<float >& operator [](int ind);
    Matrix operator * (Matrix x);
    Matrix operator - (Matrix x);
    friend ostream& operator << (ostream & out, const Matrix& x);
    friend istream& operator >> (istream &in, Matrix &x);

    Matrix solve(Matrix b);
    Matrix getInverseMatrix();
    bool calculateLU();
    float calculateNorm();
};


#endif //VM1_MATRIX_H
