//
// Created by mikhail on 02.02.19.
//

#ifndef VM1_MATRIX_H
#define VM1_MATRIX_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#define initial 1

using namespace std;

class Matrix {
public:
    vector<vector<double>> matr;
    vector<vector<double>> lu;

    vector<double >& operator [](int ind);
    Matrix operator * (Matrix x);
    vector<double> operator * (vector<double> x);
    Matrix operator - (Matrix x);
    friend ostream& operator << (ostream & out, const Matrix& x);
    friend istream& operator >> (istream &in, Matrix &x);
    double inline scalar(vector<double> &a, vector<double> &b);
    Matrix transpon();

    double calculateNorm();

    Matrix solveSeidel(Matrix b, double eps, int &num, double &norm);
    Matrix minimalResiduals(Matrix b, double eps, int &num, double &norm);
};


#endif //VM1_MATRIX_H
