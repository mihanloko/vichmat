#include <iostream>
#include <fstream>
#include <ctime>
#include "Matrix.h"

using namespace std;

int main() {
    srand(time(NULL));
    ifstream in("input.txt");
    ofstream out("output.txt");
    Matrix a;
    in >> a;
    float eps;
    in >> eps;

    vector<float> numbers;
    vector<Matrix> vectors;
    Matrix U = a.getHouseholder();
    Matrix B = U.transpon() * a * U;
    Matrix tempA = B;
    Matrix Q, R;
    int iterCount = 0;
    while (!tempA.checkEnd(eps)) {
        tempA.QR(Q, R);
        tempA = R * Q;
        iterCount++;
    }

    out << "Исходная матрица:\n" << a << endl;
    Matrix E;
    E.matr.resize(a.matr.size(), vector<double>(a.matr.size()));
    for (int i = 0; i < a.matr.size(); i++)
        E[i][i] = 1;
    for (int i = 0; i < a.matr.size(); i++) {
        out << "Собственное число " << tempA[i][i] << endl;
        Matrix ownVector = (B - E * tempA[i][i]).getX(eps);
        ownVector = U * ownVector;
        out << "Собственный вектор:\n" << ownVector;
        out << "Невязка:\n" << (a * ownVector) - (ownVector * tempA[i][i]);
    }

    return 0;
}