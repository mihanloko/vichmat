#include <iostream>
#include <fstream>
#include "Matrix.h"

using namespace std;

int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");
    Matrix a;
    in >> a;
    float eps;
    in >> eps;
    Matrix E;
    E.matr.resize(a.matr.size(), vector<double>(a.matr.size()));
    for (int i = 0; i < a.matr.size(); i++)
        E[i][i] = 1;

    vector<float> numbers;
    vector<Matrix> vectors;
//    a.getAll(eps, numbers, vectors);

    for (int i = 0; i < numbers.size(); i++) {
        out << "Собственное число " << numbers[i] << endl;
        out << "Собственный вектор\n" << vectors[i];
        out << "Невязка\n" << a * vectors[i] - vectors[i] * numbers[i] << endl << endl;
    }


    return 0;
}