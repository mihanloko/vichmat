#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "Matrix.h"

using namespace std;

int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");
    Matrix a, b;
    in >> a;
    in >> b;

    Matrix oldA = a;

    // вычисляем разложение
    bool isOk = a.calculateLU();
    if (!isOk) {
        out << "Было деление на 0 при вычислении матрицы L";
        return 0;
    }
    Matrix lu;
    lu.matr = a.lu;
    float det = 1;
    for (int i = 0; i < lu.matr.size(); i++)
        det *= lu.matr[i][i];
    out << "Определитель равен " << det << endl;
    out << "L - E + U" << endl << lu << endl;
    if (abs(det) < 0.00001) {
        cout << "Матрица вырожденная" << endl;
        return 0;
    }

    Matrix result = a.solve(b);
    out << "Решение системы" << endl << result << endl;
    out << "Проверка" << endl << a  * result - b << endl;

    Matrix inverse = a.getInverseMatrix();
    out << "Обратная матрица" << endl << inverse << endl;
    out << "Проверка" << endl << inverse * a << endl;

    float normA = a.calculateNorm(), normA1 = inverse.calculateNorm();
    out << "||A|| * ||A^-1|| = " << normA << " * " << normA1 << " = " << normA * normA1;

    return 0;
}