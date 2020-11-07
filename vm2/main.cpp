#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <chrono>
#include "Matrix.h"

using namespace std;

Matrix generateMatrix(long size) {
    Matrix triangle;
    triangle.matr.assign(size, vector<double>(size, 0));
    for (int i = 0; i < size; i++) {
        for (int j = i; j < size; j++) {
            triangle.matr[i][j] = 10 - rand() % 21;
        }
        if (triangle.matr[i][i] == 0) triangle.matr[i][i]++;
    }
    return triangle * triangle.transpon();
}

Matrix generateVector(long size) {
    Matrix mVector;
    mVector.matr.assign(size, vector<double>(1, 0));
    for (int i = 0; i < size; i++) {
        mVector.matr[i][0] = 10 - rand() % 21;
    }
    return mVector;
}

int main() {
    srand(time(NULL));
    ofstream out("output.txt");

    long size = 10;
    double eps = 0.0000001;

    Matrix a = generateMatrix(size);
    Matrix x = generateVector(size);
    Matrix b = a * x;
//    out << "x\n" << x << endl;
    out << "Размерность матрицы = " << size << "\nТочность = " << eps << endl;

    int num = 0;
    double norm;

    auto begin = chrono::high_resolution_clock::now();
    Matrix seidel = a.solveSeidel(b, eps, num, norm);
    auto elapsedTime = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count();

    out << "Метод Зейделя\nНорма невзязки = " << (b - a * seidel).calculateNorm() << endl;
    out << "Число итераций = " << num << endl;
    out << "Норма delta = " << norm << endl;
    out << "Время работы = " << elapsedTime << endl << endl;

    num = 0;
    begin = chrono::high_resolution_clock::now();
    Matrix minimalResiduals = a.minimalResiduals(b, eps, num, norm);
    elapsedTime = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count();

    out << "Метод минимальных невязок\nНорма невзязки = " << (b - a * minimalResiduals).calculateNorm() << endl;
    out << "Число итераций = " << num << endl;
    out << "Норма delta = " << norm << endl;
    out << "Время работы = " << elapsedTime << endl;

//    out << x << endl << seidel << endl << minimalResiduals;

    return 0;
}