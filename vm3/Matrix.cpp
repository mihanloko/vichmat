#include "Matrix.h"

void split(const string &str, const string &delim, vector<string> &parts) {
    size_t start, end = 0;
    while (end < str.size()) {
        start = end;
        while (start < str.size() && (delim.find(str[start]) != string::npos)) {
            start++;  // skip initial whitespace
        }
        end = start;
        while (end < str.size() && (delim.find(str[end]) == string::npos)) {
            end++; // skip to end of word
        }
        if (end - start != 0) {  // just ignore zero-length strings.
            parts.emplace_back(str, start, end - start);
        }
    }
}

vector<double> &Matrix::operator[](int ind) {
    return matr[ind];
}

ostream& operator<<(ostream &out, const Matrix &x) {
    for (const vector<double> &v: x.matr) {
        for (double f : v) {
            out.width(15);
            //out.setf(ios::left);
            out << f;
        }
        out << endl;
    }
    return out;
}

Matrix Matrix::operator*(Matrix x) {
    Matrix res;
    if (matr[0].size() != x.matr.size())
        return res;

    int rowCount = matr.size();
    int columnCount = x[0].size();
    int count = matr[0].size();
    for (int i = 0; i < rowCount; i++) {
        res.matr.push_back(vector<double>());
        for (int j = 0; j < columnCount; j++) {
            double sum = 0;
            for (int k = 0; k < count; k++) {
                sum += matr[i][k] * x[k][j];
            }
            res[i].push_back(sum);
        }
    }
    return res;
}

Matrix Matrix::operator-(Matrix x) {
    Matrix result = *this;
    for (int i = 0; i < matr.size(); i++)
        for (int j = 0; j < matr[0].size(); j++)
            result[i][j] -= x[i][j];
    return result;
}

istream& operator>>(istream &in, Matrix &x) {
    string str;
    getline(in, str);
    int row = 0;
    while (!str.empty()) {
        x.matr.push_back(vector<double>());
        vector<string> dataS;
        vector<double> dataF;
        split(str, " ", dataS);
        for (string s : dataS) {
            x[row].push_back(std::stof(s));
        }
        row++;
        getline(in, str);
    }
    return in;
}


// вычисление нормы матрицы как максимум суммы модулей строки
double Matrix::calculateNorm() {
    double maxSum = 0;
    long n = matr.size();
    for (int i = 0; i < n; ++i) {
        double res = 0;
        for (int j = 0; j < matr[i].size(); ++j) {
            res += abs(matr[i][j]);
        }
        maxSum = max(maxSum, res);
    }
    return maxSum;
}

vector<double> Matrix::operator*(vector<double> x) {

    unsigned long size = matr.size();
    vector<double> delta(size);

    for (int i = 0; i < size; i++) {
        delta[i] = 0;
        for (int j = 0; j < size; j++) {
            delta[i] += matr[i][j] * x[j];
        }
    }
    return delta;
}

double inline Matrix::scalar(vector<double> &a, vector<double> &b) {
    double  sum = 0;
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

Matrix Matrix::transpon() {
    Matrix t;
    long size = matr.size();
    t.matr.assign(size, vector<double>(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            t.matr[i][j] = matr[j][i];
        }
    }
    return t;
}

Matrix Matrix::getHouseholder() {
    Matrix E;
    E.matr.assign(matr.size(), vector<double>(matr.size(), 0));
    for (int i = 0; i < matr.size(); i++)
        E[i][i] = 1;

    int n = matr.size();
    Matrix tempA = *this;
    Matrix resU = E;
    Matrix w;
//    w.matr.resize(n, vector<double>(1, 0));
    for (int i = 0; i < n - 2; i++) {
        w.matr.resize(n, vector<double>(1, 0));
        int rowI = i + 1;
        float s = sign(tempA[rowI][i]) * tempA.col(i).tail(n - rowI).norm();
        float alpha = 1 / sqrt(2 * s * (s - tempA(rowI, i)));
        w(rowI) = alpha * (tempA(rowI, i) - s);
        for (int j = rowI + 1; j < n; ++j)
            w(j) = alpha * tempA(j, i);
        auto U = E - 2 * w * w.transpose();
        resU = resU * U;
        tempA = U * tempA * U;
    }
    return resU;


    return Matrix();
}

Matrix Matrix::operator/(double n) {
    Matrix res = *this;
    for (int i = 0; i < matr.size(); i++)
        for (int j = 0; j < matr[i].size(); j++)
            res[i][j] /= n;
    return res;
}

Matrix Matrix::operator*(double x) {
    Matrix res = *this;
    for (int i = 0; i < matr.size(); i++)
        for (int j = 0; j < matr[i].size(); j++)
            res[i][j] *= x;
    return res;
}