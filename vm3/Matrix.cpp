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
    long h = matr.size(), w = matr[0].size();
    t.matr.assign(w, vector<double>(h));
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
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
    for (int i = 0; i < n - 2; i++) {
        w.matr.assign(n, vector<double>(1, 0));
        int rowI = i + 1;
        double s = sign(tempA[rowI][i]) * norm(tempA.getTailColumn(i, rowI));
        double alpha = 1.0 / sqrt(2 * s * (s - tempA[rowI][i]));
        w[rowI][0] = alpha * (tempA[rowI][i] - s);
        for (int j = rowI + 1; j < n; ++j)
            w[j][0] = alpha * tempA[j][i];
        auto U = E - w * w.transpon() * 2;
        resU = resU * U;
        tempA = U * tempA * U;
    }
    /*for (int i = 0; i < n - 2; i++) {
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
    }*/
    return resU;
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

vector<double> Matrix::getTailColumn(int i, int s) {
    vector<double> result;
    for (int j = s; j < matr.size(); j++)
        result.push_back(matr[j][i]);
    return result;
}

double Matrix::norm(const vector<double>& column) {
    double sum = 0;
    for (double i: column)
        sum += i * i;
    return sqrt(sum);
}

void Matrix::QR(Matrix &Q, Matrix &R) {
    int n = matr.size();
    Matrix E;
    E.matr.resize(matr.size(), vector<double>(matr.size()));
    for (int i = 0; i < matr.size(); i++)
        E[i][i] = 1;
    
    R = *this;
    Q = E;
    for (int i = 0; i < n - 1; ++i) {
        int nextI = i + 1;
        double tau = atan(R[i][i] / R[nextI][i]);
        double a = sin(tau), b = cos(tau);
        for (int j = 0; j < n; j++) {
            auto firstR = R[i][j] * a + R[nextI][j] * b;
            auto secondR = -1 * R[i][j] * b + R[nextI][j] * a;
            auto firstQ = Q[j][i] * a + Q[j][nextI] * b;
            auto secondQ = -1 * Q[j][i] * b + Q[j][nextI] * a;
            R[i][j] = firstR;
            R[nextI][j] = secondR;
            Q[j][i] = firstQ;
            Q[j][nextI] = secondQ;
        }
    }
}

bool Matrix::checkEnd() {
    auto n = matr.size();
    for (int j = 0; j < n - 1; ++j)
        if (fabs(matr[j + 1][j]) >= 0.0001)
            return false;
    return true;
}
