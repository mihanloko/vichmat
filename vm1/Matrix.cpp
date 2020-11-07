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

vector<float> &Matrix::operator[](int ind) {
    return matr[ind];
}

ostream& operator<<(ostream &out, const Matrix &x) {
    for (const vector<float> &v: x.matr) {
        for (float f : v) {
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
        res.matr.push_back(vector<float>());
        for (int j = 0; j < columnCount; j++) {
            float sum = 0;
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

/**
 * Решение системы уравнений
 * @param b матрица столбец свободных членов
 * @return матрица х решение системы
 */
Matrix Matrix::solve(Matrix b) {
    Matrix result;
    int rowCount = matr.size();
    vector<float> y(rowCount);
    // вычисляем y из уравнения L * y = b
    for (int i = 0; i < rowCount; ++i) {
        float res = b[i][0];
        for (int j = 0; j < i; j++) {
            res -= lu[i][j] * y[j];
        }
        y[i] = res;
    }
    vector<float> x(rowCount);
    // вычисляем х из уравнения U * x = y
    for (int i = rowCount - 1; i > -1; --i) {
        float res = y[i];
        for (int j = i + 1; j < rowCount; j++) {
            res -= lu[i][j] * x[j];
        }
        x[i] = res / lu[i][i];
    }

    result.matr.resize(rowCount);
    for (int i = 0; i < rowCount; ++i) {
        result[i].push_back(x[i]);
    }
    return result;
}

/**
 * вычисление обратной матрицы
 * @return обратная матрица
 */
Matrix Matrix::getInverseMatrix() {
    Matrix inverse;
    long rowCount = matr.size();
    vector<vector<float>> inv(rowCount, vector<float>(rowCount));
    inverse.matr = inv;

    // все проходы с конца в начало
    for (int j = rowCount - 1; j >= 0; --j) {
        // (j;j) сначала вычисляется диагональный элемент
        inverse[j][j] = 1;
        for (int k = j + 1; k < rowCount; ++k) {
            inverse[j][j] -= lu[j][k] * inverse[k][j];
        }
        inverse[j][j] /= lu[j][j];
        //column затем столбец
        for (int i = j - 1; i >= 0; --i) {
            float res = 0;
            for (int k = i + 1; k < rowCount; k++)
                res += lu[i][k] * inverse[k][j];
            inverse[i][j] = -res / lu[i][i];
        }
        //row затем строка
        for (int i = j - 1; i >= 0; --i) {
            float res = 0;
            for (int k = i + 1; k < rowCount; k++)
                res += inverse[j][k] * lu[k][i];
            inverse[j][i] = -res;
        }
    }

    return inverse;
}

istream& operator>>(istream &in, Matrix &x) {
    string str;
    getline(in, str);
    int row = 0;
    while (!str.empty()) {
        x.matr.push_back(vector<float>());
        vector<string> dataS;
        vector<float> dataF;
        split(str, " ", dataS);
        for (string s : dataS) {
            x[row].push_back(std::stof(s));
        }
        row++;
        getline(in, str);
    }
    return in;
}

/**
 * вычисления разложения матрицы А
 * @return удалось разложить матрицу?
 */
bool Matrix::calculateLU() {
    long n = matr.size();
    lu.assign(n, vector<float>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float res = matr[i][j];
            if (i <= j) {
                for (int k = 0; k < i; k++) {
                    res -= lu[i][k] * lu[k][j];
                }
            } else {
                for (int k = 0; k < j; k++) {
                    res -= lu[i][k] * lu[k][j];
                }
                if (abs(lu[j][j]) < 0.0001)
                    return false;
                res /= lu[j][j];
            }
            lu[i][j] = res;
        }
    }
    return true;
}

// вычисление нормы матрицы как максимум суммы модулей строки
float Matrix::calculateNorm() {
    float maxSum = 0;
    long n = matr.size();
    for (int i = 0; i < n; ++i) {
        float res = 0;
        for (int j = 0; j < n; ++j) {
            res += abs(matr[i][j]);
        }
        maxSum = max(maxSum, res);
    }
    return maxSum;
}

