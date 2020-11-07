import numpy as np
import math


def input_matrix(file_name):
    with open(file_name) as f:
        matrix = np.array([list(map(float, row.split())) for row in f.readlines()])
    return matrix


def sign(num):
    return 1 if num > 0 else -1


def norm(vector):
    s = 0
    for x in vector:
        s += x * x
    return math.sqrt(s)


def householder(a):
    n = a.shape[0]
    temp_aa = a.copy()
    res_u = np.eye(n)
    for i in range(n - 2):
        w = np.zeros((4, 1))
        row_i = i + 1
        s = sign(temp_aa[row_i, i]) * norm(temp_aa[row_i:, i])
        alpha = 1. / math.sqrt(2 * s * (s - temp_aa[row_i, i]))
        w[row_i, 0] = alpha * (temp_aa[row_i, i] - s)
        for j in range(row_i + 1, n):
            w[j, 0] = alpha * temp_aa[j, i]
        u = E - 2 * w * w.transpose()
        res_u = res_u * u
        temp_aa = u * temp_aa * u
    return res_u


def check_end(matrix):
    n = matrix.shape[0]
    for j in range(n - 1):
        if math.fabs(matrix[j + 1, j]) >= 0.0001:
            return False
    return True


def QR(B, Q, R):
    n = B.shape[0]
    R = B
    Q = np.eye(n)
    for i in range(n - 1):
        next_i = i + 1
        tau = math.atan(R[i, i] / R[next_i, i])
        a = math.sin(tau)
        b = math.cos(tau)
        for j in range(n):
            first_r = R[i, j] * a + R[next_i, j] * b
            second_r = -1 * R[i, j] * b + R[next_i, j] * a
            first_q = Q[j, i] * a + Q[j, next_i] * b
            second_q = -1 * Q[j, i] * b + Q[j, next_i] * a
            R[i, j] = first_r
            R[next_i, j] = second_r
            Q[j, i] = first_q
            Q[j, next_i] = second_q
    return Q, R


if __name__ == '__main__':
    A = input_matrix("input.txt")
    # print(A)
    E = np.eye(A.shape[0])
    # print(E)
    U = householder(A)
    # print(U)
    B = U.transpose() * A * U
    print(U)
    print(B)
    temp_a = B
    iter_count = 0
    Q = None
    R = None
    while not check_end(temp_a):
        Q, R = QR(temp_a, Q, R)
        temp_a = R * Q
        iter_count += 1
    # print(iter_count)
    # print(temp_a)
