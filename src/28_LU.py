from __future__ import division, print_function
import numpy as np


def LU(mat):
    m, _ = mat.shape
    mat = mat.copy()
    for col in range(0, m - 1):
        for row in range(col + 1, m):
            if mat[row, col] != 0.0:
                lam = mat[row, col]/mat[col, col]
                mat[row, col + 1:m] = mat[row, col + 1:m] -\
                                      lam * mat[col, col + 1:m]
                mat[row, col] = lam
    return mat


A = np.array([
    [1, 1, 0, 3],
    [2, 1, -1, 1],
    [3, -1, -1, 2],
    [-1, 2, 3, -1]], dtype=float)
B = LU(A)