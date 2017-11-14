from __future__ import division, print_function
import numpy as np


def cholesky(mat):
    m, _ = mat.shape
    mat = mat.copy()
    for col in range(m):
        mat[col, col] = np.sqrt(mat[col, col] -
               mat[col, 0:col].dot(mat[col, 0:col]))
        for row in range(col + 1, m):
            mat[row, col] = (mat[row, col] -
               mat[row, 0:col].dot(mat[col, 0:col]))/mat[col, col]
    for row in range(1, m):
        mat[0:row, row] = 0
    return mat


A = np.array([
    [4, -1, 1],
    [-1, 4.25, 2.75],
    [1, 2.75, 3.5]], dtype=float)
B = cholesky(A)