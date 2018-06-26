from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def vander_mat(x):
    n = len(x)
    van = np.zeros((n, n))
    power = np.array(range(n))
    for row in range(n):
        van[row, :] = x[row]**power
    return van


def inter_coef(x):
    vand_mat = vander_mat(x)
    coef = np.linalg.solve(vand_mat, np.eye(len(x)))
    return coef


def compute_interp(x, f, x_eval):
    n = len(x)
    coef = inter_coef(x)
    f_eval = np.zeros_like(x_eval)
    for row in range(n):
        for col in range(n):
            f_eval += coef[row, col]*x_eval**row*f[col]
    return f_eval


n = 11
x = -np.cos(np.linspace(0, np.pi, n))
f = 1/(1 + 25*x**2)
x_eval = np.linspace(-1, 1, 500)
interp_f = compute_interp(x, f, x_eval)
plt.figure()
plt.plot(x_eval, 1/(1 + 25*x_eval**2))
plt.plot(x_eval, interp_f)
plt.plot(x, f, ".")
plt.ylim(0, 1.2)
plt.show()