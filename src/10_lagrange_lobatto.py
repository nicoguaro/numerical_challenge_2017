from __future__ import division
from numpy import (zeros_like, pi, cos, zeros, amax, abs,
                   array, linspace, prod)
import matplotlib.pyplot as plt


def lagrange(x_int, y_int, x_new):
    y_new = zeros_like(x_new)
    for xi, yi in zip(x_int, y_int):
        y_new += yi*prod([(x_new - xj)/(xi - xj)
                         for xj in x_int if xi!=xj], axis=0)
    return y_new


def gauss_lobatto(N, tol=1e-15):
    x = -cos(linspace(0, pi, N))
    P = zeros((N, N))  # Vandermonde Matrix
    x_old = 2
    while amax(abs(x - x_old)) > tol:
        x_old = x
        P[:, 0] = 1
        P[:, 1] = x
        for k in range(2, N):
            P[:, k] = ((2 * k - 1) * x * P[:, k - 1] -
                       (k - 1) * P[:, k - 2]) / k
        x = x_old - (x * P[:, N - 1] - P[:, N - 2]) / (N * P[:, N - 1])
        print(x)
    return array(x)


runge = lambda x: 1/(1 + 25*x**2)
x = linspace(-1, 1, 100)
x_int = linspace(-1, 1, 11)
x_int2 = gauss_lobatto(11)
x_new = linspace(-1, 1, 1000)
y_new = lagrange(x_int, runge(x_int), x_new)
y_new2 = lagrange(x_int2, runge(x_int2), x_new)
plt.plot(x, runge(x), "k")
plt.plot(x_new, y_new)
plt.plot(x_new, y_new2)
plt.legend(["Runge function",
            "Uniform interpolation",
            "Lobatto-sampling interpolation"])
plt.xlabel("x")
plt.ylabel("y")
plt.show()
