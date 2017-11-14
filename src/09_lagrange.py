from __future__ import division
from numpy import zeros_like, linspace, exp, prod
import matplotlib.pyplot as plt


def lagrange(x_int, y_int, x_new):
    y_new = zeros_like(x_new)
    for xi, yi in zip(x_int, y_int):
        y_new += yi*prod([(x_new - xj)/(xi - xj)
                         for xj in x_int if xi!=xj], axis=0)
    return y_new


x_int = linspace(-5, 5, 11)
y_int = 1/(1 + exp(-x_int))
x_new = linspace(-5, 5, 1000)
y_new = lagrange(x_int, y_int, x_new)
plt.plot(x_int, y_int, "ok")
plt.plot(x_new, y_new)
plt.show()
