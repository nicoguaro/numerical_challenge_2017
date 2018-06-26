from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def spline_coeff(x, y):
    h = x[1:] - x[:-1]
    d_up = h.copy()
    d_up[0] = 0
    d_down = h.copy()
    d_down[-1] = 0
    d_cent = np.ones_like(x)
    d_cent[1:-1] = 2*(h[:-1] + h[1:])
    mat = np.diag(d_cent) + np.diag(d_up, 1) + np.diag(d_down, -1)
    alpha = np.zeros_like(x)
    alpha[1:-1] = 3/h[1:]*(y[2:] - y[1:-1]) - 3/h[:-1]*(y[1:-1] - y[:-2])
    c = np.linalg.solve(mat, alpha)
    b = np.zeros_like(x)
    d = np.zeros_like(x)
    b[:-1] = (y[1:] - y[:-1])/h - h/3*(c[1:] + 2*c[:-1])
    d[:-1] = (c[1:] - c[:-1])/(3*h)
    return y, b, c, d


n = 20
x = np.linspace(0, 2*np.pi, n)
y = np.sin(x)
a, b, c, d = spline_coeff(x, y)
for cont in range(n - 1):
    x_loc = np.linspace(x[cont], x[cont + 1], 100)
    y_loc = a[cont] + b[cont]*(x_loc - x[cont]) +\
            c[cont]*(x_loc - x[cont])**2 +\
            d[cont]*(x_loc - x[cont])**3
    plt.plot(x_loc, y_loc, "red")
plt.plot(x, y, "o")
plt.show()
