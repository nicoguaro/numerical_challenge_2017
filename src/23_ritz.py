from __future__ import division, print_function
import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve
import matplotlib.pyplot as plt


def ritz(N, source):
    stiff_mat = np.zeros((N, N))
    rhs = np.zeros((N))
    for row in range(N):
        for col in range(N):
            numer = (2 + 2*row + 2*col + 2*row*col)
            denom = (row + col + 1) * (row + col + 2) * (row + col + 3)
            stiff_mat[row, col] = numer/denom
        fun = lambda x: x**(row + 1)*(1 - x)*source(x)
        rhs[row], _ = quad(fun, 0, 1)
    return stiff_mat, rhs


N = 3
source = lambda x: x**3
mat, rhs = ritz(N, source)
c = solve(mat, -rhs)
x = np.linspace(0, 1, 100)
y = np.zeros_like(x)
for cont in range(N):
    y += c[cont]*x**(cont + 1)*(1 - x)



#%% Plotting
plt.figure(figsize=(4, 3))
plt.plot(x, y)
plt.plot(x, x*(x**4 - 1)/20, linestyle="dashed")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.legend(["Ritz solution", "Exact solution"])
plt.tight_layout()
plt.savefig("ritz-N-3.svg")
plt.show()

