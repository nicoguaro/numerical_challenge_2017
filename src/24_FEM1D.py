from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def FEM1D(coords, source):
    N = len(coords)
    stiff_loc = np.array([[2.0, -2.0], [-2.0, 2.0]])
    eles = [np.array([cont, cont + 1]) for cont in range(0, N - 1)]
    stiff = np.zeros((N, N))
    rhs = np.zeros(N)
    for ele in eles:
        jaco = coords[ele[1]] - coords[ele[0]]
        rhs[ele] = rhs[ele] + jaco*source(coords[ele])
        for cont1, row in enumerate(ele):
            for cont2, col in enumerate(ele):
                stiff[row, col] = stiff[row, col] +  stiff_loc[cont1, cont2]/jaco
    return stiff, rhs


N = 100
fun = lambda x: x**3
x = np.linspace(0, 1, N)
stiff, rhs = FEM1D(x, fun)
sol = np.zeros(N)
sol[1:-1] = np.linalg.solve(stiff[1:-1, 1:-1], -rhs[1:-1])


#%% Plotting
plt.figure(figsize=(4, 3))
plt.plot(x, sol)
plt.plot(x, x*(x**4 - 1)/20, linestyle="dashed")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.legend(["FEM solution", "Exact solution"])
plt.tight_layout()
plt.show()


