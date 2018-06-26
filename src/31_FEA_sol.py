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


def conj_grad(A, b, x, tol=1e-8):
    r = b - A.dot(x)
    p = r
    rsq_old = r.dot(r)
    for cont in range(len(b)):
        Ap = A.dot(p)
        alpha = rsq_old / p.dot(Ap)
        x = x + alpha*p
        r = r - alpha*Ap
        rsq_new = r.dot(r)
        if np.sqrt(rsq_new) < tol:
            break
        p = r + (rsq_new / rsq_old) * p
        rsq_old = rsq_new
    return x, cont, np.sqrt(rsq_new)


fun = lambda x: x**3
N_vec = np.logspace(0.5, 3, 50, dtype=int)
err_vec = np.zeros_like(N_vec, dtype=float)
niter_vec = np.zeros_like(N_vec)
plt.figure(figsize=(8, 3))
for cont, N in enumerate(N_vec):
    x = np.linspace(0, 1, N)
    stiff, rhs = FEM1D(x, fun)
    sol = np.zeros(N)
    sol[1:-1], niter, _ = conj_grad(stiff[1:-1, 1:-1], -rhs[1:-1], rhs[1:-1])
    err = np.linalg.norm(sol - x*(x**4 - 1)/20)/np.linalg.norm(x*(x**4 - 1)/20)
    err_vec[cont] = err
    niter_vec[cont] = niter

plt.subplot(121)
plt.loglog(N_vec, err_vec)
plt.xlabel("Number of nodes")
plt.ylabel("Relative error")
plt.subplot(122)
plt.loglog(N_vec, niter_vec)
plt.xlabel("Number of nodes")
plt.ylabel("Number of iterations")
plt.tight_layout()
plt.show()


