from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def jacobi(u, update, tol=1e-6, niter=500):
    for n in range(niter):
        u_new = update(u)
        if np.linalg.norm(u_new - u) < tol:
            break
        else:
            u = u_new.copy()
    return u, n


def heat_FDM(u):
    u_new = u.copy()
    u_new[1:-1, 1:-1] = 0.25*(u_new[0:-2, 1:-1] + u_new[2:, 1:-1] +\
         u_new[1:-1, 0:-2] + u_new[1:-1, 2:])
    return u_new

    
nx = 50
ny = 50
x_vec = np.linspace(-0.5, 0.5, nx)
y_vec = np.linspace(-0.5, 0.5, ny)
u0 = np.zeros((nx, ny))
u0[:, 0] = 1 - x_vec
u0[0, :] = 1 - y_vec
nvec = [100, 1000, 10000, 100000]
for num, niter in enumerate(nvec):
    u, n = jacobi(u0, heat_FDM, tol=1e-12, niter=niter)
    plt.subplot(2, 2, num + 1)
    plt.contourf(u, cmap='hot')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.title('{} iterations'.format(n + 1))
    plt.axis('image')
plt.tight_layout()
plt.show()