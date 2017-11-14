from __future__ import division, print_function
import numpy as np


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


N = 1000
A = -np.diag(2*np.ones(N)) + np.diag(np.ones(N-1), -1) +\
    np.diag(np.ones(N-1), 1)
b = np.ones(N)
x0 = np.ones(N)
x, niter, accu = conj_grad(A, b, x0)