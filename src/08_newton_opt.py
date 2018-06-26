from __future__ import division, print_function
from numpy import array
from numpy.linalg import norm, solve


def newton_opt(fun, grad, hess, x, niter=50, gtol=1e-8, verbose=False):
    msg = "Maximum number of iterations reached."
    g = grad(x)
    for cont in range(niter):
        if verbose:
            print("n: {}, x: {}, g: {}".format(cont, x, g))
        x = x - solve(hess(x), g)
        g = grad(x)
        if norm(g) < gtol:
            msg = "Extremum found with desired accuracy."
            break
    return x, fun(x), msg


def rosen(x):
    return (1 - x[0])**2 + 100*(x[1] - x[0]**2)**2


def rosen_grad(x):
    return array([
        -2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]**2),
        200*(x[1] - x[0]**2)])

def rosen_hess(x):
    return array([[-400*(x[1]-x[0]**2) + 800*x[0]**2 + 2, -400*x[0]],
                 [-400*x[0], 200]])


print(newton_opt(rosen, rosen_grad, rosen_hess, [2.0, 1.0], verbose=True))
