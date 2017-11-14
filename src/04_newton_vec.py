from __future__ import division, print_function
from numpy import array
from numpy.linalg import solve, norm, det


def newton(fun, jaco, x, niter=50, ftol=1e-12, verbose=False):
    msg = "Maximum number of iterations reached."
    for cont in range(niter):
        J = jaco(x)
        f = fun(x)
        if det(J) < ftol:
            x = None
            msg = "Derivative near to zero."
            break
        if verbose:
            print("n: {}, x: {}".format(cont, x))
        x = x - solve(J, f)
        if norm(f) < ftol:
            msg = "Root found with desired accuracy."
            break
    return x, msg


def fun(x):
    return array([x[0] + 2*x[1] - 2, x[0]**2 + 4*x[1]**2 - 4])


def jaco(x):
    return array([
            [1, 2],
            [2*x[0], 8*x[1]]])


print(newton(fun, jaco, array([1.0, 10.0])))
