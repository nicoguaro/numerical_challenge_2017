from __future__ import division, print_function
from numpy import array, outer, dot
from numpy.linalg import solve, norm, det


def broyden(fun, jaco, x, niter=50, ftol=1e-12, verbose=False):
    msg = "Maximum number of iterations reached."
    J = jaco(x)
    for cont in range(niter):
        if det(J) < ftol:
            x = None
            msg = "Derivative near to zero."
            break
        if verbose:
            print("n: {}, x: {}".format(cont, x))
        f_old = fun(x)
        dx = -solve(J, f_old)
        x = x + dx
        f = fun(x)
        df = f - f_old
        J = J + outer(df - dot(J, dx), dx)/dot(dx, dx)
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


print(broyden(fun, jaco, [1.0, 2.0]))
