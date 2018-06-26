from __future__ import division, print_function
from numpy import abs, cos, sin

def newton(fun, grad, x, niter=50, ftol=1e-12, verbose=False):
    msg = "Maximum number of iterations reached."
    for cont in range(niter):
        if abs(grad(x)) < ftol:
            x = None
            msg = "Derivative near to zero."
            break
        if verbose:
            print("n: {}, x: {}".format(cont, x))
        x = x - fun(x)/grad(x)
        if abs(fun(x)) < ftol:
            msg = "Root found with desired accuracy."
            break
    return x, msg


def fun(x):
    return cos(x) - x


def grad(x):
    return -sin(x) - 1.0


print(newton(fun, grad, 1.0))
