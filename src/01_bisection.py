from __future__ import division, print_function
from numpy import log2, ceil, abs, cos

def bisection(fun, a, b, xtol=1e-6, ftol=1e-12):
    if fun(a) * fun(b) > 0:
        c = None
        msg = "The function should have a sign change in the interval."
    else:
        nmax = int(ceil(log2((b - a)/xtol)))
        for cont in range(nmax):
            c = 0.5*(a + b)
            if abs(fun(c)) < ftol:
                msg = "Root found with desired accuracy."
                break
            elif fun(a) * fun(c) < 0:
                b = c
            elif fun(b) * fun(c) < 0:
                a = c
            msg = "Maximum number of iterations reached."
    return c, msg

def fun(x):
    return cos(x) - x**2

print(bisection(fun, 0.0, 1.0))

