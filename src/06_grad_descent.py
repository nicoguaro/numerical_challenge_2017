from __future__ import division, print_function
from numpy import array
from numpy.linalg import norm


def grad_descent(fun, grad, x, niter=50, gtol=1e-8, verbose=False):
    msg = "Maximum number of iterations reached."
    g_old = grad(x)
    gamma = 0.1
    for cont in range(niter):
        if verbose:
            print("n: {}, x: {}, g: {}".format(cont, x, g_old))
        dx = -gamma*g_old
        x = x + dx
        g = grad(x)
        dg = g - g_old
        g_old = g
        gamma = dx.dot(dg)/dg.dot(dg)
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

    
print(grad_descent(rosen, rosen_grad, [2.0, 1.0], verbose=True))
