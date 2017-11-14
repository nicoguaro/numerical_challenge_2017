from __future__ import division, print_function
from numpy import linspace, cos, pi, exp, sum


def clenshaw_curtis(fun, N=50, a=-1, b=1):
    nodes = linspace(-1, 1, N + 1)*pi
    jac = 0.5*(b - a)
    tfun = lambda x: fun(a + jac*(x + 1))
    inte = 0
    for k in range(0, N + 1, 2):
        coef = 1/N*(tfun(1) + tfun(-1)*(-1)**k +\
            2*sum(tfun(cos(nodes[1:-1]))*cos(k*nodes[1:-1])))
        if k == 0:
            inte += coef
        elif k == N:
            inte += coef/(1 - N**2)
        else:
            inte += 2*coef/(1 - k**2)
    return inte*jac


N = 100
fun = lambda x: exp(-x**2)
print(clenshaw_curtis(fun, N=N, a=0, b=3))
