from __future__ import division, print_function
import numpy as np


def simps(y, x=None, dx=1.0):
    n = len(y)
    if x is None:
        inte = np.sum(y[0:n-2:2] + 4*y[1:n-1:2] + y[2:n:2])*dx
    else:
        dx = x[1::2] - x[:-1:2]
        inte = np.dot(y[0:n-2:2] + 4*y[1:n-1:2] + y[2:n:2], dx)
    return inte/3


n = 21
x = np.linspace(0, 10, n)
y = np.sin(x)
print(simps(y, x=x))
print(1 - np.cos(10))
