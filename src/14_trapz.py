from __future__ import division, print_function
import numpy as np


def trapz(y, x=None, dx=1.0):
    if x is None:
        inte = 0.5*dx*(2*np.sum(y[1:-1]) + y[0] + y[-1])
    else:
        dx = x[1:] - x[:-1]
        inte = 0.5*np.dot(y[:-1] + y[1:], dx)
    return inte


dx = 0.001
x = np.arange(0, 10, dx)
y = np.sin(x)
print(trapz(y, dx=dx))
print(1 - np.cos(10))
