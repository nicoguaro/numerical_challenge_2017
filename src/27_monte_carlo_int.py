from __future__ import division, print_function
import numpy as np


def monte_carlo_int(fun, N, low, high, args=()):
    ndims = len(low)
    pts = np.random.uniform(low=low, high=high, size=(N, ndims))
    V = np.prod(np.asarray(high) - np.asarray(low))
    return V*np.sum(fun(pts, *args))/N


def circ(x, rad):
    return 0.5*(1 - np.sign(x[:, 0]**2 + x[:, 1]**2 - rad**2))


N = 1000000
low = [-1, -1]
high = [1, 1]
rad = 1
inte = monte_carlo_int(circ, N, low, high, args=(rad,))