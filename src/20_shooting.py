from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.integrate import odeint


def shooting(dydx, x, x0, xf, shoot=None):
    if shoot is None:
        shoot = np.random.uniform(-20, 20)
    F = lambda s, x0, xf, x: odeint(dydx, [x0, s], x)[-1, 0] - xf
    shoot = newton(F, shoot, args=(x0, xf, x))
    y = odeint(dydx, [x0, shoot], x)
    return y[:, 0], shoot


func = lambda y, t: [y[1], 1.5*y[0]**2]
x = np.linspace(0, 1, 1000)
y, shoot = shooting(func, x, 4, 1, shoot=-5)
plt.plot(x, y)
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.show()