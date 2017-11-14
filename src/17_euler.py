from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def euler(dydt, y0, t, args=()):
    ndof = len(y0)
    ntimes = len(t)
    y = np.zeros((ndof, ntimes))
    y[:, 0] = y0
    for cont in range(1, ntimes):
        h = t[cont] - t[cont - 1]
        y[:, cont] = y[:, cont - 1] + h*dydt(y[:, cont - 1], t[cont], *args)
    return y


def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return np.array(dydt)


b = 0.25
c = 5.0
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 101)
y = euler(pend, y0, t, args=(b, c))
plt.plot(t, y[0, :])
plt.plot(t, y[1, :])
plt.xlabel(r"$t$")
plt.legend([r"$\theta(t)$", r"$\omega(t)$"])
plt.show()