from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def RK4(dydt, y0, t, args=()):
    ndof = len(y0)
    ntimes = len(t)
    y = np.zeros((ndof, ntimes))
    y[:, 0] = y0
    for cont in range(1, ntimes):
        h = t[cont] - t[cont - 1]
        k1 = dydt(y[:, cont - 1], t[cont], *args)
        k2 = dydt(y[:, cont - 1] + 0.5*h*k1, t[cont] + 0.5*h, *args)
        k3 = dydt(y[:, cont - 1] + 0.5*h*k2, t[cont] + 0.5*h, *args)
        k4 = dydt(y[:, cont - 1] + h*k3, t[cont] + h, *args)
        y[:, cont] = y[:, cont - 1] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return y


def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return np.array(dydt)


b = 0.25
c = 5.0
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 1001)
y = RK4(pend, y0, t, args=(b, c))
plt.plot(t, y[0, :])
plt.plot(t, y[1, :])
plt.xlabel(r"$t$")
plt.legend([r"$\theta(t)$", r"$\omega(t)$"])
plt.show()