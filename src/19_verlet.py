from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def verlet(force, x0, v0, m, t, args=()):
    ndof = len(x0)
    ntimes = len(t)
    x = np.zeros((ndof, ntimes))
    dt = t[1] - t[0]
    x[:, 0] = x0
    F = force(x0, v0, m, t, *args)
    x[::2, 1] = x0[::2] + v0[::2]*dt + 0.5*F[::2]*dt**2/m
    x[1::2, 1] = x0[1::2] + v0[1::2]*dt + 0.5*F[1::2]*dt**2/m
    for cont in range(2, ntimes):
        dt = t[cont] - t[cont - 1]
        v = (x[:, cont - 1] - x[:, cont - 2])/dt
        acel = force(x[:, cont - 1], v, m, t, *args)
        acel[::2] = acel[::2]/m
        acel[1::2] = acel[1::2]/m
        x[:, cont] = 2*x[:, cont - 1] - x[:, cont - 2] + acel*dt**2
    return x


spring_force = lambda x, v, m, t, k: -k*x
x0 = np.array([1.0, 0.0])
v0 = np.array([0.0, 1.0])
m = np.array([1.0])
t = np.linspace(0, 10.0, 1000)
k = 1.0
x = verlet(spring_force, x0, v0, m, t, args=(k,))

#%% Plot
plt.figure(figsize=(6, 3))
plt.subplot(121)
plt.plot(t, x[0, :])
plt.plot(t, x[1, :])
plt.subplot(122)
plt.plot(x[0, :], x[1, :])
plt.show()