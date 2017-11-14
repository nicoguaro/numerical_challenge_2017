from __future__ import division, print_function
import numpy as np
from scipy.linalg import eigh as eigsh
import matplotlib.pyplot as plt


def beam_FDM_eigs(L, N): 
    x = np.linspace(0, L, N)
    dx = x[1] - x[0]
    stiff = np.diag(6*np.ones(N - 2)) -\
            np.diag(4*np.ones(N - 3), -1) - np.diag(4*np.ones(N - 3), 1) +\
            np.diag(1*np.ones(N - 4), 2) + np.diag(1*np.ones(N - 4), -2)
    vals, vecs = eigsh(stiff/dx**4)     
    return vals, vecs, x


N = 1001
nvals = 20
nvecs = 4
vals, vecs, x = beam_FDM_eigs(1.0, N)

#%% Plotting
num = np.linspace(1, nvals, nvals)
plt.rcParams["mathtext.fontset"] = "cm"
plt.figure(figsize=(8, 3))
plt.subplot(1, 2, 1)
plt.plot(num, np.sqrt(vals[0:nvals]), "o")
plt.xlabel(r"$N$")
plt.ylabel(r"$\omega\sqrt{\frac{\lambda}{EI}}$")
plt.subplot(1, 2 ,2)
for k in range(nvecs):
    vec = np.zeros(N)
    vec[1:-1] = vecs[:, k]
    plt.plot(x, vec, label=r'$n=%i$'%(k+1))

plt.xlabel(r"$x$")
plt.ylabel(r"$w$")
plt.legend(ncol=2, framealpha=0.8, loc=1)
plt.tight_layout()
plt.show()
