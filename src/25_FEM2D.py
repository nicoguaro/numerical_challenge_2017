from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

def assem(coords, elems, source):
    ncoords = coords.shape[0]
    stiff = np.zeros((ncoords, ncoords))
    rhs = np.zeros((ncoords))
    for el_cont, elem in enumerate(elems):
        stiff_loc, jaco = local_stiff(coords[elem])
        rhs[elem] += jaco*np.mean(source[elem])
        for row in range(3):
            for col in range(3):
                row_glob, col_glob = elem[row], elem[col]
                stiff[row_glob, col_glob] += stiff_loc[row, col]
    return stiff, rhs


def local_stiff(coords):
    dHdr = np.array([[-1, 1, 0], [-1, 0, 1]])
    jaco = dHdr.dot(coords)
    stiff = 0.5*np.array([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]])
    return stiff/np.linalg.det(jaco), np.linalg.det(jaco)


sq3 = np.sqrt(3)
coords = np.array([
        [sq3, -1],
        [0, 0],
        [2*sq3, 0],
        [0, 2],
        [2*sq3, 2],
        [sq3, 3],
        [sq3, 1]])
elems = np.array([
        [1, 0, 6],
        [0, 2, 6],
        [2, 4, 6],
        [4, 5, 6],
        [5, 3, 6],
        [3, 1, 6]])
source = -np.ones(7)
stiff, rhs = assem(coords, elems, source)
free = range(6)
sol = np.linalg.solve(stiff[np.ix_(free, free)], rhs[free])
sol_c = np.zeros(coords.shape[0])
sol_c[free] = sol
plt.tricontourf(coords[:, 0], coords[:, 1], sol_c, cmap="hot")
plt.colorbar()
plt.axis("image")
plt.savefig("FEM2D.svg")
plt.show()