from __future__ import division, print_function
import numpy as np
from numpy import array
from numpy.linalg import norm


def nelder_mead_step(fun, verts, alpha=1, gamma=2, rho=0.5,
                     sigma=0.5):
    nverts, _ = verts.shape
    f = np.apply_along_axis(fun, 1, verts)
    # 1. Order
    order = np.argsort(f)
    verts = verts[order, :]
    f = f[order]
    # 2. Calculate xo, the centroid"
    xo = verts[:-1, :].mean(axis=0)
    # 3. Reflection
    xr = xo + alpha*(xo - verts[-1, :])
    fr = fun(xr)
    if f[0]<=fr and fr<f[-2]:
        new_verts = np.vstack((verts[:-1, :], xr))
    # 4. Expansion
    elif fr<f[0]:
        xe = xo + gamma*(xr - xo)
        fe = fun(xe)
        if fe < fr:
            new_verts = np.vstack((verts[:-1, :], xe))
        else:
            new_verts = np.vstack((verts[:-1, :], xr))
    # 5. Contraction
    else:
        xc = xo + rho*(verts[-1, :] - xo)
        fc = fun(xc)
        if fc < f[-1]:
            new_verts = np.vstack((verts[:-1, :], xc))
    # 6. Shrink
        else:
            new_verts = np.zeros_like(verts)
            new_verts[0, :] = verts[0, :]
            for k in range(1, nverts):
                new_verts[k, :] = sigma*(verts[k,:] - verts[0,:])

    return new_verts


def nelder_mead(fun, x, niter=200, atol=1e-8, verbose=False):
    msg = "Maximum number of iterations reached."
    f_old = fun(x.mean(0))
    for cont in range(niter):
        if verbose:
            print("n: {}, x: {}".format(cont, x.mean(0)))
        x = nelder_mead_step(fun, x)
        df = fun(x.mean(0)) - f_old
        f_old = fun(x.mean(0))
        if norm(df) < atol:
            msg = "Extremum found with desired accuracy."
            break
    return x.mean(0), f_old, msg


def rosen(x):
    return (1 - x[0])**2 + 100*(x[1] - x[0]**2)**2


x = array([[1, 0],
           [1, 1],
           [2, 0]])
print(nelder_mead(rosen, x))
