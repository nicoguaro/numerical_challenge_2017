from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def vander_mat(x):
    n = len(x)
    van = np.zeros((n, n))
    power = np.array(range(n))
    for row in range(n):
        van[row, :] = x[row]**power
    return van


def conf_vander_mat(x):
    n = len(x)
    conf_van = np.zeros((2*n, 2*n))
    power = np.array(range(2*n))
    for row in range(n):
        conf_van[row, :] = x[row]**power
        conf_van[row + n, :] = power*x[row]**(power - 1)
    return conf_van


def inter_coef(x, inter_type="lagrange"):
    if inter_type == "lagrange":
        vand_mat = vander_mat(x)
    elif inter_type == "hermite":
        vand_mat = conf_vander_mat(x)
    coef = np.linalg.solve(vand_mat, np.eye(vand_mat.shape[0]))
    return coef


def compute_interp(x, f, x_eval, df=None):
    n = len(x)
    if df is None:
        coef = inter_coef(x, inter_type="lagrange")
    else:
        coef = inter_coef(x, inter_type="hermite")
    f_eval = np.zeros_like(x_eval)
    nmat = coef.shape[0]
    for row in range(nmat):
        for col in range(nmat):
            if col < n or nmat == n:
                f_eval += coef[row, col]*x_eval**row*f[col]
            else:
                f_eval += coef[row, col]*x_eval**row*df[col - n]
    return f_eval


case = "comparison"
if case == "run":
    n = 7
    x = -np.cos(np.linspace(0, np.pi, n))
    f = lambda x: 1/(1 + 25*x**2)
    df = lambda x: -50*x/(1 + 25*x**2)**2
    x_eval = np.linspace(-1, 1, 500)
    interp_f = compute_interp(x, f(x), x_eval, df=df(x))
    plt.plot(x_eval, f(x_eval))
    plt.plot(x_eval, interp_f)
    plt.plot(x, f(x), ".")
    plt.ylim(0, 1.2)
    plt.show()

elif case == "comparison":
    n_dof = np.array(range(1, 20))
    error_herm = np.zeros(19)
    error_lag = np.zeros(19)
    for cont, n in enumerate(n_dof):
        f = lambda x: 1/(1 + 25*x**2)
        df = lambda x: -50*x/(1 + 25*x**2)**2
        x = -np.cos(np.linspace(0, np.pi, n))
        x2 = -np.cos(np.linspace(0, np.pi, 2*n))
        x_eval = np.linspace(-1, 1, 500)
        herm = compute_interp(x, f(x), x_eval, df=df(x))
        lag = compute_interp(x2, f(x2), x_eval)
        fun = f(x_eval)
        error_herm[cont] = np.linalg.norm(fun - herm)/np.linalg.norm(fun)
        error_lag[cont] = np.linalg.norm(fun - lag)/np.linalg.norm(fun)

    plt.plot(2*n_dof, error_lag)
    plt.plot(2*n_dof, error_herm)
    plt.xlabel("Number of degrees of freedom")
    plt.ylabel("Relative error")
    plt.legend(["Lagrange", "Hermite"])
    plt.show()