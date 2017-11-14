from __future__ import division, print_function
import numpy as np
from numpy import log, sin, cos, arctan2, pi, mean, arctan
from numpy.linalg import norm, solve
import matplotlib.pyplot as plt


def assem(coords, elems):
    nelems = elems.shape[0]
    Gmat = np.zeros((nelems, nelems))
    Fmat = np.zeros((nelems, nelems))
    for ev_cont, elem1 in enumerate(elems):
        print(coords[elem1[0]], coords[elem1[1]])
        for col_cont, elem2 in enumerate(elems):
            pt_col = mean(coords[elem2], axis=0)
            if ev_cont == col_cont:
                L = norm(coords[elem1[1]] - coords[elem1[0]])
                Gmat[ev_cont, ev_cont] = - L/(2*pi)*(log(L/2) - 1)
                Fmat[ev_cont, ev_cont] = - 0.5
            else:
                Gij, Fij = influence_coeff(elem1, coords, pt_col)
                Gmat[ev_cont, col_cont] = Gij
                Fmat[ev_cont, col_cont] = Fij
    return Gmat, Fmat


def influence_coeff(elem, coords, pt_col):
    r_A = coords[elem[0]] - pt_col
    r_B = coords[elem[1]] - pt_col
    theta_A = arctan2(r_A[1], r_A[0])
    theta_B = arctan2(r_B[1], r_B[0])
    G_coeff = r_B[1]*(log(norm(r_B)) - 1) + theta_B*r_B[0] -\
              (r_A[1]*(log(norm(r_A)) - 1) + theta_A*r_A[0])
    F_coeff = theta_B - theta_A
    return -G_coeff/(2*pi), F_coeff/(2*pi)


def eval_sol(ev_coords, coords):
    nelems = elems.shape[0]
    Gmat = np.zeros((nelems, nelems))
    Fmat = np.zeros((nelems, nelems))
    for ev_cont, elem1 in enumerate(elems):
        L = norm(coords[elem1[1]] - coords[elem1[0]])
        for col_cont, elem2 in enumerate(elems):
            pt_col = mean(coords[elem2], axis=0)
            if ev_cont == col_cont:
                Gmat[ev_cont, ev_cont] = - L/(2*pi)*(log(L/2) - 1)
                Fmat[ev_cont, ev_cont] = - 0.5
            else:
                Gmat[ev_cont, col_cont], Fmat[ev_cont, col_cont] = \
                    influence_coeff(elem1, coords, pt_col)

nelems = 3
rad = 1.0
theta =  np.linspace(0, 2*pi, nelems, endpoint=False)
coords = rad * np.vstack((cos(theta), sin(theta))).T
elems = np.array([[cont, (cont + 1)%nelems] for cont in range(nelems)])
Gmat, Fmat = assem(coords, elems)
u_boundary = np.ones_like(theta)
q_boundary = solve(Gmat, Fmat.dot(u_boundary))
