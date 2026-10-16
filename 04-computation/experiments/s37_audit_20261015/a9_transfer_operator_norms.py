#!/usr/bin/env python3
"""S37 audit, part 9: THM-4658 statement 5 -- the branch norms are exact, but the sum of squared branch norms phi(sigma)
is not a bound for the norm of P*_rho = -2^-rho U - (3/2)^rho V.  Exact finite-level matrices: P* maps functions of
y mod 3^Q to functions of y mod 3^(Q+1); with normalised Haar measure ||P*|_Q|| = sqrt(lambda_max(M^T M)/3).
Reports, for real rho in {1/2, 3/4, 1}, the norm on level-Q functions and on mean-zero level-Q functions, against
phi(sigma) and the triangle bound 2^-sigma (1 + 3^(sigma - 1/2)).
"""
import numpy as np


def matrix(Q, rho):
    M3 = 3 ** Q
    M = np.zeros((3 * M3, M3))
    for y in range(3 * M3):
        M[y, (2 * y) % M3] += -2.0 ** (-rho)
        if y % 3 == 2:
            M[y, ((2 * y - 1) // 3) % M3] += -(1.5 ** rho)
    return M


for rho in (0.5, 0.75, 1.0):
    phi = 4 ** (-rho) * (1 + 3 ** (2 * rho - 1))
    tri = 2 ** (-rho) * (1 + 3 ** (rho - 0.5))
    row = []
    for Q in range(1, 7):
        M = matrix(Q, rho)
        G = M.T @ M / 3.0
        nrm = np.sqrt(np.linalg.eigvalsh(G).max())
        n = M.shape[1]
        Pm = np.eye(n) - np.ones((n, n)) / n          # projector onto mean-zero functions
        nrm0 = np.sqrt(np.linalg.eigvalsh(Pm @ G @ Pm).max())
        row.append((Q, round(float(nrm), 4), round(float(nrm0), 4)))
    print(f'rho = {rho}: phi = {phi:.4f}, sqrt(phi) = {np.sqrt(phi):.4f}, triangle bound {tri:.4f}; '
          f'(Q, ||P*|| on level Q, on mean-zero level Q): {row}')
