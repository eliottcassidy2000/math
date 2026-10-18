#!/usr/bin/env python3
"""Delsarte / Schrijver LP bound on the independence number of G_q = Cay(F_(q^2), mu_(q+1)) (mac-mini-2026-10-18).
By the symmetry z -> lambda z (lambda in mu_(q+1)) an optimal F may be taken constant on norm circles N(z) = c, so the
LP has q variables f(c), c in F_q (f(0) = F(0) = 1, f(1) = 0 on the connection set, f >= 0), and the Fourier condition
F_hat(w) = 1 + sum_(c != 0) f(c) K(c N(w)) >= 0 for all w != 0, where K(m) = sum_(N(z) = m) e(Tr(z)/q).
alpha(G_q) <= max sum_z F(z) = 1 + (q+1) sum_c f(c). Then chi >= q^2 / alpha.
The dual certificate (multipliers y_d >= 0 for the Fourier constraints) is printed and re-verified in exact-margin
floating point; the K(m) are real (closed under z -> -z) and computed with |error| < 1e-9."""
import sys, cmath, math
import numpy as np
from scipy.optimize import linprog
q = int(sys.argv[1])
r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
circ = {m: [] for m in range(q)}
for a in range(q):
    for b in range(q):
        circ[(a * a - r * b * b) % q].append((a, b))
assert all(len(circ[m]) == q + 1 for m in range(1, q))
K = {m: sum(cmath.exp(2j * math.pi * ((2 * a) % q) / q) for (a, b) in circ[m]).real for m in range(1, q)}   # Tr(a + b t) = 2a
cs = [c for c in range(2, q)]            # variables f(c), c != 0, 1
# maximize sum f(c)  s.t. for each d in 1..q-1:  -sum_c f(c) K(c d) <= 1 ; f >= 0
A = [[-K[(c * d) % q] for c in cs] for d in range(1, q)]
res = linprog(c=[-1.0] * len(cs), A_ub=A, b_ub=[1.0] * (q - 1), bounds=[(0, None)] * len(cs), method='highs')
alpha_lp = 1 + (q + 1) * (-res.fun)
print('q = %d: Delsarte/Schrijver LP bound alpha <= %.6f ; q^2/6 = %.3f, q^2/7 = %.3f ; chi >= q^2/alpha_LP = %.4f' % (q, alpha_lp, q * q / 6, q * q / 7, q * q / alpha_lp))
# dual: minimize sum_d y_d subject to sum_d y_d * K(c d) >= 1 ... (for each c); bound = 1 + (q+1) * sum y_d
dres = linprog(c=[1.0] * (q - 1), A_ub=[[-K[(c * d) % q] for d in range(1, q)] for c in cs], b_ub=[-1.0] * len(cs), bounds=[(0, None)] * (q - 1), method='highs')
y = dres.x
viol = max(1.0 - sum(y[d - 1] * K[(c * d) % q] for d in range(1, q)) for c in cs)
print('dual objective 1 + (q+1) sum y = %.6f ; max dual constraint violation %.3e (bound valid after scaling y by 1/(1 - viol))' % (1 + (q + 1) * y.sum(), max(viol, 0)))
yb = y / (1 - max(viol, 0) - 1e-9) if viol > -1 else y
print('certified alpha <= %.6f -> chi >= %.4f' % (1 + (q + 1) * yb.sum(), q * q / (1 + (q + 1) * yb.sum())))
