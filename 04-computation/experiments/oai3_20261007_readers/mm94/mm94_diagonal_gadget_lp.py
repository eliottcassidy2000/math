#!/usr/bin/env python3
"""Least diagonal growth under (S),(B),(C) plus ONE diagonal-only gadget P(a, ceil(u a)) >= m P(a,a).
Prediction (mm94 tangent/fake-profile argument): growth exponent kappa = 2(m-1)/(u-1), i.e. the method
proves at best omega <= 3(u-1)/(2(m-1)).  Flattening ranks force u(u+1) >= 2 m^2 for such a gadget."""
import numpy as np, scipy.sparse as sp, math
from scipy.optimize import linprog
def least_diag(N, m, u, ns):
    idx = {}
    for a in range(1, N + 1):
        for b in range(a, N + 1):
            idx[(a, b)] = len(idx)
    I = lambda a, b: idx[(min(a, b), max(a, b))]
    rows, cols, vals, rhs, r = [], [], [], [], 0
    def add(coefs, bound):
        nonlocal r
        for (a, b), c in coefs:
            rows.append(r); cols.append(I(a, b)); vals.append(c)
        rhs.append(bound); r += 1
    for a in range(1, N + 1):
        for b in range(2, N):
            add([((a, b + 1), 1.0), ((a, b - 1), 1.0), ((a, b), -2.0)], 0.0)
    for a in range(1, N + 1):
        B = math.ceil(u * a)
        if B <= N:
            add([((a, a), float(m)), ((a, B), -1.0)], 0.0)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, len(idx)))
    Ae = sp.csr_matrix((np.ones(N), (np.arange(N), [I(1, b) for b in range(1, N + 1)])), shape=(N, len(idx)))
    be = np.arange(1, N + 1, dtype=float)
    out = {}
    for n in ns:
        c = np.zeros(len(idx)); c[I(n, n)] = 1.0
        res = linprog(c, A_ub=A, b_ub=np.zeros(r), A_eq=Ae, b_eq=be, bounds=(0, None), method="highs")
        out[n] = res.fun
    return out
for m, u in ((3, 4.0), (2, 2.5), (2, 2.4), (3, 3.8)):
    ns = (10, 20, 30)
    N = math.ceil(u * 30) + 1
    v = least_diag(N, m, u, ns)
    k1 = math.log(v[20] / v[10]) / math.log(2); k2 = math.log(v[30] / v[20]) / math.log(1.5)
    pred = 2 * (m - 1) / (u - 1)
    flat_ok = u * (u + 1) >= 2 * m * m
    print(f"m={m} u={u}: min P(n,n) = {[round(v[n], 3) for n in ns]}; local exponents {k1:.4f}, {k2:.4f}; "
          f"predicted kappa = {pred:.4f} -> omega <= {3/pred:.4f}; flattening-consistent: {flat_ok}")
