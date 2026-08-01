#!/usr/bin/env python3
"""Towards a Gregory-Newton certificate: which constraints force non-compensation?

Adjacent compensation needs  c^(2m) = -c^(m) * (1+u)^(2m)  to sit inside its
box.  Writing S_t = sum_j c_j binom(2m, j+t), the constraint at index 2m+t is
|S_t| <= N^(2m)_(2m+t), and the three indices around the middle are the
tightest the row has:

    t = -1 :  |S_-1| <= binom(2m,1)  = 2m
    t =  0 :  |S_0 | <= 2                      (the two-word middle class)
    t = +1 :  |S_1 | <= binom(2m,1)  = 2m

together with the shell-m box |c_j| <= N^(m)_j, the parity c_j = N^(m)_j (2),
and sum_j c_j = 0 (forced, THM/HYP-9075 sec 2).

QUESTION.  Do those THREE constraints alone already force c = 0?  If so the
certificate to be made positive in m is a 3-row object, not the whole 4m-row
system -- exactly the shape THM-2922 makes positive in its translation
parameter.  Tested here by LP over the RELAXED (rational) box: if even the
relaxation forces c = 0, the integer statement follows a fortiori and the
certificate is a genuine one.
"""
import sys
from math import comb

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from amm12592_shell_imbalance_module import N_vector


def middle_rows(m):
    """Rows S_t = sum_j c_j binom(2m, j+t) and their bounds, t = -1,0,1."""
    n = 2 * m + 1                      # c_0..c_2m
    rows, bnds = [], []
    N2 = N_vector(2 * m)
    for t in (-1, 0, 1):
        row = [comb(2 * m, j + t) if 0 <= j + t <= 2 * m else 0 for j in range(n)]
        rows.append(row)
        bnds.append(N2[2 * m + t])
    return rows, bnds


def forced_zero(m, extra_rows=0):
    """Is c = 0 the only point of the relaxed polytope?  Maximise |c_j| for each j."""
    n = 2 * m + 1
    N1 = N_vector(m)
    rows, bnds = middle_rows(m)
    if extra_rows:
        N2 = N_vector(2 * m)
        for t in list(range(-1 - extra_rows, -1)) + list(range(2, 2 + extra_rows)):
            idx = 2 * m + t
            if 0 <= idx <= 4 * m:
                rows.append([comb(2 * m, j + t) if 0 <= j + t <= 2 * m else 0
                             for j in range(n)])
                bnds.append(N2[idx])
    A_ub, b_ub = [], []
    for row, b in zip(rows, bnds):
        A_ub.append(row);        b_ub.append(b)
        A_ub.append([-v for v in row]); b_ub.append(b)
    A_eq = [[1.0] * n]
    b_eq = [0.0]
    bounds = [(-N1[j], N1[j]) for j in range(n)]
    best = 0.0
    for j in range(n):
        if N1[j] == 0:
            continue
        for s in (1, -1):
            c = [0.0] * n
            c[j] = -s
            r = linprog(c, A_ub=np.array(A_ub, float), b_ub=np.array(b_ub, float),
                        A_eq=np.array(A_eq, float), b_eq=np.array(b_eq, float),
                        bounds=bounds, method="highs")
            if r.status == 0:
                best = max(best, abs(r.x[j]))
    return best


if __name__ == "__main__":
    print("Do the three middle rows (t = -1,0,1) alone force c = 0?")
    print("  (LP over the RELAXED box; max_j max |c_j| on the polytope)")
    print()
    print("    m    3 rows      5 rows      7 rows")
    for m in (2, 4, 8, 16, 32, 64):
        v3 = forced_zero(m, 0)
        v5 = forced_zero(m, 1)
        v7 = forced_zero(m, 2)
        print(f"  {m:4d}   {v3:9.4f}   {v5:9.4f}   {v7:9.4f}")
    print()
    print("""READING.  A value of 0 in a column means those rows alone pin c = 0 over
the rationals, hence a fortiori over the integers, and the certificate needed
for all dyadic m is that small fixed row-set -- the analogue of THM-2922's
fixed degree-seven Macaulay chart.  A positive value means the relaxation has
room and more rows (or the integrality) are doing the work, so no fixed finite
row-set of this kind can be the certificate.""")
