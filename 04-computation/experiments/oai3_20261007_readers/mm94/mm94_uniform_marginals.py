#!/usr/bin/env python3
"""Does the support {(i,j,i+j)} of C(a,b) carry a probability law with all three marginals uniform?
If yes, every quantum functional F_theta (CVZ) takes the value a^th1 b^th2 (a+b-1)^th3 on C(a,b), so all
known spectral points have the flattening profile P = sqrt(ab(a+b-1)) on polynomial multiplication."""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
bad = []
for a in range(1, 31):
    for b in range(a, 31):
        n = a * b
        I = np.repeat(np.arange(a), b); J = np.tile(np.arange(b), a); K = I + J
        rows, cols, rhs = [], [], []
        r = 0
        for m, arr, size in ((0, I, a), (1, J, b), (2, K, a + b - 1)):
            for v in range(size):
                for p in np.nonzero(arr == v)[0]:
                    rows.append(r); cols.append(p)
                rhs.append(1.0 / size); r += 1
        A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(r, n))
        res = linprog(np.zeros(n), A_eq=A, b_eq=np.array(rhs), bounds=(0, None), method="highs")
        if res.status != 0:
            bad.append((a, b))
print(f"uniform-marginal law exists for all 1 <= a <= b <= 30 except: {bad}")
# explicit law for a = 2 (Collatz slice): P(0,j) = (b-j)/(b(b+1)), P(1,j) = (j+1)/(b(b+1))
from fractions import Fraction as Fr
ok = True
for b in range(1, 60):
    P = {(0, j): Fr(b - j, b * (b + 1)) for j in range(b)}
    P.update({(1, j): Fr(j + 1, b * (b + 1)) for j in range(b)})
    ok &= all(sum(v for (i, j), v in P.items() if i == x) == Fr(1, 2) for x in (0, 1))
    ok &= all(sum(v for (i, j), v in P.items() if j == y) == Fr(1, b) for y in range(b))
    ok &= all(sum(v for (i, j), v in P.items() if i + j == z) == Fr(1, b + 1) for z in range(b + 1))
print(f"explicit a=2 law P(0,j)=(b-j)/(b(b+1)), P(1,j)=(j+1)/(b(b+1)) has uniform marginals for b<60: {ok}")
