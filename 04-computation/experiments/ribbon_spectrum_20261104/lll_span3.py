#!/usr/bin/env python3
"""lll_span3.py -- LLL-reduce the integer basis of the observable module of span 3 (rank 62) to obtain
small-entry transition matrices; writes span3_intrep.txt in the same format as intrep.py."""
import sys
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
from intrep import build, vecmat, in_module_coeffs, echelon
from ribbon import letters
w = 3
basis, A, y0, hcoef, order, states, mats, final = build(w, verbose=False)
r = len(basis); d = len(order); L = letters(w)
print("rank", r, "largest basis entry", max(abs(x) for row in basis for x in row), flush=True)
dm = DomainMatrix([[ZZ(x) for x in row] for row in basis], (r, d), ZZ)
red = dm.lll()
rows = [[int(x) for x in row] for row in red.to_Matrix().tolist()]
print("after LLL: largest basis entry", max(abs(x) for row in rows for x in row), flush=True)
# re-echelon for the coefficient solver (echelon form keeps the module, may regrow entries), so instead
# solve coefficients by exact linear algebra on the LLL basis
from sympy import Matrix, Rational
Bm = Matrix(rows)                      # r x d
BBt = Bm * Bm.T
BBt_inv = BBt.inv()
def coeffs(v):
    c = BBt_inv * (Bm * Matrix(v))
    assert all(x.q == 1 for x in c), "non-integer coefficient: not in module?"
    assert Bm.T * c == Matrix(v)
    return [int(x) for x in c]
phi = [0] * d
for j in final: phi[j] = 1
Anew = {a: [coeffs(vecmat(f, mats[a], d)) for f in rows] for a in L}
y0n = [f[0] for f in rows]
hc = coeffs(phi)
print("largest |entry| of the reduced matrices:", max(abs(x) for a in L for row in Anew[a] for x in row), flush=True)
with open("span3_intrep.txt", "w") as f:
    f.write(f"{r} {len(L)}\n")
    f.write(" ".join(map(str, y0n)) + "\n")
    f.write(" ".join(map(str, hc)) + "\n")
    for a in L:
        mask = sum(1 << (s - 1) for s in a)
        f.write(f"{mask}\n")
        for k in range(r):
            f.write(" ".join(map(str, Anew[a][k])) + "\n")
print("written span3_intrep.txt")
