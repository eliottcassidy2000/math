#!/usr/bin/env python3
"""exploration: spectra of the letter matrices of the minimal representation; affine-span dimensions by suffix."""
import itertools
import sys
from fractions import Fraction

import sympy as sp

import minrep
from lattice_tools import rank_q

rep = minrep.build()
R = rep["R"]
A = rep["A"]
L0 = rep["L0"]
x = sp.symbols("x")
for c in range(4):
    Mx = sp.Matrix(A[c])
    print(f"letter {c}: characteristic polynomial", sp.factor(Mx.charpoly(x).as_expr()), " rank", Mx.rank())
for d in [(0, 1), (1, 2), (0, 3), (1, 2, 0), (2, 3, 3, 0)]:
    Mx = sp.eye(R)
    for c in d:
        Mx = Mx * sp.Matrix(A[c])
    print(f"block {d}: characteristic polynomial", sp.factor(Mx.charpoly(x).as_expr()), " rank", Mx.rank())

# direction space ker rho: differences of left vectors
words = [w for n in range(0, 4) for w in itertools.product(range(4), repeat=n)]
lv = {w: rep["l_of"](w) for w in words}
base = lv[()]
diffs = [[a - b for a, b in zip(lv[w], base)] for w in words if w]
# basis of the direction space
basis = []
for dv in diffs:
    if rank_q(basis + [dv]) > len(basis):
        basis.append(dv)
print("dimension of the affine span of all left vectors:", len(basis))


def img_rank(word):
    Mx = [[int(i == j) for j in range(R)] for i in range(R)]
    rows = basis
    for c in word:
        rows = [[sum(r[i] * A[c][i][j] for i in range(R)) for j in range(R)] for r in rows]
    return rank_q(rows)


for k in range(1, 5):
    tab = {}
    for s in itertools.product(range(4), repeat=k):
        tab[s] = img_rank(s)
    vals = sorted(set(tab.values()))
    print(f"suffix length {k}: dimensions of the affine spans of l(w s): "
          + ", ".join(f"{v} (x{sum(1 for t in tab.values() if t == v)})" for v in vals))
    if k <= 2:
        for s, v in sorted(tab.items()):
            print("    ", "".join(map(str, s)), v)
