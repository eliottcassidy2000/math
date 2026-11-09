#!/usr/bin/env python3
"""rankdrop.py -- where does the minimal dimension of the span-2 packing series u -> I(Omega(T_u), x) drop?
Generic rank over Q(x) of the module of functionals phi.M_v(x) (20-state automaton), the polynomial determinant of a
nonsingular 14 x 14 minor, its factorisation, and the exact rank at every root (x = 0 and x = 2 expected)."""
import sys, os, itertools, random
sys.path.insert(0, '../ribbon_spectrum_20261104')
import sympy as sp
from ribbon import letters
from ribbon_packing import packing_system
x = sp.symbols('x')
w = 2
order, states, mats, init, final = packing_system(w)
d = len(order); L = letters(w)
M = {a: sp.zeros(d, d) for a in L}
for a in L:
    for (j, i), (c0, c1) in mats[a].items(): M[a][j, i] = c0 + c1 * x
phi = sp.zeros(1, d); phi[0, final] = 1
rows = [phi]; frontier = [phi]
for depth in range(1, 5):
    frontier = [f * M[a] for f in frontier for a in L]
    rows += frontier
F = sp.Matrix.vstack(*rows)
F = F.applyfunc(sp.expand)
# generic rank via a random rational specialisation (sanity) then exact
rk_generic = F.subs(x, sp.Rational(7, 3)).rank()
print("rows", F.shape, "generic rank (x = 7/3):", rk_generic)
for xv in [-3, -2, -1, 0, sp.Rational(1, 2), 1, sp.Rational(3, 2), 2, sp.Rational(5, 2), 3, 4, 5, 6, 7, 10]:
    print(f"  rank at x = {xv}: {F.subs(x, xv).rank()}")
# an independent set of 14 rows at a generic point, then a 14 x 14 minor with nonzero polynomial determinant
Fq = F.subs(x, sp.Rational(7, 3))
_, piv_rows = Fq.T.rref()   # independent rows
R = F.extract(list(piv_rows), list(range(d)))
_, piv_cols = R.subs(x, sp.Rational(7, 3)).rref()
S = R.extract(list(range(R.shape[0])), list(piv_cols))
det = sp.factor(S.det())
print("a nonsingular 14x14 minor has determinant", det)
roots = sp.roots(sp.Poly(sp.expand(S.det()), x))
print("its roots:", roots)
for rt in roots:
    print(f"  exact rank at root x = {rt}: {F.subs(x, rt).rank()}")
# the relation at x = 2: left kernel of the 14 independent rows
K = R.subs(x, 2).T.nullspace()
print("relations among the 14 generic basis functionals at x = 2:", len(K))
