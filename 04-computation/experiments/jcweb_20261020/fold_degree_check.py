#!/usr/bin/env python3
"""Corollary R' check on THM-1300: (i) generic target: number of complex preimages under F (expect 3 = generic degree of the
quotient map); (ii) a target on the fixed hyperplane X1 = 0 with Y != 0: preimages come from the fixed plane x = 0 (one) and
from orbits over {l = 0} where the orbit map x -> (A/x, B/x^2) is 1-to-1 generically (gcd(1,2) = 1) -- expect 1 + 1 = 2... 
(iii) the special target (-1/4, 0, 0) with Y = 0: 1 + 2 = 3 (fold x -> B/x^2)."""
import sympy as sp
x, y, z = sp.symbols('x y z')
u = 1 + x*y
F = [u**3*z + y**2*u*(4 + 3*x*y), y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y), 2*x - 3*x**2*y - x**3*z]
def preimages(target):
    eqs = [sp.expand(F[i] - target[i]) for i in range(3)]
    G = sp.groebner(eqs, x, y, z, order='lex')
    # count solutions via the univariate polynomial in z (last) if zero-dimensional
    polys = list(G.exprs)
    uni = [p for p in polys if p.free_symbols <= {z}]
    sols = sp.solve(eqs, [x, y, z], dict=True)
    return len(sols), [tuple(sp.nsimplify(s[v]) for v in (x, y, z)) for s in sols]
for name, tgt in [("generic (3, 5, 7)", (3, 5, 7)), ("generic (1, 2, -1)", (1, 2, -1)), ("fixed-plane target (2, 1, 0)", (2, 1, 0)), ("fixed-plane target (5, -2, 0)", (5, -2, 0)), ("the triple collision (-1/4, 0, 0)", (sp.Rational(-1, 4), 0, 0)), ("fixed-plane Y = 0 target (1, 0, 0)", (1, 0, 0))]:
    n, sols = preimages(tgt)
    print(f"{name}: {n} preimage(s)", "" if n > 3 else sols)
