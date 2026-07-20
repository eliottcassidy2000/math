#!/usr/bin/env python3
"""even_n_perron_deficiency_kps_S128c119.py -- kind-pasteur-2026-07-20-S128c119

THM-1555 I.c says the Perron deficiency (n-1)/2 - rho EQUALS the total height of the
spectrum above the line Re = -1/2.  At EVEN n no regular tournament exists, so the
deficiency is strictly positive and its MINIMUM is a genuine constant of each even n.
This identifies those constants exactly.

Method: exhaustive over all tournaments at n = 4, 6 (up to the cap), take the argmax of
rho, then get its characteristic polynomial with INTEGER coefficients (via the
integer-valued elementary symmetric functions of a 0/1 matrix) and factor it over Q.
"""
import sys
from itertools import product
import numpy as np
import sympy as sp

lam = sp.Symbol('lam')

def from_bits(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1: A[i, j] = 1
            else:               A[j, i] = 1
            idx += 1
    return A

for n, cap in ((4, 1 << 6), (5, 1 << 10), (6, 1 << 15), (7, 200000)):
    best, bestA = -1.0, None
    for bits in range(min(1 << (n * (n - 1) // 2), cap)):
        A = from_bits(bits, n)
        r = max(np.linalg.eigvals(A.astype(float)).real)
        if r > best:
            best, bestA = r, A
    M = sp.Matrix(bestA.tolist())
    cp = sp.Poly(M.charpoly(lam).as_expr(), lam)
    fac = sp.factor(cp.as_expr())
    scores = sorted(bestA.sum(axis=1).tolist())
    # minimal polynomial of the Perron root = the irreducible factor vanishing at rho
    minpoly = None
    for f, _ in sp.factor_list(cp.as_expr())[1]:
        rts = sp.Poly(f, lam).nroots()
        if any(abs(complex(r) - best) < 1e-8 for r in rts):
            minpoly = sp.expand(f)
    print("n = %d  exhaustive: %s" % (n, (1 << (n*(n-1)//2)) <= cap))
    print("   max rho          = %.10f      (n-1)/2 = %s" % (best, sp.Rational(n-1, 2)))
    print("   deficiency       = %.10f" % ((n - 1) / 2 - best))
    print("   score sequence   = %s" % scores)
    print("   char poly        = %s" % cp.as_expr())
    print("   factored         = %s" % fac)
    print("   min poly of rho  = %s   (degree %d)"
          % (minpoly, sp.Poly(minpoly, lam).degree() if minpoly is not None else -1))
    print()
