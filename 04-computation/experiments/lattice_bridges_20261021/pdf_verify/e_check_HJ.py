#!/usr/bin/env python3
"""Check Prop. 2.2 properties for reconstructed H_J: integral over O_F, det = 7^6, H_J = sqrt(-7) H_1 with H_1 in GL_12(O_F)."""
import json, sys
from fractions import Fraction as Fr
import sympy as sp
s = sp.sqrt(-7)
for cand in map(int, sys.argv[1:]):
    raw = json.load(open(f'HJ_cand{cand}.json'))
    HJ = sp.Matrix(12, 12, lambda i, j: (sp.Rational(raw[i][j][0]) + sp.Rational(raw[i][j][1])*s)/2)
    herm = sp.simplify(HJ - HJ.H) == sp.zeros(12)
    d = sp.nsimplify(sp.expand(HJ.det(method='berkowitz')))
    H1 = (HJ/s).applyfunc(lambda e: sp.nsimplify(sp.expand(e)))
    def in_OF(e):
        a, b = sp.re(e), sp.im(e)/sp.sqrt(7)
        A, Bq = 2*a, 2*b
        return A.is_integer and Bq.is_integer and (A - Bq) % 2 == 0
    integral1 = all(in_OF(e) for e in H1)
    d1 = sp.nsimplify(sp.expand(H1.det(method='berkowitz')))
    print(f"candidate {cand}: Hermitian {herm}; det H_J = {d} (7^6 = {7**6}); H_1 = H_J/sqrt(-7) has entries in O_F: {integral1}; det H_1 = {d1}; skew-Hermitian: {sp.simplify(H1 + H1.H) == sp.zeros(12)}")
