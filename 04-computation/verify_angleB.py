#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL verification of ANGLE B dissociation-closure claim.
Independent exact code (fractions). Does NOT import the claim's artifact.
"""
import sys, itertools
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb

# ---------------------------------------------------------------
# 1. g(t) per THM-534 (copy ONLY the math definition, re-derive y_r myself)
# ---------------------------------------------------------------
def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            val = F((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = F(-(t-2)*(t-3)*(t-6), 36)
        else:
            val = F((t-3)*(t-4), 12)
        g.append(val)
    return g

# ---------------------------------------------------------------
# 2. Re-derive y_r from g: g(t) = sum_r y_r * C(t,r).  This is the
#    inverse binomial (finite difference) transform: y_r = sum_{j=0}^{r} (-1)^{r-j} C(r,j) g(j).
# ---------------------------------------------------------------
def y_coeffs(k):
    g = g_poly(k)
    y = []
    for r in range(7):
        s = F(0)
        for j in range(r+1):
            s += (-1)**(r-j) * comb(r, j) * g[j]
        y.append(s)
    return y

# sanity: reconstruct g(t) from y to confirm y is correct
def check_y(k):
    g = g_poly(k); y = y_coeffs(k)
    for t in range(7):
        recon = sum(y[r]*comb(t, r) for r in range(7))
        assert recon == g[t], (k, t, recon, g[t])
    return y

# ---------------------------------------------------------------
# 3. The claimed closed form for the fully-independent limit:
#    L_inf = sum_r y_r C(6,r) (1-r/7)^(k-1)
#    Derivation check: S_r = E[C(N,r)] = sum_{|A|=r} J(A) ; in the fully
#    independent limit each of the k-1 free generators avoids r sectors
#    with prob (1-r/7) independently, J(A)=(1-r/7)^(k-1), and there are
#    C(6,r) sets A. So S_r^inf = C(6,r)(1-r/7)^(k-1), L_inf=sum y_r S_r^inf.
# ---------------------------------------------------------------
def L_inf(k):
    y = y_coeffs(k)
    return sum(y[r] * comb(6, r) * F(7-r, 7)**(k-1) for r in range(7))

claimed = {8: F(40573,823543), 9: F(2616535,17294403), 10: F(8830611,40353607)}
caps = {8: 0.3815, 9: 0.49426, 10: 0.6044}

print("=== y_r reconstruction and L_inf closed form ===")
for k in [8,9,10]:
    y = check_y(k)
    Li = L_inf(k)
    print(f"k={k}: y = {[str(v) for v in y]}")
    print(f"   L_inf computed = {Li} = {float(Li):.6f}")
    print(f"   L_inf claimed  = {claimed[k]} = {float(claimed[k]):.6f}  match={Li==claimed[k]}")
    print(f"   cap={caps[k]}  margin={caps[k]-float(Li):.4f}")
    print()
