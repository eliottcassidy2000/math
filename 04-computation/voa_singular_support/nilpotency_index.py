#!/usr/bin/env python3
"""Nilpotency index of the stress tensor in Zhu's C_2-algebra R_V (Deb-Meneghelli-Rastelli 2025):
smallest n with T^n in grI (Li degree 0).  For L_{-1/2}(sl2) T = Sugawara ~ Cas = e f + h^2/4;
for the small N=4 algebra at c=-9 the Sugawara null identifies x_L with Cas, so test Cas^n as well as x_L^n."""
import sys
from fractions import Fraction
from lie_voa_ss import *
from small_n4 import SmallN4, L as NL, E as NE, H as NH, Fm as NF, G1, G2, B1, B2

def polymul(ss, a, b):
    out = {}
    for ma, ca in a.items():
        # multiply monomial-by-monomial using mult_var repeatedly (even variables only here)
        cur = dict((mb, cb * ca % ss.P) for mb, cb in b.items())
        for var in ma:
            cur = ss.mult_var(cur, var)
        for k, v in cur.items():
            out[k] = (out.get(k, 0) + v) % ss.P
    return {k: v for k, v in out.items() if v}

def member(ss, vec, D, p, ch):
    g = ss.grI.get((D, p, ch), [])
    r0 = ss._rank(g, D, p, ch)
    r1 = ss._rank(g + [vec], D, p, ch)
    return r1 == r0

P = 2147483647
inv4 = pow(4, -1, P)
# ---- L_{-1/2}(sl2)
alg = AffineSL2(Fraction(-1, 2), P)
sv = find_singular(alg, 4, (4,))
ss = SingularSupport(alg, 8, sv, verbose=False)
ss.run(gen_modes=[(1, 0), (1, 1), (1, 2)], extra_ops=[(0, 2), (0, 0)])
e0, h0, f0 = (1, 0), (1, 1), (1, 2)
Cas = {(f0, e0): 1, (h0, h0): inv4}
pw = {(): 1}
for n in range(1, 5):
    pw = polymul(ss, pw, Cas)
    D = 2 * n
    if D > 8: break
    print(f"L_-1/2(sl2): Cas^{n} in grI (weight {D}, Li 0)? {member(ss, pw, D, 0, (0,))}")
# ---- small N=4, c=-9
alg4 = SmallN4(-9, P)
sv4 = find_singular(alg4, 2, (0,))
ss4 = SingularSupport(alg4, 6, sv4, verbose=False)
h = Fraction(1, 2)
ss4.run(gen_modes=[(h, G1), (h, G2), (h, B1), (h, B2), (1, NE), (1, NH), (1, NF)], extra_ops=[(0, NF), (0, NE)])
xL = {((2, NL),): 1}
pw = {(): 1}
for n in range(1, 4):
    pw = polymul(ss4, pw, xL)
    D = 2 * n
    print(f"small N=4 c=-9: T^{n} (x_L^{n}) in grI (weight {D}, Li 0)? {member(ss4, pw, D, 0, (0,))}")
