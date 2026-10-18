#!/usr/bin/env python3
"""Schur index of 4d SU(2) with N_f = 4 (8 half-hypers in (2, 8_v) of SU(2) x SO(8)):
  total (unflavoured) and the SO(8) zero-weight part, via the factorisation over the four fugacities z_i:
  hyper PE = prod_i F(z_i),  F(z) = prod_{n>=0, eps, delta = +-1} 1/(1 - s^(2n+1) a^eps z^delta),  s = q^(1/2);
  zero-weight part uses G(s,a) = [z^0] F(z) to the 4th power.  Haar (1/2)(1-a^2)(1-a^-2); vector PE
  prod_{n>=1} (1-q^n)^2 (1-q^n a^2)^2 (1-q^n a^-2)^2."""
import sys
from collections import defaultdict
Dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
S = 2 * Dmax  # max power of s

def mul(A, B):
    out = defaultdict(int)
    for (s1, a1, z1), c1 in A.items():
        for (s2, a2, z2), c2 in B.items():
            if s1 + s2 <= S:
                out[(s1 + s2, a1 + a2, z1 + z2)] += c1 * c2
    return {k: v for k, v in out.items() if v}

def geom(spow, apow, zpow):
    out = {}
    kk = 0
    while kk * spow <= S:
        out[(kk * spow, kk * apow, kk * zpow)] = 1
        kk += 1
    return out

# F(z) with z tracked
F = {(0, 0, 0): 1}
for n in range(0, S):
    if 2 * n + 1 > S: break
    for eps in (1, -1):
        for dl in (1, -1):
            F = mul(F, geom(2 * n + 1, eps, dl))
G = defaultdict(int)       # [z^0] F
Ftot = defaultdict(int)    # F at z = 1
for (sp_, a, z), c in F.items():
    Ftot[(sp_, a, 0)] += c
    if z == 0:
        G[(sp_, a, 0)] += c
G = dict(G); Ftot = dict(Ftot)
def power(X, m):
    R = {(0, 0, 0): 1}
    for _ in range(m):
        R = mul(R, X)
    return R
Vec = {(0, 0, 0): 1}
for n in range(1, Dmax + 1):
    for ap in (0, 2, -2):
        for _ in range(2):
            fac = {(0, 0, 0): 1, (2 * n, ap, 0): -1}
            Vec = mul(Vec, fac)
haar = {(0, 0, 0): 1, (0, 2, 0): -0.5, (0, -2, 0): -0.5}
def schur(H4):
    tot = mul(mul(Vec, H4), {(0, 0, 0): 2, (0, 2, 0): -1, (0, -2, 0): -1})
    res = defaultdict(int)
    for (sp_, a, z), c in tot.items():
        if a == 0:
            res[sp_] += c
    return [res.get(2 * D, 0) // 2 for D in range(Dmax + 1)], [res.get(2 * D + 1, 0) for D in range(Dmax)]
tot, odd = schur(power(Ftot, 4))
zero, odd0 = schur(power(G, 4))
print("SU(2) Nf=4 Schur index (unflavoured):", tot, " odd powers:", set(odd))
print("SO(8) zero-weight part:              ", zero, " odd powers:", set(odd0))
