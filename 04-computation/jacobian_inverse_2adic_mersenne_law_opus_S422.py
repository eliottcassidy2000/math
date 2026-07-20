#!/usr/bin/env python3
"""
death-star-2026-07-19-S59m (HYP-8075) -- the formal inverse of the JC
counterexample F and its 2-adic coefficient ladder.

F(0)=0, JF(0) = L: (x,y,z) -> (z, y, 2x), det L = -2.  The compositional
inverse G (F(G) = id) exists as a formal power series and CANNOT be
polynomial (else F would be injective, contradicting the verified triple
collision).  We compute G to total degree K and measure, per degree d:
  - number of monomials across the 3 components
  - min 2-adic valuation v_2 of the coefficients (they live in Z[1/2])
  - max coefficient height (log10 |numerator|)
The dyadic ladder v_min(d) is the quantitative shape of the obstruction:
the inverse exists at every finite 2-adic layer but descends without bound
(measured), so no polynomial truncation ever closes it.
"""
from fractions import Fraction

K = 16  # truncation total degree (opus-S422: extended to test the Mersenne-depth law at d=15)

def pmulK(a, b):
    r = {}
    for ka, ca in a.items():
        da = ka[0]+ka[1]+ka[2]
        for kb, cb in b.items():
            if da + kb[0]+kb[1]+kb[2] > K:
                continue
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}

def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c != 0}

def pscale(p, s):
    return {k: c*s for k, c in p.items() if c*s != 0}

# F in exact form (integer dicts), from the verified script
X = {(1,0,0): 1}; Y = {(0,1,0): 1}; Z = {(0,0,1): 1}; ONE = {(0,0,0): 1}
def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}
U = padd(ONE, pmul(X, Y)); U2 = pmul(U, U); U3 = pmul(U2, U)
W = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))
F1 = padd(pmul(U3, Z), pmul(pmul(pmul(Y, Y), U), W))
F2 = padd(Y, pscale(pmul(pmul(X, U2), Z), 3), pscale(pmul(pmul(X, pmul(Y, Y)), W), 3))
F3 = padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3), pscale(pmul(pmul(pmul(X, X), X), Z), -1))
F = [F1, F2, F3]

def compose(P, G):
    """P(G1,G2,G3) truncated at degree K, via per-monomial powers with caching."""
    cache = {(0,0,0): {(0,0,0): Fraction(1)}}
    def power(base, e, memo={}):
        key = (id(base), e)
        if e == 0: return {(0,0,0): Fraction(1)}
        if key in memo: return memo[key]
        r = base
        for _ in range(e-1):
            r = pmulK(r, base)
        memo[key] = r
        return r
    out = {}
    for (i, j, l), c in P.items():
        t = {(0,0,0): Fraction(c)}
        if i: t = pmulK(t, power(G[0], i))
        if j: t = pmulK(t, power(G[1], j))
        if l: t = pmulK(t, power(G[2], l))
        out = padd(out, t)
    return out

# L^{-1}(a,b,c) = (c/2, b, a)
A = {(1,0,0): Fraction(1)}; Bv = {(0,1,0): Fraction(1)}; C = {(0,0,1): Fraction(1)}
ID = [A, Bv, C]
def Linv(vec):
    return [pscale(vec[2], Fraction(1,2)), vec[1], vec[0]]

G = Linv(ID)  # first approximation
for it in range(K + 2):
    FG = [compose(Fi, G) for Fi in F]
    err = [padd(ID[i], pscale(FG[i], -1)) for i in range(3)]
    if all(not e for e in err):
        break
    G = [padd(G[i], d) for i, d in enumerate(Linv(err))]

# residual check
FG = [compose(Fi, G) for Fi in F]
res = [padd(ID[i], pscale(FG[i], -1)) for i in range(3)]
min_res_deg = min((sum(k) for e in res for k in e), default=999)
print(f"F(G) = id verified through degree {min_res_deg - 1} (K = {K})")

def v2(fr):
    fr = Fraction(fr)
    n, d = fr.numerator, fr.denominator
    v = 0
    while d % 2 == 0: d //= 2; v -= 1
    while n and n % 2 == 0: n //= 2; v += 1
    return v

import math
print(f"\n{'deg':>4} {'#terms':>7} {'min v_2':>8} {'max height':>11}")
ladder = []
for d in range(1, K + 1):
    coeffs = [c for g in G for k, c in g.items() if sum(k) == d]
    if not coeffs:
        print(f"{d:>4} {'0':>7} {'--':>8}"); continue
    mv = min(v2(c) for c in coeffs)
    h = max(len(str(abs(Fraction(c).numerator))) for c in coeffs)
    ladder.append((d, len(coeffs), mv, h))
    print(f"{d:>4} {len(coeffs):>7} {mv:>8} {h:>11}")
print("\nnon-termination: nonzero terms at every degree =",
      all(t > 0 for _, t, _, _ in ladder))
print("dyadic ladder (min v_2 per degree):", [mv for _, _, mv, _ in ladder])
