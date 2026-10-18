#!/usr/bin/env python3
"""Check the claimed C_2 coefficients of W_3(3,q) (q = 7, 8, 10, 11, 13) at two primes, and CRT-reconstruct the q = 13
x^6 coefficient."""
from fractions import Fraction
from aud_algebras import W3, LW, WW
from vamod import rowspace_basis
import math

def li0_rel(p, q, w, P):
    c = 2 * (1 - Fraction(12 * (p - q) ** 2, p * q))
    va = W3(c, P)
    G = va.gram_all(2 * w)
    I = va.radical(G)
    B = va.basis(2 * w)
    order = sorted(B, key=lambda m: (va.mono_li(m), m))
    basisI = rowspace_basis(I[(2 * w, 0)], order, P)
    out = []
    for v in basisI:
        sym = {m: cc for m, cc in v.items() if va.mono_li(m) == 0}
        if not sym: continue
        terms = {}
        for m, cc in sym.items():
            nx = sum(1 for X in m if X[0] == LW); ny = sum(1 for X in m if X[0] == WW)
            terms[(nx, ny)] = cc
        out.append(terms)
    return out

claims = {  # (q, weight): (leading monomial, {monomial: claimed coefficient})
    (7, 5): ((1, 1), {}),
    (7, 6): ((0, 2), {(3, 0): Fraction(-7, 39)}),
    (8, 6): ((0, 2), {(3, 0): Fraction(-16, 279)}),
    (8, 7): ((2, 1), {}),
    (10, 8): ((1, 2), {(4, 0): Fraction(-2, 123)}),
    (10, 9): ((0, 3), {(3, 1): Fraction(-8, 41)}),
    (11, 9): ((0, 3), {(3, 1): Fraction(-11, 138)}),
    (11, 10): ((2, 2), {(5, 0): Fraction(-11, 1035)}),
    (13, 11): ((1, 3), {(4, 1): Fraction(-13, 448)}),
    (13, 12): ((0, 4), {(3, 2): Fraction(-13, 56), (6, 0): Fraction(-28529, 3278)}),
}
vals = {}
for P in (2147483647, 2147483629):
    for (q, w), (lead, rest) in claims.items():
        rels = li0_rel(3, q, w, P)
        assert len(rels) == 1, (q, w, len(rels))
        t = rels[0]
        inv = pow(t[lead], P - 2, P)
        tn = {k: v * inv % P for k, v in t.items()}
        ok = all(tn.get(k, 0) == (Fraction(v).numerator * pow(Fraction(v).denominator, P - 2, P)) % P for k, v in rest.items())
        ok = ok and set(tn) == {lead} | set(rest)
        print(f'P={P} q={q} weight {w}: claimed coefficients', 'MATCH' if ok else f'MISMATCH  computed={tn}')
        vals.setdefault((q, w), []).append((P, tn))
# CRT + rational reconstruction for q=13 weight 12 x^6 coefficient
(P1, t1), (P2, t2) = vals[(13, 12)]
for mono in ((6, 0), (3, 2)):
    a1, a2 = t1[mono], t2[mono]
    M = P1 * P2
    x = (a1 + P1 * ((a2 - a1) * pow(P1, -1, P2) % P2)) % M
    # rational reconstruction
    r0, r1, s0, s1 = M, x, 0, 1
    B = math.isqrt(M // 2)
    while r1 > B:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    fr = Fraction(r1, s1)
    print('q=13 weight 12, coefficient of x^%d y^%d (y^4-normalised) by CRT of two primes:' % mono, fr)
