#!/usr/bin/env python3
"""Recon B for S49: the B5 limit constant on the canonical family v = (1..13).

B5 = sum_{d<=5} (-1)^d S_d.  Compute B5/(q-1) at several large coprime q and
extract the limit constant; positivity => the a-priori live floor on the
canonical family is one explicit-threshold statement away (the nucleus demo).
Also S2/(q-1) vs the assembled pair law (49 locked closed forms + 29 sparse
branch laws) as an independent cross-check.
"""
from math import gcd, comb
from fractions import Fraction as F

v = list(range(1, 14))

def band_fail(vi, q, p):
    r = (vi * p) % q
    return 14 * r < q or 14 * r > 13 * q

rows = []
for q in [2003, 5003, 10007, 20011]:   # primes > 13 => coprime to everything
    S = [0] * 6
    live = 0
    for p in range(1, q):
        c = sum(1 for vi in v if band_fail(vi, q, p))
        for d in range(6):
            S[d] += comb(c, d)
        live += (c == 0)
    B5 = sum((-1) ** d * S[d] for d in range(6))
    rows.append((q, F(B5, q - 1), F(live, q - 1), F(S[2], q - 1)))
    print(f"q={q}: B5/(q-1) = {B5}/{q-1} = {B5/(q-1):+.6f}   live/(q-1) = {live/(q-1):.6f}"
          f"   S2/(q-1) = {S[2]/(q-1):.6f}")

# assembled pair law for S2 limit
tot = F(0)
for a in range(1, 14):
    for b in range(a + 1, 14):
        g = gcd(a, b); ip, jp = a // g, b // g
        M = max(ip, jp)
        if ip + jp <= 13:
            tot += F(1, 7 * M)
        else:
            tot += F(1, 7 * M) + 2 * F(ip + jp - 14, 14 * ip * jp)
print(f"assembled S2 limit = {tot} = {float(tot):.6f}   equilibrium 78/49 = {float(F(78,49)):.6f}")
print(f"S2 excess over equilibrium: {tot - F(78,49)} = {float(tot - F(78,49)):+.6f}")
