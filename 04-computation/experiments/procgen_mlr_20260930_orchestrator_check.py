#!/usr/bin/env python3
"""Orchestrator audit of lane `mlr` (multiplicative lonely runner), written from the note's statements; the
lane's code was not read.

  1. Proposition K: kappa of the 3-smooth boxes B(J,K) = {2^j 3^k : j < J, k < K}, exactly (max over critical
     times m/(v+v'), m/|v-v'|, (2m+1)/2v), for J <= 5, K <= 3: 1/2 (J=1), 1/3 (K=1, J>=2), 1/4 (J=2, K>=2),
     1/5 (J>=3, K>=2).
  2. Theorem S census: over all rationals r/N in (0,1), N <= 700, the lonely value I(r/N) = min over the
     x2,x3-orbit (residues reachable from r by multiplying by 2 and 3 mod N) of the centred distance; the
     points with I >= 1/14 number exactly 88, with denominators {5,7,10,11,13,14,26,28,33,52,56}, and the
     values >= 1/14 are exactly {1/5, 1/7, 1/10, 1/11, 1/13, 1/14}.
  3. The discrete LRC exceptions for the box (3,2) = {1,2,3,4,6,12}: among D coprime to 6, 5 <= D <= 3000,
     max_h min_v ||h v/D|| < 1/7 exactly for D in {11, 17, 37}.
"""
from fractions import Fraction as F
from math import gcd


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def dist(t, v):
    x = (t * v) % 1
    return min(x, 1 - x)


def kappa(V):
    cands = set()
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            for d in (V[i] + V[j], abs(V[i] - V[j])):
                for m in range(1, d):
                    cands.add(F(m, d))
        for m in range(V[i]):
            cands.add(F(2 * m + 1, 2 * V[i]))
    return max(min(dist(t, v) for v in V) for t in cands)


tab = {}
for J in range(1, 6):
    for K in range(1, 4):
        V = sorted(2 ** j * 3 ** k for j in range(J) for k in range(K))
        if len(V) < 2:
            continue
        tab[(J, K)] = kappa(V)
exp = lambda J, K: F(1, 2) if J == 1 else (F(1, 3) if K == 1 else (F(1, 4) if J == 2 else F(1, 5)))
check(all(v == exp(J, K) for (J, K), v in tab.items()), "Proposition K: kappa(B(J,K)) = 1/2, 1/3, 1/4, 1/5 by the stated cases, for all J <= 5, K <= 3 (exact)")


def I_value(r, N):
    seen = {r % N}
    stack = [r % N]
    while stack:
        x = stack.pop()
        for y in ((2 * x) % N, (3 * x) % N):
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return min(F(min(y, N - y), N) for y in seen)


pts = set()
vals = set()
for N in range(2, 701):
    for r in range(1, N):
        if gcd(r, N) != 1:
            continue
        v = I_value(r, N)
        if v >= F(1, 14):
            pts.add(F(r, N))
            vals.add(v)
dens = sorted(set(p.denominator for p in pts))
check(len(pts) == 88 and dens == [5, 7, 10, 11, 13, 14, 26, 28, 33, 52, 56] and vals == {F(1, 5), F(1, 7), F(1, 10), F(1, 11), F(1, 13), F(1, 14)},
      f"Theorem S census (all r/N, N <= 700): {len(pts)} points with I >= 1/14, denominators {dens}, values {sorted(vals)}")

box = [1, 2, 3, 4, 6, 12]
exc = []
for D in range(5, 3001):
    if gcd(D, 6) != 1:
        continue
    best = max(min(min((h * v) % D, D - (h * v) % D) for v in box) for h in range(1, D))
    if 7 * best < D:          # best/D < 1/7
        exc.append(D)
check(exc == [11, 17, 37], f"discrete LRC for the box (3,2): exceptions among D coprime to 6, 5 <= D <= 3000: {exc}")
