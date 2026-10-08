#!/usr/bin/env python3
"""Rank-1 two-valued translation-only maps (multipliers 1 and mu, mu = 1 mod d): does y, y+1 coalesce a.s. with a
T^(-1/2) tail?  Regimes: mu <= d^2 (THM-4581-type level weight exists) and mu > d^2 (no balancing weight).
Exact pair chain from (M, e) = (1, 1), fresh uniform digits.  Usage: NSAMP TMAX"""
import random, math, sys
from fractions import Fraction as Fr
MAPS = {'Z3 (1,1,4) mu=4<9': (3, [1, 1, 4]), 'Z5 (1,1,1,1,26) mu=26>25': (5, [1, 1, 1, 1, 26]),
        'Z5 (1,1,1,1,6) mu=6<25': (5, [1, 1, 1, 1, 6]), 'Z3 (1,4,4) mu=4': (3, [1, 4, 4])}
def modd(q, d): return (q.numerator * pow(q.denominator, -1, d)) % d
nsamp, tmax = int(sys.argv[1]), int(sys.argv[2])
rnd = random.Random(8)
for name, (d, m) in MAPS.items():
    r = [(-m[i] * i) % d for i in range(d)]
    Lam = sum(math.log(x / d) for x in m) / d
    times = []
    for _ in range(nsamp):
        M, e = Fr(1), Fr(1); t = 0
        while t < tmax and not (M == 1 and e == 0):
            j = rnd.randrange(d); i = (modd(M, d) * j + modd(e, d)) % d
            e = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / d; M = M * Fr(m[i], m[j]); t += 1
        times.append(t if (M == 1 and e == 0) else None)
    row = []
    T = 16
    while T <= tmax:
        q = sum(1 for x in times if x is None or x > T) / nsamp; row.append(f"T={T}: q={q:.4f} sqrtT*q={math.sqrt(T)*q:.2f}"); T *= 4
    print(f"{name}: Lambda={Lam:+.3f}  " + "  ".join(row), flush=True)
