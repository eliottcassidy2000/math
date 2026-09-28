#!/usr/bin/env python3
"""Extras for the STICKY note: unconditional aliquot drift (all n, even n), the
free-cofactor proposition of the parallel session's tenth note (independent check),
and the class-entry rate P(v2 s(n) >= 2 | v2 n = 1) by size.
Run: python 04-computation/experiments/collatz_sticky_20260927_extra.py
"""
import math
import random
from collections import Counter

import numpy as np
from sympy import factorint


def v2(x):
    return (x & -x).bit_length() - 1 if x else 10 ** 9


def sigma_sieve(N):
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        s[d::d] += d
    return s


N = 2 * 10 ** 6
sig = sigma_sieve(N)
all_d, even_d, odd_d = [], [], []
for n in range(2, N + 1):
    s = int(sig[n]) - n
    if s <= 0:
        continue
    x = math.log2(s / n)
    all_d.append(x)
    (even_d if n % 2 == 0 else odd_d).append(x)
print(f"unconditional aliquot drift E[log2(s(n)/n)] over n <= {N}: all {sum(all_d)/len(all_d):+.3f}, even {sum(even_d)/len(even_d):+.3f}, odd {sum(odd_d)/len(odd_d):+.3f}")
print(f"P(s(n) > n) over n <= {N}: {sum(1 for n in range(2, N+1) if int(sig[n]) > 2*n)/N:.4f}")

# entry into growth classes from a = 1 by size: P(v2 sigma(m) = 1 | m odd) in [N/2, N)
print("entry rate P(v2 s(n) >= 2 | v2 n = 1) = P(v2 sigma(m) = 1) for odd m in [M/2, M):")
for M in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
    c = t = 0
    for m in range(M // 2 + 1, M, 2):
        t += 1
        c += (v2(int(sig[m])) == 1)
    print(f"   M={M:>8}: {c/t:.4f}  x ln M = {c/t*math.log(M):.3f}")
for M in (10 ** 7, 10 ** 8, 10 ** 9):
    rng = random.Random(M)
    c = t = 0
    for _ in range(20000):
        m = rng.randrange(M // 2, M) | 1
        tot = 0
        for p, e in factorint(m).items():
            if e % 2 == 1:
                tot += v2(p + 1) + v2(e + 1) - 1
        t += 1
        c += (tot == 1)
    print(f"   M={M:>8}: {c/t:.4f}  x ln M = {c/t*math.log(M):.3f} (sampled)")

# free-cofactor proposition (parallel tenth note, Proposition 3): m = 2^(K+1) t - 1, t odd
ok = 0
for K in range(1, 13):
    for t in range(1, 100, 2):
        m = 2 ** (K + 1) * t - 1
        x = m
        for _ in range(K):
            y = 3 * x + 1
            assert v2(y) == 1
            x = y // 2
        assert x == 2 * 3 ** K * t - 1
        y = 3 * x + 1
        assert v2(y) == 1 + v2(3 ** (K + 1) * t - 1)
        ok += 1
print(f"free cofactor: run of exactly K ones, U^K(2^(K+1)t-1) = 2*3^K t - 1, next valuation 1 + v2(3^(K+1) t - 1): {ok} cases OK")
lte = all(v2(3 ** (K + 1) - 1) == (1 if (K + 1) % 2 == 1 else 2 + v2(K + 1)) for K in range(1, 400))
print(f"LTE: v2(3^(K+1) - 1) = 1 (K+1 odd) or 2 + v2(K+1) (K+1 even) for K < 400: {lte}")
# Terras exactness: P(v2 = 1 | v1 = 1) exactly 1/2: n = 7 mod 8 within n = 3 mod 4
print("Terras exactness: (v1, v2) of odd n is a function of n mod 8: " +
      str(all(v2(3 * ((3 * n + 1) >> v2(3 * n + 1)) + 1) == 1 for n in range(7, 10 ** 5, 8)) and
          all(v2(3 * ((3 * n + 1) >> v2(3 * n + 1)) + 1) != 1 for n in range(3, 10 ** 5, 8))))
print("DONE")
