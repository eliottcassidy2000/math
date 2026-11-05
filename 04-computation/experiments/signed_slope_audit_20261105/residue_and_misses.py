#!/usr/bin/env python3
"""(a) Empirical check of the STRONG_SOURCE_INTERFACE statement for the moduli used by the both-odd recipes:
    for every even M in a list and every pair of odd residues (D0, U0) mod M there are many odd n < 2^22 with
    (g(n), n) = (D0, U0) mod M (the orbit points n odd are the strong one-apex sources; D = g(n), U = n).
(b) The smallest-|target| misses of the small-witness probe (n < 2^24) for the class (i)/(ii) slopes, to judge the
    literal pre-registered predictions P1/P2 ("except possibly a few of the largest")."""
import numpy as np
NB = 24; top = 1 << NB
g = np.zeros(top + 1, dtype=np.int64); g[1] = 1
j = 1
while (1 << j) <= top:
    lo = 1 << j; hi = min(1 << (j + 1), top + 1)
    ev = np.arange(lo, hi, 2, dtype=np.int64); g[ev] = g[ev // 2]
    od = np.arange(lo + 1, hi, 2, dtype=np.int64); mm = (od + 1) // 2; g[od] = g[mm] + mm - (mm % 2)
    j += 1
nodd = np.arange(3, 1 << 22, 2, dtype=np.int64); godd = g[3:1 << 22:2]
print("(a) residue pairs (g(n) mod M, n mod M), odd n < 2^22, M even: missing odd pairs and the minimum count per pair")
for M in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 30, 36, 50]:
    key = (godd % M) * M + (nodd % M)
    counts = np.bincount(key, minlength=M * M)
    missing = [(d0, u0) for d0 in range(1, M, 2) for u0 in range(1, M, 2) if counts[d0 * M + u0] == 0]
    mincount = min(counts[d0 * M + u0] for d0 in range(1, M, 2) for u0 in range(1, M, 2))
    print(f"  M={M}: odd pairs {((M//2)**2)}, missing {len(missing)} {missing[:5]}, minimum multiplicity {mincount}")
print("(b) smallest-|t| misses in [-999, 999] of the right parity with n < 2^24")
nodd = np.arange(3, top, 2, dtype=np.int64); godd = g[3:top:2]
for (a, b) in [(2,3), (4,5), (6,7), (8,9), (4,7), (6,11), (8,11), (8,13), (10,11), (10,13), (3,5), (5,7), (5,9), (7,9), (7,11), (9,11), (7,13), (9,13)]:
    F = a * nodd - b * godd
    att = set(F[(F >= -999) & (F <= 999)].tolist())
    par = (a - b) % 2
    missed = sorted([t for t in range(-999, 1000) if t % 2 == par and t not in att], key=abs)
    print(f"  {a}/{b}: {len(missed)} misses; smallest |t|: {missed[:10]}")
