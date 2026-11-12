#!/usr/bin/env python3
"""Value sets of affine squares (two disjoint arcs) whose four members are strong, n = 5, 6 (exhaustive, exact).
Used for the small cases of Theorem A (b): does a strong affine cube contain all odd numbers of some [X, 2X]?"""
import sys
import numpy as np
from collections import Counter
from arc_cube import table
for n in (5, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    H, S = table(n, pi, pj)
    na = len(pairs); N = 1 << na
    cnt = Counter()
    m = np.arange(N)
    for e in range(na):
        for f in range(e + 1, na):
            if set(pairs[e]) & set(pairs[f]): continue
            be, bf = 1 << e, 1 << f
            base = m[(m & (be | bf)) == 0]
            h0 = H[base]; h1 = H[base | be]; h2 = H[base | bf]; h3 = H[base | be | bf]
            ok = (h0 - h1 - h2 + h3 == 0) & (S[base] == 1) & (S[base | be] == 1) & (S[base | bf] == 1) & (S[base | be | bf] == 1)
            for a, b, c, d in zip(h0[ok], h1[ok], h2[ok], h3[ok]):
                cnt[tuple(sorted({int(a), int(b), int(c), int(d)}))] += 1
    print(f"n = {n}: strong affine squares by value set (set: count), sets with 4 distinct values listed first")
    four = sorted((k, v) for k, v in cnt.items() if len(k) == 4)
    print("   4 distinct values:", four[:60], " total kinds:", len(four))
    consec = [(k, v) for k, v in four if k[3] - k[0] == 6]
    print("   4 consecutive odd values:", consec)
    cover = [(k, v) for k, v in consec if k[0] <= 2 * ((k[0] - 1)) and k[3] >= 2 * (k[0] - 1) - 1]
    print("   of these, sets containing all odd numbers of [X, 2X] with X = min - 1:", cover)
    print("   3 distinct values:", sorted((k, v) for k, v in cnt.items() if len(k) == 3)[:40])
