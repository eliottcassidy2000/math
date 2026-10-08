#!/usr/bin/env python3
"""Depth spectrum of deletion certificates on the Mersenne residual line (odd K): the post-run Terras depth s at which
2^K - 1 merges at equal time with 2^(K-1) - 1 (D = 1) or 2^(K-3) - 1 (D = 3), s measured from the source's run end x
(Terras time K-1). Also the descent depth (first time the orbit drops below n) for comparison."""
import sys, math
K0, K1 = int(sys.argv[1]), int(sys.argv[2])
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
out = open(f'depth_spectrum_{K0}_{K1}.txt', 'w')
for K in range(K0, K1 + 1, 2):
    n = (1 << K) - 1
    x = 2 * 3**(K-1) - 1                     # source at Terras time K-1
    res = {}
    for D in (1, 3):
        y = 2 * 3**(K-1-D) - 1               # child run end, aligned (shift D)
        u, v, k = x, y, D
        s = 0; cap = 12 * K
        while s < cap:
            if u == v and k == 0: break
            if v == 1 and u == 1: s = None; break
            k += (u & 1) - (v & 1); u, v = T(u), T(v); s += 1
        res[D] = s if (s is not None and s < cap and u == v) else None
    # descent depth: first Terras time t > K-1 with orbit value < n
    u = x; s = 0
    while u >= n: u = T(u); s += 1
    out.write(f"{K} {res[1]} {res[3]} {s}\n")
    if K % 200 == 1: out.flush()
out.close()
