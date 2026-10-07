#!/usr/bin/env python3
"""Audit A, item 6(c): absorption by step K of the chain from (0,1) [n vs n+1] and from (1,0) [n vs 3n] depends only
on n mod 2^K and coincides with an actual equal-time integer merge; unmerged residue fractions."""
import random
from a3_drift import step
def T(x): return x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1
rng = random.Random(4581)
for start, name, f in (((0, 1), 'n vs n+1', lambda n: n + 1), ((1, 0), 'n vs 3n', lambda n: 3 * n)):
    for K in (8, 12, 16):
        unmerged = 0; dis = 0
        for r in range(2 ** K):
            res = []
            for M in (rng.getrandbits(40) + 1, rng.getrandbits(40) + 1):
                n = r + (M << K)
                k, N = start; u, v = f(n), n; ab = None
                for t in range(K + 1):
                    if (k, N) == (0, 0): ab = t; break
                    if t == K: break
                    k, N = step(k, N, v & 1); u, v = T(u), T(v)
                # actual integer merge by time K
                a, b = f(n), n; im = None
                for t in range(K + 1):
                    if a == b: im = t; break
                    a, b = T(a), T(b)
                if (ab is None) != (im is None) or (ab is not None and ab != im): dis += 1
                res.append(ab is None)
            if res[0] != res[1]: dis += 1
            unmerged += res[0]
        print(f"[{name}] K={K:2d}: unmerged residues {unmerged}/{2**K} = {unmerged/2**K:.5f}; disagreements {dis}")
