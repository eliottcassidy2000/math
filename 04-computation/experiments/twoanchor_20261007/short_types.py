#!/usr/bin/env python3
"""Fixed end-of-run chain states for the short run types j = 3..6 (and j >= 7 for D = 3..6) (mac-mini-2026-10-07-twoanchor).
For each type and D, the state (k, e) of the D-chain at Terras time 2J after the run end x must be the same for all K, t."""
import random
from collections import defaultdict
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
rnd = random.Random(1)
table = defaultdict(set)
for j in (3, 4, 5, 6, 7, 8, 11, 12):
    J = (j - 1)//2
    for _ in range(150):
        K = rnd.randint(9, 40)
        M = 1 << (j - 1)
        t0 = pow(3, -(K-1), M)
        while True:
            t = t0 + M * rnd.getrandbits(80)
            if t & 1 and v2(3**(K-1)*t - 1) == j - 1: break
        x = 2*3**(K-1)*t - 1
        for D in range(1, 7):
            y = (x + 1)//3**D - 1
            u, v, k = x, y, D
            for s in range(2*J):
                k += (u & 1) - (v & 1); u, v = T(u), T(v)
            table[(j, D)].add((k, u - 3**k * v if k >= 0 else None))
for j in (3, 4, 5, 6, 7, 8, 11, 12):
    row = []
    for D in range(1, 7):
        st = table[(j, D)]
        row.append(str(next(iter(st))) if len(st) == 1 else f"{len(st)} states")
    print(f"j={j:2d} J={(j-1)//2}:", " | ".join(f"D={D}: {r}" for D, r in zip(range(1, 7), row)))
