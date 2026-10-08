#!/usr/bin/env python3
"""HYP-9241 test at R = 3: does P(y, y+1, y+2, y+3 pairwise unmerged by T) equal the product of the six pair survival
probabilities P(y+i, y+j unmerged by T) (independent pair events)? Haar-like y: random 2048-bit integers."""
import random, sys, itertools, math
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
TMAX = 1024; cps = [32, 64, 128, 256, 512, 1024]
pairs = list(itertools.combinations(range(4), 2))
rnd = random.Random(20261008)
alive_pair = {p: {c: 0 for c in cps} for p in pairs}
alive_all = {c: 0 for c in cps}; alive_3 = {c: 0 for c in cps}
for it in range(N):
    y = rnd.getrandbits(2048) | (1 << 2047)
    xs = [y + r for r in range(4)]
    dead = set()
    for s in range(1, TMAX + 1):
        xs = [T(x) for x in xs]
        for p in pairs:
            if p not in dead and xs[p[0]] == xs[p[1]]: dead.add(p)
        if s in alive_all:
            for p in pairs:
                if p not in dead: alive_pair[p][s] += 1
            if not dead: alive_all[s] += 1
            if not any(p in dead for p in pairs if p[1] <= 2): alive_3[s] += 1
        if len(dead) == 6: break
print(f"N = {N}")
for c in cps:
    prod6 = 1.0
    for p in pairs: prod6 *= alive_pair[p][c] / N
    prod3 = (alive_pair[(0,1)][c]/N) * (alive_pair[(1,2)][c]/N) * (alive_pair[(0,2)][c]/N)
    q3, q2 = alive_all[c] / N, alive_3[c] / N
    print(f"T={c:5d}: q2={q2:.4f} prod3={prod3:.4f} ratio={q2/prod3 if prod3 else float('nan'):.3f} | q3={q3:.5f} prod6={prod6:.5f} ratio={q3/prod6 if prod6 else float('nan'):.3f} (events {alive_all[c]})")
