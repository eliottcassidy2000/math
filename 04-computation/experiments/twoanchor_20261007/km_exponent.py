#!/usr/bin/env python3
"""Karlin-McGregor test for Collatz translation partners (mac-mini-2026-10-07-twoanchor).
For Haar-random y (4096-bit random integers; parities exact 2-adically for T << 4096 steps), estimate
 q1(T) = P(y, y+1 not merged by Terras time T)            (THM-4593: ~ T^(-1/2))
 q2(T) = P(y, y+1, y+2 pairwise unmerged by time T)       (vicious-walker prediction: ~ T^(-3/2))
 q3(T) = P(y..y+3 pairwise unmerged by time T)            (prediction: ~ T^(-3))
Merges are equal-time coincidences of orbit values."""
import random, sys, math
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
TMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
checkpoints = [16, 32, 64, 128, 256, 512, 1024, 2048]
rnd = random.Random(99)
surv = {R: {c: 0 for c in checkpoints} for R in (1, 2, 3)}
for it in range(N):
    y = rnd.getrandbits(4096) | (1 << 4095)
    xs = [y + r for r in range(4)]
    alive = {1: True, 2: True, 3: True}
    for s in range(1, TMAX + 1):
        xs = [T(x) for x in xs]
        # pairwise distinctness among first R+1
        if alive[3] and len(set(xs)) < 4: alive[3] = False
        if alive[2] and len(set(xs[:3])) < 3: alive[2] = False
        if alive[1] and xs[0] == xs[1]: alive[1] = False
        if s in surv[1]:
            for R in (1, 2, 3):
                if alive[R]: surv[R][s] += 1
        if not any(alive.values()): break
for R in (1, 2, 3):
    row = [(c, surv[R][c] / N) for c in checkpoints]
    sl = []
    for (c1, p1), (c2, p2) in zip(row, row[1:]):
        if p1 > 0 and p2 > 0: sl.append(round(-math.log(p2/p1)/math.log(2), 2))
    print(f"R={R} survival", [(c, round(p, 5)) for c, p in row], "local slopes", sl, flush=True)
