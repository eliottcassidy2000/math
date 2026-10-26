#!/usr/bin/env python3
"""Referee: A11 countdown = Hydra map x -> floor(3x/2) on x = a0 + b + 18 from the restart seeds
(b = 0 with a0 in {0,1}: x0 = 18, 19;  C(6,0,3,3): x0 = 27;  C(6,3,k,.): x0 = 30).
At a d = 3 end after K rounds: a0 = x_{K-1} mod 2, b = x_K - 18 - a0; HALT iff a0 = 0 and b parity = 1 - (c mod 2).
Measure, along each seed orbit, P(x_{K-1} even), the joint law of (x_{K-1}, x_K) parities, and the halting
frequency for c0 = 0, 1 (c0 = c mod 2)."""
import sys, time
Kmax = int(float(sys.argv[1])) if len(sys.argv) > 1 else 300000
for x0 in (18, 19, 27, 30):
    t0 = time.time()
    x = x0
    prev = x & 1
    cnt = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for K in range(1, Kmax + 1):
        x = x + (x >> 1)
        cur = x & 1
        cnt[(prev, cur)] += 1
        prev = cur
    n = Kmax
    pe = (cnt[(0, 0)] + cnt[(0, 1)]) / n
    h0 = cnt[(0, 1)] / n   # c0 = 0: need a0 = 0 (x_{K-1} even) and b odd (x_K odd)
    h1 = cnt[(0, 0)] / n   # c0 = 1: need x_{K-1} even and x_K even
    print(f'seed {x0}: K <= {n}: P(x_(K-1) even) = {pe:.4f}; joint (prev,cur): ' +
          ', '.join(f'{k}:{v / n:.4f}' for k, v in cnt.items()) +
          f'; halting freq c0=0: {h0:.4f}, c0=1: {h1:.4f}, mean {(h0 + h1) / 2:.4f} ({time.time() - t0:.0f}s)', flush=True)
