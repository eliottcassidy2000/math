#!/usr/bin/env python3
"""false_peaks_cd_tower_kps_S128c22.py -- kind-pasteur S128 cont.22.
THE CLIMB ASYMMETRY: down-moves (dc3 = -1) always exist off the transitive (THM-854); do UP-moves
(dc3 = +1: an arc u->v with d_u = d_v + 2) always exist below c3max? Census of FALSE PEAKS (local
maxima of c3 under single flips that are not global) n = 4, 5, 6 -- the Morse obstruction of the
+8 climb. Plus the CD-tower level-count table."""
import sys
from math import comb
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)
for n in [4, 5, 6]:
    m = comb(n, 2)
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    c3max = (n**3 - n) // 24 if n % 2 == 1 else (n**3 - 4*n) // 24
    peaks = Counter(); tot_below = 0; no_up = 0
    for mask in range(1 << m):
        B = [[False] * n for _ in range(n)]
        for k, (u, v) in enumerate(pairs):
            if (mask >> k) & 1: B[u][v] = True
            else: B[v][u] = True
        d = [sum(r) for r in B]
        c3 = comb(n, 3) - sum(comb(x, 2) for x in d)
        if c3 == c3max: continue
        tot_below += 1
        up = any(B[u][v] and d[u] == d[v] + 2 for u in range(n) for v in range(n) if u != v)
        if not up:
            no_up += 1
            peaks[c3] += 1
    print("n=%d: below-max tournaments %d ; FALSE PEAKS (no +1 move) %d ; by c3 %s ; c3max=%d"
          % (n, tot_below, no_up, dict(peaks), c3max))
print()
print("CD-tower level counts (populated E4 levels = c3max+1): n=2,3,5,9,17 ->",
      [ (n, ((n**3-n)//24 if n%2 else (n**3-4*n)//24) + 1) for n in [3,5,9,17] ])
print("Moser rank note: A000127 = sum_{k<=4} C(n-1,k) -- truncation order 4 = last solvable S_k")
print("DONE")
