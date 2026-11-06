#!/usr/bin/env python3
"""flip_stats.py (klein-2026-11-06): statistics of the Hamiltonian-path gradient along the flip edges of the tournament
metagraph, from the hflip outputs flips_<n>.txt (lines "F h c3 strong d_1 ... d_m" and "G e_1 ... e_m", d = change of h,
e = change of c3 for the reversal of each arc).
usage: flip_stats.py dir nmin nmax"""
import sys
from collections import defaultdict
D = sys.argv[1]; NMIN = int(sys.argv[2]); NMAX = int(sys.argv[3])
for n in range(NMIN, NMAX + 1):
    tot = 0; same = defaultdict(int); up_up = up_tot = down_down = down_tot = lev_up = lev_tot = 0
    dc3 = defaultdict(int); maxratio = 0.0; zero_by_level = defaultdict(lambda: [0, 0]); dh_values = defaultdict(set)
    maxrel = defaultdict(float)
    with open(f"{D}/flips_{n}.txt") as fh:
        F = None
        for line in fh:
            a = line.split()
            if not a: continue
            if a[0] == 'F': F = a; continue
            if a[0] != 'G' or F is None: continue
            h = int(F[1]); k = int(F[2]); ds = list(map(int, F[4:])); es = list(map(int, a[1:]))
            for d, e in zip(ds, es):
                tot += 1; dc3[e] += 1
                zero_by_level[k][1] += 1
                if d == 0: zero_by_level[k][0] += 1
                if e > 0:
                    up_tot += 1; up_up += (d > 0)
                elif e < 0:
                    down_tot += 1; down_down += (d < 0)
                else:
                    lev_tot += 1; lev_up += (d > 0)
                r = abs(d) / h
                if r > maxratio: maxratio = r
                if r > maxrel[k]: maxrel[k] = r
            F = None
    print(f"order {n}: {tot} flips; level change e = c3' - c3 distribution: {dict(sorted(dc3.items()))}")
    print(f"  level-raising flips (e > 0): {up_tot}, of which h rises: {up_up} ({up_up / up_tot:.4f}); level-lowering flips: {down_tot}, h falls: {down_down} ({down_down / down_tot:.4f}); level-preserving: {lev_tot}, h rises: {lev_up}")
    print(f"  max |dh| / h over all flips: {maxratio:.3f}; by level: " + " ".join(f"{k}:{maxrel[k]:.2f}" for k in sorted(maxrel)))
    print("  fraction of h-preserving flips by level: " + " ".join(f"{k}:{z[0] / z[1]:.3f}" for k, z in sorted(zero_by_level.items())))
