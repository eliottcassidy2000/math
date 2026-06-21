#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_L7_tight_cluster_ratio_THREADA_opus.py   (THREAD A, opus)

L7 SOLE GAP -- TIGHT two-cluster, the REAL danger model.

The genuine L7 worry is NOT clusters at very different scales (those separate;
measS7 collapses, huge margin -- verified in lrc14_L7_genuine_2cluster).  The
worry is two TIGHT clusters at bases f1 < f2 with a GAP, where rho=f2/f1 is NEAR 1:
then E ~ one big consecutive run (the consec maximizer, ~0.327 at k=8) and we must
show measS7 stays <= cap as rho moves off 1.

MODEL:  cluster1 = {f1, f1+1, ..., f1+s-1}  (s tight consecutive ints at base f1)
        cluster2 = {f2, f2+1, ..., f2+s-1}  (s tight consecutive ints at base f2)
        E = cluster1 U cluster2,  k=2s,  rho = f2/f1 (center-ish ratio).
        GAP between clusters: f2 - (f1+s-1) >= 1 (a real gap; else they MERGE to consec).

To sweep rho continuously we hold s fixed and vary (f1, f2) with f2 > f1+s-1 (real gap).
Because measS7 is DILATION-INVARIANT, only the *shape* matters; large f1 gives fine
ratio resolution.  We sweep f1 up to F1MAX and report measS7 vs rho.

Confirms:
  (a) at rho->1+ (gap=1, smallest separation) measS7 approaches the consec value;
  (b) the actual WORST (max) measS7 over genuine-gap tight 2-clusters and its margin;
  (c) where measS7 drops well below cap as rho grows (the ~2.15 separation point).

EXACT Fractions.  cap_k: 8:2243/5880, 9:1979/4004, 10:55/91, 11:66/91, 12:6/7, 13:1.
"""
import sys
from math import gcd, comb
from fractions import Fraction as F
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(int(e) for e in E if e != 0))
    bps = set([F(0), F(1)])
    for e in E:
        ae = abs(e)
        for m in range(0, 7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = set(int(((e * xm) % 1) * 7) for e in E)
        if len(secs) == 7: total += x1 - x0
    return total

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
       11: F(66, 91), 12: F(6, 7), 13: F(1)}

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

def tight_pair(f1, f2, s):
    return tuple(sorted(set(range(f1, f1 + s)) | set(range(f2, f2 + s))))

def consec_meas(k):
    return measS7(list(range(1, k + 1)))

if __name__ == "__main__":
    F1MAX = 40
    for s in [4, 5]:
        k = 2 * s
        c = CAP[k]
        cm = consec_meas(k)
        print("=" * 78)
        print(f"  TIGHT 2-cluster  E={{f1..f1+{s-1}}} U {{f2..f2+{s-1}}}, k={k}")
        print(f"  cap_{k}={float(c):.5f} ({c}),  consec_{k} measS7={float(cm):.5f} ({cm})")
        print("=" * 78)
        rows = []
        for f1 in range(1, F1MAX + 1):
            # real gap: f2 >= f1+s  (gap = f2-(f1+s-1) >= 1). sweep rho up to ~3.2
            for f2 in range(f1 + s, int((f1 + s) * 3.3) + 2):
                E = tight_pair(f1, f2, s)
                if len(E) != k:    # overlap (shouldn't happen with gap>=1) - skip
                    continue
                if not primitive(E):
                    continue
                m = measS7(E)
                # define rho by cluster centers
                c1 = F(2 * f1 + s - 1, 2); c2 = F(2 * f2 + s - 1, 2)
                rho = c2 / c1
                gap = f2 - (f1 + s - 1)
                rows.append((m, rho, gap, f1, f2, E))
        # WORST
        bm, brho, bgap, bf1, bf2, bE = max(rows, key=lambda r: r[0])
        print(f"  scanned {len(rows)} genuine-gap shapes (f1<= {F1MAX}).")
        print(f"  >>> WORST tight 2-cluster: measS7={float(bm):.6f} ({bm})")
        print(f"      rho={float(brho):.4f}  gap={bgap}  f1={bf1},f2={bf2}  E={list(bE)}")
        print(f"      cap_{k}={float(c):.6f}  MARGIN={float(c-bm):+.6f} ({c-bm})  "
              f"{'BREACH' if bm>c else 'OK'}")
        print(f"      (consec measS7={float(cm):.6f}; worst-2cluster {'<' if bm<cm else '>='} consec)")
        # TOP 8
        rows.sort(reverse=True)
        print("  TOP 8 by measS7:")
        for (m, rho, gap, f1, f2, E) in rows[:8]:
            print(f"      measS7={float(m):.5f}  rho={float(rho):.3f} gap={gap} f1={f1},f2={f2}  E={list(E)}")
        # trend: max measS7 as a function of rho (bucket rho to 2 decimals, take max)
        from collections import defaultdict
        bucket = defaultdict(lambda: F(0))
        for (m, rho, gap, f1, f2, E) in rows:
            key = round(float(rho), 2)
            if m > bucket[key]: bucket[key] = m
        print("  [max measS7 per rho bucket (0.01 res), rho in (1,3.2)]:")
        last_over = None
        for key in sorted(bucket):
            m = bucket[key]
            tag = "OVER" if m > c else "ok"
            print(f"      rho~{key:.2f}: max measS7={float(m):.5f}  cap-m={float(c-m):+.5f} [{tag}]")
        print()
