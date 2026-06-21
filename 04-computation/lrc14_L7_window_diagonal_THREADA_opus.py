#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_L7_window_diagonal_THREADA_opus.py   (THREAD A, opus)

THE (1, 2.15) WINDOW, resolved.  The danger of a balanced 2-cluster is NOT a big ratio
(those separate -> measS7 collapses).  The danger is rho NEAR 1: then the line
x->(frac(f1 x),frac(f2 x)) hugs the DIAGONAL of the 2-torus, the two clusters read the
SAME torus coordinate, and the union behaves like ONE cluster of k offsets near a single
scale -- i.e. like consec_k.  As rho moves away from 1 (off-diagonal slope), the clusters
decorrelate and measS7 drops.

This script pins WHERE measS7 drops below the consec_k value as rho grows, by holding the
2-cluster at the SAME base scale (so it never merges into a literal AP) and sweeping the
ratio CONTINUOUSLY via large denominators.

CONCRETE MODEL (avoids the merged-AP degeneracy):
   cluster1 = {N, N+1, ..., N+s-1}                 (s tight ints at scale ~N)
   cluster2 = {round(rho*N), ..., round(rho*N)+s-1} (s tight ints at scale ~rho*N)
   with N LARGE so rho is resolved finely and the two clusters are SEPARATE integer sets
   (real internal gap, never an AP).
We sweep rho over a fine grid in (1, 3) and report max measS7 per rho, and the
consec_k reference.  The claim "above ~2.15 they separate" = measS7 stays low for all rho,
but the PEAK band (where measS7 is largest) sits near rho->1.

EXACT Fractions.
"""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict
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

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

def is_AP(E):
    E = sorted(E); d = E[1]-E[0]
    return all(E[i+1]-E[i] == d for i in range(len(E)-1))

if __name__ == "__main__":
    for s in [4, 5]:
        k = 2*s; c = CAP[k]
        cm = measS7(list(range(1, k+1)))
        print("="*78)
        print(f"  DIAGONAL/WINDOW probe, k={k}, cap={float(c):.5f}, consec={float(cm):.5f}")
        print(f"  cluster2 base = round(rho*N), N fixed; both clusters size {s}; rho in (1,3)")
        print("="*78)
        # for each rho = a/b (b<=12), use N a multiple of b so rho*N is exact integer
        results = []  # (rho_float, max measS7 over N, exampleE)
        ratios = sorted(set(F(a,b) for b in range(1,13) for a in range(b+1, 3*b+1) if gcd(a,b)==1))
        for rho in ratios:
            a, b = rho.numerator, rho.denominator
            best = F(0); bestE = None
            # N must make rho*N integer => N multiple of b. sweep a few N to vary phase/overlap.
            for nmul in range(1, 9):
                N = b * nmul
                base2 = (a * N) // b
                if base2 < N + s:        # clusters would overlap/merge -> require real gap
                    continue
                E = tuple(sorted(set(range(N, N+s)) | set(range(base2, base2+s))))
                if len(E) != k or not primitive(E):
                    continue
                if is_AP(E):             # exclude merged-AP (that's consec, not a 2-cluster)
                    continue
                m = measS7(E)
                if m > best: best, bestE = m, E
            if bestE is not None:
                results.append((float(rho), best, bestE))
        # report: peak band
        results.sort(key=lambda r: -float(r[1]))
        print("  TOP 12 (rho, max measS7 over phase) -- the PEAK band:")
        for (rf, m, E) in results[:12]:
            print(f"     rho={rf:.4f}  measS7={float(m):.5f} ({m})  margin={float(c-m):+.5f}  E={list(E)}")
        # trend sorted by rho
        results.sort(key=lambda r: r[0])
        print("  [max measS7 vs rho, all separated 2-clusters]:")
        prev_above_consec = True
        crossing = None
        for (rf, m, E) in results:
            above = m >= cm
            tag = ">=consec" if above else "<consec"
            print(f"     rho={rf:.3f}: measS7={float(m):.5f}  (cap-m={float(c-m):+.4f}) [{tag}]")
        # where does the separated max fall below, say, 0.9*consec? (the practical separation pt)
        thr = cm * F(9,10)
        sep_pt = None
        for (rf, m, E) in sorted(results, key=lambda r: r[0]):
            if m < thr and sep_pt is None and rf > 1.3:
                sep_pt = rf
        print(f"  consec={float(cm):.5f}; 0.9*consec={float(thr):.5f}; "
              f"first rho>1.3 with separated-max<0.9consec: {sep_pt}")
        # absolute worst separated and its margin
        bm = max(results, key=lambda r: r[1])
        print(f"  >>> WORST separated (over all rho): measS7={float(bm[1]):.6f} ({bm[1]}) "
              f"at rho={bm[0]:.4f}  margin={float(c-bm[1]):+.6f}  E={list(bm[2])}")
        print()
