#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_L7_balanced_2cluster_window_THREADA_opus.py   (THREAD A, opus)

L7 = THE SOLE GAP: BALANCED MULTI-CLUSTER. The remaining LRC(14) gap after
[L1-L3 reduction] + [L4 single-far collar] + [L5,L6] is the balanced two-cluster:

    E = cluster1 (offsets near f1) UNION cluster2 (offsets near f2),   ratio rho = f2/f1.

measS7(E) = meas{ x in [0,1) : { floor(7 frac(e x)) : e in E } = Z/7 }
          = the joint distribution of x -> (frac(f1 x), frac(f2 x)) covering Z/7.

Geometrically x->(frac(f1 x),frac(f2 x)) traces a line of slope rho=f2/f1 on the
2-torus; the cover is the line's visit to the surjective region.

This script:
 (1) Verifies the (1, 2.15) ratio-window claim: for balanced two-cluster shapes E
     (two consecutive blocks, block1 at base f1, block2 at base f2), compute measS7
     over the ratio rho = f2/f1 and confirm it drops below cap above ~2.15, and find
     the actual WORST ratio (max measS7) and its margin to cap.
 (2) Records the worst balanced two-cluster shape and its exact measS7 vs cap.

All measures are EXACT Fractions (dilation-invariant: scale E by 1/gcd; measS7 uses
exact breakpoints at multiples of 1/(7e)).

cap_k (from THM-534 finite cert, k=8..13):
    8:2243/5880, 9:1979/4004, 10:55/91, 11:66/91, 12:6/7, 13:1
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce
sys.path.insert(0, '04-computation')
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

# ---- exact measS7 (matches lrc14_sector_certificate exactly) ----
def measS7(E):
    E = sorted(set(int(e) for e in E if e != 0))
    bps = set([F(0), F(1)])
    for e in E:
        ae = abs(e)
        for m in range(0, 7 * ae + 1):
            bps.add(F(m, 7 * ae))
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
    g = reduce(gcd, [abs(e) for e in E if e != 0], 0)
    return g == 1

def reduce_E(E):
    g = reduce(gcd, [abs(e) for e in E if e != 0], 0) or 1
    return tuple(sorted(set(e // g for e in E)))

# ---- balanced two-cluster builder ----
# cluster1 = consecutive block of size s starting at a1: {a1, a1+1, ..., a1+s-1}
# cluster2 = consecutive block of size s starting at a2: {a2, ..., a2+s-1}
# total k = 2s.  base f1 ~ a1 (center), f2 ~ a2.  ratio rho ~ a2/a1.
def two_block(a1, a2, s):
    return tuple(sorted(set(list(range(a1, a1 + s)) + list(range(a2, a2 + s)))))


if __name__ == "__main__":
    print("=" * 78)
    print("(1) BALANCED 2-CLUSTER: two consecutive blocks, ratio sweep, vs cap")
    print("=" * 78)
    print("    E = {a1..a1+s-1} U {a2..a2+s-1};  k=2s; rho=a2/a1 (block-base ratio)")
    print()

    for s in [4, 5]:          # block size; k=2s = 8, 10
        k = 2 * s
        c = CAP[k]
        print(f"  --- block size s={s}, k={k}, cap_{k}={float(c):.5f} ({c}) ---")
        # anchor block1 at a1=1 (cluster around small base), sweep a2 to sweep ratio
        # ratio rho = a2/a1 ; to get fine ratio resolution, also vary a1.
        best = None
        rows = []
        for a1 in range(1, 9):
            # block2 must be disjoint-ish; sweep a2 from a1+1 upward (rho>1)
            for a2 in range(a1 + 1, a1 * 6 + s + 4):
                E = two_block(a1, a2, s)
                if len(E) != k:        # overlap collapsed the size; skip (not balanced 2-cluster)
                    continue
                if not primitive(E):
                    continue
                m = measS7(E)
                rho = F(a2, a1)
                over = m > c
                rows.append((m, float(rho), a1, a2, E, over))
                if best is None or m > best[0]:
                    best = (m, rho, a1, a2, E)
        rows.sort(reverse=True)
        print(f"    scanned {len(rows)} balanced shapes. TOP 8 by measS7:")
        for (m, rho, a1, a2, E, over) in rows[:8]:
            flag = "  *** OVER CAP ***" if over else ""
            print(f"      measS7={float(m):.5f}  rho~{float(rho):.3f}  a1={a1},a2={a2}  E={list(E)}{flag}")
        # how measS7 trends with ratio: bucket by rounded ratio at fixed a1=1
        print(f"    [trend at a1=1] (rho=a2):")
        for a2 in range(2, 14):
            E = two_block(1, a2, s)
            if len(E) != k or not primitive(E): continue
            m = measS7(E)
            tag = "OVER" if m > c else "ok"
            print(f"        rho={a2}: measS7={float(m):.5f}  margin_to_cap={float(c-m):+.5f}  [{tag}]")
        bm, brho, ba1, ba2, bE = best
        print(f"    >>> WORST (max measS7) balanced 2-cluster: measS7={float(bm):.6f} ({bm})")
        print(f"        rho={brho}~{float(brho):.4f}  a1={ba1},a2={ba2}  E={list(bE)}")
        print(f"        cap_{k}={float(c):.6f}  MARGIN={float(c-bm):+.6f} ({c-bm})  "
              f"{'<<< BREACH' if bm>c else 'OK (under cap)'}")
        print()
