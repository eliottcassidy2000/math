#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_L7_genuine_2cluster_THREADA_opus.py   (THREAD A, opus)

L7 SOLE GAP -- GENUINE balanced two-cluster (clusters SEPARATED, not abutting).

A *genuine* balanced two-cluster has TWO clusters at DIFFERENT SCALES f1, f2 with
a real gap.  The clean model: pick a small base "cluster shape" C (a few consecutive
offsets), and dilate it to two scales:

    E = f1 * C  UNION  f2 * C ,   rho = f2/f1  (>1).

Then x -> (frac(f1 x), frac(f2 x)) is a line of slope rho on the 2-torus; cluster1
reads the frac(f1 x) coordinate at fine sub-sector resolution (because C spreads it),
cluster2 reads frac(f2 x).  The cover (all 7 sectors hit) = the line's visit to the
surjective region.

NB the earlier two-ABUTTING-block model degenerated to consec (blocks touch -> single
run).  Here clusters are forced apart by the scale ratio rho.  When rho is an INTEGER
and C={1..s}, f1*C and f2*C interleave; the genuine separated case is rho NON-integer
or C not starting at 1.

We sweep rho = f2/f1 over rationals in (1, 3) at fine resolution and report measS7 vs
cap, locate the worst ratio, and confirm the (1,2.15) separation claim.

EXACT Fractions throughout (E is an integer set after clearing denominators of f1,f2).
cap_k: 8:2243/5880, 9:1979/4004, 10:55/91, 11:66/91, 12:6/7, 13:1
"""
import sys, itertools
from math import comb, gcd
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

def build_2cluster(C, p, q):
    """E = p*C U q*C ; rho=q/p.  Clear denominators not needed (C,p,q ints)."""
    s1 = set(p * c for c in C)
    s2 = set(q * c for c in C)
    E = tuple(sorted((s1 | s2) - {0}))
    return E

if __name__ == "__main__":
    print("=" * 78)
    print("(1) GENUINE balanced 2-cluster  E = p*C U q*C ,  rho=q/p  swept over (1,3)")
    print("=" * 78)

    # base cluster shapes C (each gives a cluster of |C| offsets when dilated)
    cluster_shapes = {
        "C={1,2,3,4}_s4": (1, 2, 3, 4),
        "C={1,2,3,4,5}_s5": (1, 2, 3, 4, 5),
        "C={2,3,4,5}_s4off": (2, 3, 4, 5),
    }

    # ratio grid: rho = q/p with small p (so scales stay moderate), q>p, gcd handled
    # sweep p in 1..7, q in p+1..ceil(2.6 p)+1 to cover rho up to ~2.6, plus a few higher
    for cname, C in cluster_shapes.items():
        k = 2 * len(C)
        if k not in CAP:
            continue
        c = CAP[k]
        print(f"\n  ### {cname}: |C|={len(C)} -> k=2|C|={k}, cap_{k}={float(c):.5f} ({c}) ###")
        rows = []
        for p in range(1, 9):
            for q in range(p + 1, int(p * 3.2) + 2):
                if gcd(p, q) != 1:
                    continue  # rho already covered by reduced form (dilation-invariant)
                E = build_2cluster(C, p, q)
                if len(E) != k:
                    # the two dilated clusters overlapped: not a clean balanced 2-cluster
                    continue
                if not primitive(E):
                    continue
                m = measS7(E)
                rho = F(q, p)
                rows.append((m, rho, p, q, E))
        rows.sort(key=lambda r: (-r[0], r[1]))
        print(f"    {len(rows)} clean separated shapes. TOP 10 by measS7:")
        for (m, rho, p, q, E) in rows[:10]:
            flag = "  *** OVER CAP ***" if m > c else ""
            print(f"      measS7={float(m):.5f}  rho={rho}~{float(rho):.3f}  (p={p},q={q})  E={list(E)}{flag}")
        # margin curve as function of rho (sorted by rho), at smallest p giving each rho band
        print(f"    [measS7 vs rho, sorted by rho]:")
        for (m, rho, p, q, E) in sorted(rows, key=lambda r: r[1]):
            tag = "OVER" if m > c else "ok"
            print(f"        rho={float(rho):.3f} (p{p}/q{q}): measS7={float(m):.5f}  cap-m={float(c-m):+.5f} [{tag}]")
        if rows:
            bm, brho, bp, bq, bE = max(rows, key=lambda r: r[0])
            print(f"    >>> WORST genuine 2-cluster: measS7={float(bm):.6f} ({bm}) at rho={brho}~{float(brho):.4f}")
            print(f"        E={list(bE)}  cap_{k}={float(c):.6f}  MARGIN={float(c-bm):+.6f} ({c-bm})  "
                  f"{'BREACH' if bm>c else 'OK'}")
