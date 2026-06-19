#!/usr/bin/env python3
"""
LRC(14) ANGLE subtorus-relation-lattice -- SPREAD-BOUND STABILITY (kps-S5).

THE crux of turning B(k) into a FINITE check: is there a spread cap D(k) beyond which
widening the window does NOT lower the bounded minimum of mu?  We test, for k=7..12,
the perforated near-AP minimum as the window width W grows:
   m(k,W) = min over E subset {0,...,W}, |E|=k, 0,W in E, of mu(E).
If m(k,W) STABILIZES (stops decreasing) at some W* = D(k), the global min over ALL
spreads equals m(k,W*) -- a finite check -- PROVIDED widening never helps beyond W*.
We report m(k,W) for W = k-1, k, ..., up to a cap, and flag the W* where it stabilizes.

This is the empirical content of the spread bound.  (A rigorous D(k) still needs the
effective-equidistribution input; this script measures whether the bound exists and
how large it is.)

EXACT Fractions.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from functools import reduce

G0 = F(2, 7)

def mu_exact(E):
    E = sorted(set(int(e) for e in E)); k = len(E)
    if k == 1: return F(1)
    diffs = {E[i]-E[j] for i in range(k) for j in range(k) if E[i]-E[j] > 0}
    bps = {F(0), F(1)}
    for d in diffs:
        for t in range(0, d+1): bps.add(F(t, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        fr = [F(E[i])*mid - (F(E[i])*mid).__floor__() for i in range(k)]
        order = sorted(range(k), key=lambda i: fr[i])
        n = [(F(E[i])*mid).__floor__() for i in range(k)]
        cross = {a, b}
        for r in range(k):
            i1, i2 = order[r], order[(r+1)%k]; wrap = 1 if r == k-1 else 0
            slope = E[i2]-E[i1]; const = -n[i2]+n[i1]+wrap
            if slope != 0:
                xc = (G0-const)/slope
                if a < xc < b: cross.add(xc)
        cross = sorted(cross)
        for u, v in zip(cross, cross[1:]):
            if u == v: continue
            mm = (u+v)/2
            P = sorted(F(E[i])*mm - n[i] for i in range(k))
            gaps = [P[r+1]-P[r] for r in range(k-1)] + [P[0]+1-P[-1]]
            if max(gaps) > G0: total += (v-u)
    return total

def min_over_window(k, W):
    """min mu over E={0}U combo U {W}, combo subset {1..W-1}, |E|=k, primitive."""
    best, bestE = None, None
    inner = k - 2          # interior points besides 0 and W
    if inner < 0:
        return None, None
    for combo in combinations(range(1, W), inner):
        E = (0,) + combo + (W,)
        if reduce(gcd, E) != 1: continue
        m = mu_exact(list(E))
        if best is None or m < best:
            best, bestE = m, E
    return best, bestE

if __name__ == "__main__":
    print("="*72)
    print("SPREAD-BOUND STABILITY: m(k,W) = min mu over E in {0..W}, |E|=k, 0,W in E.")
    print("Does it stabilize (stop decreasing) at a finite W* = D(k)?")
    print("="*72)
    # keep the combinatorics feasible: C(W-1,k-2).  We go to W where this stays <~5e5.
    caps = {7:22, 8:20, 9:19, 10:18, 11:17, 12:17}
    for k in (7, 8, 9, 10, 11, 12):
        Wmax = caps[k]
        print(f"\nk={k}:")
        running = None; runW = None; stableSince = None
        from math import comb
        for W in range(k-1, Wmax+1):
            if comb(W-1, k-2) > 1500000:
                print(f"   W={W:2d}: (skipped, C={comb(W-1,k-2)} too large)")
                continue
            m, E = min_over_window(k, W)
            if m is None:
                continue
            improved = (running is None) or (m < running)
            if improved:
                running, runW = m, W
                stableSince = W
            tag = "  <-- new min" if improved else ""
            print(f"   W={W:2d}: m = {str(m):>16s} = {float(m):.6f}{tag}", flush=True)
        print(f"   => k={k}: best m = {float(running):.6f} = {running} first at W={runW} "
              f"(plateau from W>={runW})")
    print("\nDONE.")
