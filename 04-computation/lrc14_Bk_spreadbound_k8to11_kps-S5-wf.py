#!/usr/bin/env python3
"""
LRC(14) ANGLE subtorus-relation-lattice -- SPREAD-BOUND PLATEAU k=8..11 (kps-S5).

Companion to the k=7 result (window-min minimized at W*=8, plateau after).  Here we
confirm the SAME plateau phenomenon for k=8,9,10,11: the window-minimum
   m(k,W) = min over E in {0,...,W}, |E|=k, {0,W} subset E, of mu(E)
is minimized at a SMALL W* (bounded spread ~ k+small) and m(k,W) STOPS DECREASING for
W > W* (it rises back).  We go a few steps past the apparent minimum to exhibit the
plateau, keeping the combinatorics feasible.

EXACT Fractions.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd, comb
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
    best, bestE = None, None
    inner = k - 2
    for combo in combinations(range(1, W), inner):
        E = (0,) + combo + (W,)
        if reduce(gcd, E) != 1: continue
        m = mu_exact(list(E))
        if best is None or m < best:
            best, bestE = m, E
    return best, bestE

if __name__ == "__main__":
    print("SPREAD-BOUND PLATEAU  m(k,W) for k=8..11 (exhibit the bounded W*):")
    # Wmax: a handful of steps past expected min (~k+3), capped by combinatorics.
    plan = {8:(7,16), 9:(8,16), 10:(9,16), 11:(10,16)}
    for k in (8, 9, 10, 11):
        Wlo, Whi = plan[k]
        print(f"\nk={k}:")
        running = None; runW = None; runE = None
        for W in range(Wlo, Whi+1):
            if comb(W-1, k-2) > 1200000:
                print(f"   W={W:2d}: (skipped C={comb(W-1,k-2)})"); continue
            m, E = min_over_window(k, W)
            if m is None: continue
            improved = (running is None) or (m < running)
            if improved: running, runW, runE = m, W, E
            print(f"   W={W:2d}: m = {str(m):>16s} = {float(m):.6f}"
                  f"{'  <-- new min' if improved else ''}", flush=True)
        print(f"   => k={k}: min m = {running} = {float(running):.6f} at W*={runW}, "
              f"E={runE} (spread {runW})")
    print("\nDONE.")
