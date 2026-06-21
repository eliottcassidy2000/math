#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_span_monotone_apriori_macmini.py  (mac-mini 2026-06-21, THREAD D audit)

THE MISSING A-PRIORI BOUND for the GLOBAL consec-max claim.

The consec-max claim is: over ALL primitive k-sets of positive integers (UNBOUNDED span),
consec [1..k] is the global measS7-maximizer. Every exhaustive check is span-bounded.
For the global claim we need an a-priori reason wide sets cannot beat consec.

CANDIDATE a-priori bound (testable): the MAX measS7 over primitive k-sets of span exactly s
is NON-INCREASING in s for s >= k-1 (consec has minimal span s=k-1). If true, only finitely
many spans can possibly beat consec, reducing the global claim to a finite check.

This script computes, for k=8 and k=9, the function
   g(s) = max{ measS7(E) : E primitive, min(E)=1, max(E)=1+s, |E|=k }
for s = k-1, k, k+1, ... up to a feasible cap, and checks monotonicity + that
g(k-1) = measS7(consec) is the global max.

FAST: fixes min=1, max=1+s, chooses the middle k-2 from {2..s}. For k=8 that's C(s-1,6),
manageable to s~20.  Reports g(s) and the argmax witness.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def measS7(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e * mid) % 1) for e in E)) == P:
            tot += (b - a)
    return tot

def g_of_span(k, s):
    """max measS7 over primitive k-sets with min=1, max=1+s."""
    best = Fr(-1); bestE = None
    lo, hi = 2, s  # middle elements from {2..s} (so max=1+s), choose k-2
    for mid in itertools.combinations(range(2, 1 + s), k - 2):
        E = (1,) + mid + (1 + s,)
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1: continue
        v = measS7(E)
        if v > best: best = v; bestE = E
    return best, bestE

def main():
    print("#"*80)
    print("# SPAN-MONOTONE A-PRIORI BOUND for global consec-max")
    print("#"*80)

    for k, smax in [(8, 18), (9, 16)]:
        print(f"\n{'='*60}\nk = {k}   (consec span = {k-1})\n{'='*60}")
        consec = list(range(1, k + 1)); pc = measS7(consec)
        print(f"  consec [1..{k}] measS7 = {float(pc):.6f}  (span {k-1})")
        gs = []
        for s in range(k - 1, smax + 1):
            g, E = g_of_span(k, s)
            gs.append((s, g, E))
            flag = ""
            if g > pc + Fr(1, 10**12): flag = "  <<< BEATS CONSEC"
            print(f"  span={s:2d}: g(s)={float(g):.6f}  witness={E}{flag}")
        # monotonicity from s=k-1
        vals = [float(g) for _, g, _ in gs]
        mono = all(vals[i] >= vals[i+1] - 1e-12 for i in range(len(vals)-1))
        amax = gs[vals.index(max(vals))]
        print(f"  -> g(s) non-increasing in s? {mono}")
        print(f"  -> global max over computed spans at span={amax[0]} = {float(amax[1]):.6f} "
              f"({'consec' if amax[0]==k-1 else 'NOT consec span!'})")

    print("\n=== INTERPRETATION ===")
    print("If g(s) is non-increasing for s>=k-1, the global consec-max claim reduces to")
    print("a FINITE span window: only spans where g(s) could still exceed consec need check,")
    print("and monotonicity guarantees none beyond a cutoff can. This is the missing")
    print("a-priori ingredient. If g(s) is NON-monotone, the global claim needs more care.")
    print("\nDONE.")

if __name__ == "__main__":
    main()
