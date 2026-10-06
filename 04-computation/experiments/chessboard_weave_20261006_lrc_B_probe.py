#!/usr/bin/env python3
"""Task B probe (chessboard weave, LRC reading), 2026-10-06.  Exact arithmetic.

Is max tau stable in B?  First lonely time census on larger universes
(tau only, no q):  n=5 B=48, n=6 B=36, n=7 B=26 (all primitive n-subsets of {1..B}).
Also tests the observed construction "tight set T, remove one speed s, add a multiple
k(n+1) of n+1" (the added speed kills the tight witnesses a/(n+1)).
Reproduce: python3 chessboard_weave_20261006_lrc_B_probe.py > chessboard_weave_20261006_lrc_B_probe.out
"""
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from chessboard_weave_20261006_lrc_core import safe_components, primitive_sets, tau, check

T0 = time.time()
TOPK = 6
for n, B in [(5, 48), (6, 36), (7, 26)]:
    t1 = time.time()
    N = n + 1
    top = {}
    best_no1 = (Fraction(0), [])
    num = 0
    running = {}
    for v in primitive_sets(n, B):
        num += 1
        comps, D = safe_components(v)
        check(comps, v)
        t = Fraction(comps[0][0], D)
        check(2 * t < 1, v)
        if v[0] >= 2:
            if t > best_no1[0]:
                best_no1 = (t, [v])
            elif t == best_no1[0] and len(best_no1[1]) < 6:
                best_no1[1].append(v)
        if t in top:
            if len(top[t]) < 8:
                top[t].append(v)
        elif len(top) < TOPK or t > min(top):
            top[t] = [v]
            if len(top) > TOPK:
                del top[min(top)]
        m = v[-1]
        if m not in running or t > running[m][0]:
            running[m] = (t, v)
    print("=" * 78)
    print(f"n={n}, B={B}: {num} primitive sets ({time.time() - t1:.1f}s); all have 0 < tau < 1/2")
    for k in sorted(top, reverse=True):
        print(f"   tau = {k} = {float(k):.6f}: {top[k]}")
    print(f"   max tau with vmin >= 2: {best_no1[0]} = {float(best_no1[0]):.6f} at {best_no1[1]}")
    run, out = Fraction(0), []
    for m in range(1, B + 1):
        if m in running and running[m][0] > run:
            run = running[m][0]
            out.append(f"B>={m}: {run} {running[m][1]}")
    print("   record progression of max tau as B grows:")
    for s in out[-6:]:
        print("     ", s)

print("=" * 78)
print("Construction 'tight set T minus one speed s plus a multiple k(n+1)', k(n+1) <= 60:")
TIGHT = [(1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 3, 4, 7), (1, 2, 3, 4, 5), (1, 3, 4, 5, 9),
         (1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]
for T in TIGHT:
    n = len(T)
    N = n + 1
    best = (Fraction(0), None)
    for s in T:
        rest = [x for x in T if x != s]
        for k in range(1, 60 // N + 1):
            m = k * N
            if m in rest:
                continue
            v = tuple(sorted(rest + [m]))
            t = tau(v)
            if t > best[0]:
                best = (t, v)
    print(f"   T = {T}: best tau = {best[0]} = {float(best[0]):.6f} at {best[1]}")
print(f"\nTotal time {time.time() - T0:.1f}s")
