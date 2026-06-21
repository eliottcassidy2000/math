#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- BROADENED genuine-wide maximizer (refines HYP-2805 further).

The synthesis verification found that at k=10 a SPLIT far pair {16,29} on a dilated base
{0,2,..,14} gives p0=154283/341040=0.452390 > the tight dilated doublet {15,16}=265/588=0.450680
> the consec doublet=0.442517. So the genuine-wide maximizer is NOT even an adjacent doublet.

This script: for k=8..12, exhaustive over (k-2)-subset bases of [0,14] containing 0 (NO
base-primitivity filter) x ALL far PAIRS {fa,fb}, 15<=fa<fb<=FHI (not just adjacent), filter
genuine(FULL primitive, >=2 far, not single-perturbation-reducible). Report the true max,
the binding margin, and CRITICALLY whether p0 < cap ALWAYS holds (the LRC requirement).

This is the honest 'is the doublet the genuine-wide max, and what is the true binding margin'
answer for question (3) of the synthesis.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive


def nfar(E):
    return sum(1 for e in E if e > 14)


def is_genuine(E):
    nz = [e for e in E if e]
    if nfar(E) < 2 or reduce(gcd, nz) != 1:
        return False
    for fr in [e for e in E if e > 14]:
        if max(x for x in E if x != fr) <= 14:
            return False
    return True


def adjacent(fa, fb):
    return fb == fa + 1


def broad_max(k, FHI):
    bs = k - 2
    cap = CAP[k]
    best = F(-1); bestE = None
    best_adj = F(-1); bestE_adj = None   # restricted to adjacent doublets (THM-564 family)
    capfail = 0
    farpairs = [(a, b) for a in range(15, FHI + 1) for b in range(a + 1, FHI + 1)]
    for rest in combinations(range(1, 15), bs - 1):
        B = (0,) + rest
        for (fa, fb) in farpairs:
            E = tuple(sorted(set(B) | {fa, fb}))
            if len(E) != k or not is_genuine(E):
                continue
            pv = p0_fast(E)
            if pv >= cap:
                capfail += 1
            if pv > best:
                best, bestE = pv, E
            if adjacent(fa, fb) and pv > best_adj:
                best_adj, bestE_adj = pv, E
    return best, bestE, best_adj, bestE_adj, cap, capfail


def main():
    print("=" * 100)
    print("BROADENED genuine-wide maximizer: split far PAIRS, dilated bases  (kind-pasteur kpswf9)")
    print("=" * 100)
    # FHI per k: keep search tractable. k=8,9 cheap -> FHI=30. k=10 FHI=30. k=11,12 FHI=24.
    FHI = {8: 30, 9: 30, 10: 30, 11: 24, 12: 24}
    print(f"  {'k':>3} {'broad max p0':>16} {'~':>9} {'cap':>9} {'margin':>9} {'>=.16':>6} "
          f"{'<cap':>5} | {'adj-doublet max':>16} {'~':>9} {'adj margin':>10}  binding-config")
    worst = (F(2), None); allcap = True
    rows = []
    for k in (8, 9, 10, 11, 12):
        best, bestE, best_adj, bestE_adj, cap, capfail = broad_max(k, FHI[k])
        m = cap - best; madj = cap - best_adj
        allcap = allcap and (best < cap) and (capfail == 0)
        if m < worst[0]:
            worst = (m, k)
        split = not (len([e for e in bestE if e > 14]) == 2 and
                     sorted(e for e in bestE if e > 14)[1] == sorted(e for e in bestE if e > 14)[0] + 1)
        rows.append((k, best, bestE, m, best_adj, madj, split, capfail))
        print(f"  {k:>3} {str(best):>16} {float(best):>9.6f} {float(cap):>9.6f} {float(m):>9.6f} "
              f"{str(m>=F(16,100)):>6} {str(best<cap):>5} | {str(best_adj):>16} {float(best_adj):>9.6f} "
              f"{float(madj):>10.6f}  {bestE}{'  [SPLIT pair]' if split else '  [adjacent]'}")
    print()
    print(f"  p0 < cap at EVERY k (LRC requirement), 0 cap-violations? {allcap}")
    print(f"  WORST robust margin = {float(worst[0]):.6f} at k={worst[1]}  "
          f"{'(>= 0.16)' if worst[0]>=F(16,100) else '(< 0.16, but > 0 => < cap holds)'}")
    print()
    print("  KEY: at k=10 the binding config is a DILATED base + (possibly SPLIT) far pair.")
    print("  The adjacent doublet (THM-564 family) is NOT the genuine-wide max; the true")
    print("  binding margin is even tighter than the adjacent-doublet 0.1537.")
    print("=" * 100)


if __name__ == "__main__":
    main()
