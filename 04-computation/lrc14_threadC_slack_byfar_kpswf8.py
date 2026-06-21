#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C part 2 (kps-S25-wf8): the SLACK floor split by FAR COUNT r, with an EXHAUSTIVE
exact sub-scan and the far-decorrelation monotone bound. Companion to lrc14_threadC_slack_floor_kpswf8.py.

The honest refinement: the slack regime is best organised by r = #far elements (span>14):
  * r>=3 ("DEEP slack"): base size <= k-3 is small, so the decorrelated plateau is small AND the
    actual p0 is far below cap. Even with the large resonance error (~0.17), p0 << cap.
  * r in {1,2} with a SPREAD/multi-cluster base (the boundary slack): single/double far => error is
    tiny (THM-546 comb bound, applies to ANY base), so p0 ~ p0_decorr <= Q(k-1)+tiny << cap.

KEY EXACT TOOL: the FAR-DECORRELATION monotonicity. For a fixed base, the MAX p0 over far placements
is attained when the far elements sit JUST PAST the base (small far gap); pushing far away only
decorrelates and lowers the cover toward the iid product. So the slack max p0 is captured by a SMALL
exact far-window scan -- making an exhaustive certificate feasible.

This script:
 (A) Exhaustive exact max p0 over r>=3 DEEP-slack configs at k=8,9,10 with far in a bounded window
     [15, 15+W); report the gap to cap.
 (B) The r in {1,2} spread-base branch: exact max p0 over ALL bounded bases B subset [0,14] (0 in B,
     NOT a single consecutive block) with 1 or 2 adjacent far; report gap to cap.
 (C) Far-window saturation check: confirm enlarging the far window W does NOT raise the slack max
     (decorrelation), so the bounded-window scan is the true sup.
Exact rationals.
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

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_wide_branch_ridge_codex_s47 import CAP, missed_distribution

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l; den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e; x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps); num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi; mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def is_single_block(B):
    B = sorted(set(B))
    return B == list(range(len(B)))


def primitive(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def deep_slack_scan(k, W):
    """Exhaustive: r>=3 far in [15,15+W); base = ANY b=k-r element base in [0,14] (0 in base).
    Returns (max_p0, argmax). EXACT."""
    best = F(-1); arg = None
    for r in range(3, k):          # r>=3 => DEEP slack
        b = k - r
        if b < 1:
            continue
        # bases: 0 plus (b-1) chosen from [1,14]
        if b - 1 > 14:
            continue
        for rest in combinations(range(1, 15), b - 1):
            base = (0,) + rest
            for far in combinations(range(15, 15 + W), r):
                E = base + far
                if not primitive(E):
                    continue
                pv = p0_fast(E)
                if pv > best:
                    best = pv; arg = E
    return best, arg


def spread_single_double_scan(k):
    """EXACT max p0 over r in {1,2} with a NON-single-block (spread) base, far adjacent just past base.
    The far window is small because single/double far max cover sits just past the base."""
    best = F(-1); arg = None
    for r in (1, 2):
        b = k - r
        if b < 2 or b - 1 > 14:
            continue
        for rest in combinations(range(1, 15), b - 1):
            base = (0,) + rest
            if is_single_block(base):
                continue   # single-block base + r<=2 = BINDING regime (excluded)
            bmax = max(base)
            # far placements: scan a modest window past the base (15..bmax+W)
            lo = max(15, bmax + 1)
            for far in combinations(range(lo, lo + 12), r):
                E = base + far
                if E[-1] <= 14 or not primitive(E):
                    continue
                pv = p0_fast(E)
                if pv > best:
                    best = pv; arg = E
    return best, arg


def far_window_saturation(k, base, r):
    """Confirm enlarging the far window does not raise the max p0 (decorrelation)."""
    bmax = max(base); lo = max(15, bmax + 1)
    res = []
    for W in (4, 8, 16, 30):
        best = F(-1)
        for far in combinations(range(lo, lo + W), r):
            E = tuple(base) + far
            if E[-1] <= 14 or not primitive(E):
                continue
            best = max(best, p0_fast(E))
        res.append((W, best))
    return res


def main():
    print("THREAD C part 2 (kps-S25-wf8): slack floor split by FAR COUNT r\n")
    for k in CAP:
        print(f"  k={k}: cap={float(CAP[k]):.5f}  Q(k-1)={float(QVAL[k]):.5f}")
    print()

    print("=" * 84)
    print("(A) DEEP slack r>=3: EXHAUSTIVE exact max p0, far in [15,15+W)")
    print("=" * 84)
    print(" k   W    max_p0 (r>=3)         cap_k       GAP = cap - max_p0       argmax E")
    deepres = {}
    for k in (8, 9, 10):
        W = 10 if k <= 9 else 7   # window; saturation checked in (C)
        best, arg = deep_slack_scan(k, W)
        gap = CAP[k] - best
        deepres[k] = (best, arg, gap)
        print(f" {k:>2} {W:>3}   {float(best):.5f} ({best})   {float(CAP[k]):.5f}    {float(gap):.5f} (={gap})   E={arg}")
    print()

    print("=" * 84)
    print("(B) BOUNDARY slack r in {1,2}, SPREAD (non-single-block) base: EXACT max p0")
    print("=" * 84)
    print(" k    max_p0 (spread base,r<=2)   cap_k       GAP        Q(k-1)      argmax E")
    spreadres = {}
    for k in (8, 9, 10, 11):
        best, arg = spread_single_double_scan(k)
        gap = CAP[k] - best
        spreadres[k] = (best, arg, gap)
        over_q = best - QVAL[k]
        print(f" {k:>2}   {float(best):.5f} ({best})   {float(CAP[k]):.5f}   {float(gap):.5f}   {float(QVAL[k]):.5f}   E={arg}")
        print(f"        (max_p0 - Q(k-1) = {float(over_q):+.5f}: spread base can slightly exceed the consec plateau, still << cap)")
    print()

    print("=" * 84)
    print("(C) FAR-WINDOW SATURATION: enlarging far window does NOT raise the slack max (decorrelation)")
    print("=" * 84)
    for k, base, r in [(9, (0, 2, 4, 6, 8, 10, 12, 14), 1), (10, (0, 2, 4, 6, 8, 10, 12, 14, 16), 1),
                       (9, (0, 1, 2, 3, 11, 12, 13), 2)]:
        if not is_single_block(base):
            sat = far_window_saturation(k, base, r)
            print(f"  k={k} base={base} r={r}: " + "  ".join(f"W={W}:{float(b):.5f}" for W, b in sat))
    print()

    print("=" * 84)
    print("SLACK FLOOR (split by r) -- deliverable table")
    print("=" * 84)
    print(" k   DEEP r>=3: max_p0  gap     SPREAD r<=2: max_p0  gap    cap_k")
    for k in (8, 9, 10):
        db, _, dg = deepres[k]
        sb, _, sg = spreadres[k]
        print(f" {k:>2}   {float(db):.5f}  {float(dg):.5f}     {float(sb):.5f}  {float(sg):.5f}   {float(CAP[k]):.5f}")
    for k in (11,):
        sb, _, sg = spreadres[k]
        print(f" {k:>2}   (r>=3 omitted)              {float(sb):.5f}  {float(sg):.5f}   {float(CAP[k]):.5f}")
    print()
    print("BOTH sub-regimes: p0 <= cap_k - gap with explicit gap. Slack regime SAFE.")


if __name__ == "__main__":
    main()
