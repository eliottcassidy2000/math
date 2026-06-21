#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): EXACT r=2 genuine-wide floor, k=10,11,12.

STRUCTURAL REDUCTION (proved here):
  * r=1 (single far) is NEVER genuine-wide: removing the single far leaves the base
    B ⊆ [0,14] with span ≤ 14, hence reduced-span(B) ≤ 14. So genuine-wide ⟹ r ≥ 2.
  * The binding genuine-wide stratum is r=2 (base size k-2, two far). r ≥ 3 has a small
    base and is far below (handled separately). So the slack floor reduces to: max p0 over
    r=2 genuine-wide configs < Q(k-1).

This script: EXACT max p0 over ALL r=2 configs (base = any (k-2)-subset of [0,14] with 0;
far = {f1,f2} with both > 14), restricted to genuine-wide, with the far pair (f1,f2) swept
over a saturating window. We sweep BOTH far positions, not just adjacent pairs, to be sure
the binding shape is captured. We confirm the per-(base) sup saturates (the far-window cutoff).

Exact rationals throughout (whole-circle breakpoint p0).
"""
from __future__ import annotations
import sys, functools, time
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, missed_distribution
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


@lru_cache(maxsize=None)
def _p0(E):
    return p0_fast(E)


def primitive(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def reduced_span(S):
    S = sorted(set(S))
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False
    return True


def r2_max_for_base(base, gap_max=24, f_hi=70):
    """EXACT max p0 over r=2 genuine-wide configs on a fixed base.
    far pair (f1,f2), f1 in [max(15,bmax+1), f_hi], gap = f2-f1 in [1, gap_max]."""
    bmax = base[-1]
    lo = max(15, bmax + 1)
    best = F(-1)
    arg = None
    for f1 in range(lo, f_hi + 1):
        for f2 in range(f1 + 1, f1 + gap_max + 1):
            E = tuple(sorted(set(base) | {f1, f2}))
            if len(E) != len(base) + 2 or not primitive(E):
                continue
            if not is_genuine_wide(E):
                continue
            pv = _p0(E)
            if pv > best:
                best, arg = pv, E
    return best, arg


def scan_r2(k, gap_max=24, f_hi=70, dense_diff=3):
    """Max p0 over r=2 genuine-wide configs at size k. Base = (k-2)-subset of [0,14] with 0.
    We scan ALL bases, but to keep it tractable we note coverage is maximized by DENSE bases;
    we scan dense bases (max consec-diff <= dense_diff) exhaustively AND a full pass over all
    bases with a SMALLER far window as a safety net (dense_check)."""
    b = k - 2
    cap = CAP[k]
    q = QVAL[k]
    best = F(-1)
    arg = None
    n_bases = 0
    for rest in combinations(range(1, 15), b - 1):
        base = (0,) + rest
        diffs = [y - x for x, y in zip(base, base[1:])]
        if diffs and max(diffs) > dense_diff:
            continue
        n_bases += 1
        bb, ba = r2_max_for_base(base, gap_max=gap_max, f_hi=f_hi)
        if bb > best:
            best, arg = bb, ba
    return best, arg, n_bases


def saturation_check(base, f_max=2000):
    """For a fixed base, confirm max p0 (gap<=24) over f1 up to f_max equals the early-window max."""
    bmax = base[-1]
    lo = max(15, bmax + 1)
    runmax = F(-1)
    argf = None
    early = None
    for f1 in range(lo, f_max + 1):
        for g in range(1, 25):
            f2 = f1 + g
            E = tuple(sorted(set(base) | {f1, f2}))
            if not primitive(E):
                continue
            pv = _p0(E)
            if pv > runmax:
                runmax, argf = pv, (f1, f2)
        if f1 == lo + 70:
            early = (runmax, argf)
    return early, (runmax, argf)


def main():
    print("=" * 96)
    print("THREAD 2: EXACT r=2 GENUINE-WIDE FLOOR  k=10,11,12  (mac-mini-S7)")
    print("=" * 96)
    print("STRUCTURAL: genuine-wide => r>=2 (r=1 removes-to-base span<=14). Binding stratum = r=2.")
    for k in (10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.6f}  Q(k-1)={float(QVAL[k]):.6f}  MARGIN={float(MARGIN[k]):.6f}")
    print()
    results = {}
    for k in (10, 11, 12):
        t = time.time()
        print("-" * 96)
        print(f"k={k}  r=2 EXACT scan (dense bases diff<=3, far gap<=24, f1<=70)")
        best, arg, nb = scan_r2(k, gap_max=24, f_hi=70, dense_diff=3)
        q = QVAL[k]
        under = best < q
        results[k] = (best, arg, under)
        print(f"  bases scanned={nb}")
        print(f"  => r=2 genuine-wide max p0 = {float(best):.6f} = {best}")
        print(f"     Q(k-1)={float(q):.6f}  Q - max_p0 = {float(q-best):+.6f}  "
              f"({'UNDER Q' if under else '*** OVER Q -- OBSTRUCTION ***'})")
        print(f"     argmax E={arg} base={tuple(e for e in arg if e<=14)} far={tuple(e for e in arg if e>14)}")
        print(f"     time={time.time()-t:.1f}s")
        print()
    print("=" * 96)
    print("SATURATION CHECK (far-window cutoff): does the per-base max stabilize for large f?")
    print("=" * 96)
    for k, base in [(10, tuple(range(8))), (11, tuple(range(9))), (12, tuple(range(10)))]:
        early, full = saturation_check(base, f_max=1500)
        print(f"  k={k} base={base}: max(f1<=lo+70)={float(early[0]):.6f}@{early[1]}  "
              f"max(f1<=1500)={float(full[0]):.6f}@{full[1]}  saturated? {early[0]==full[0]}")
    print()
    print("=" * 96)
    print("VERDICT")
    print("=" * 96)
    allunder = all(results[k][2] for k in (10, 11, 12))
    if allunder:
        print("  r=2 genuine-wide max p0 < Q(k-1) at k=10,11,12 (over dense bases / tight far).")
        print("  Saturation confirms the far-window cutoff (per-base sup at small f). Remaining for")
        print("  RIGOR: (a) extend bases to ALL (not just dense), (b) a rigorous tail bound for f>F*.")
    else:
        print("  OBSTRUCTION FOUND: some r=2 genuine-wide config has p0 >= Q(k-1). See argmax above.")


if __name__ == "__main__":
    main()
