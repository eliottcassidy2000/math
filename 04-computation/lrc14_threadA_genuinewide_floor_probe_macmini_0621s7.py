#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): probe the GENUINE-WIDE slack floor at k=10,11,12.

GOAL (deliverable): upgrade kps HYP-2788's genuine-wide slack floor
    every genuine-wide config E (span>14, remove-any-one-element STILL has reduced-span>14)
    has p0(E) < Q(k-1)
from VERIFIED (~8000 adversarial samples) to PROVED, at k=10,11,12, OR find the precise
obstruction.

This first script CHARACTERIZES the binding stratum precisely so the proof can target it:
  - r-stratified EXACT max p0 over genuine-wide configs, far window scanned over a full
    "period candidate" to find the true per-base sup.
  - identifies which (r, base-shape) is binding (closest to Q) at each k.

DEFINITIONS:
  E primitive, 0 in E, |E|=k, span=max(E)>14.
  base B = E ∩ [0,14] (contains 0), far F = {e>14}, r=|F|.
  genuine-wide := remove-ANY-single-element leaves reduced-span(rest)>14.
    reduced-span(S)=(max-min)/gcd(consecutive diffs).
  p0(E)=meas{x: E*x hits all 6 inner sectors}, EXACT rational (whole-circle breakpoint method).
  Q(k-1)=boundary_value_direct(consec_{k-1},1) = the decorrelated single-far plateau floor.

This is a PROBE, not the certificate. It maps the danger.
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
    """remove-any-one-element STILL leaves reduced-span > 14."""
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False
    return True


def Ple(B, r):
    """Deterministic upper bound: p0(B u F) <= sum_{t<=r} prof_t(B) (far covers <= r sectors)."""
    prof = missed_distribution(tuple(sorted(set(B))))
    return sum(prof[t] for t in range(min(r, 6) + 1))


# ----------------------------------------------------------------------------
# r-stratified exact scan. For each (base B, far shape S, far position f), build
# E = B u {f+s: s in S}. Sweep f over a window to find the per (B,S) sup.
# We restrict to bases that are dense enough to matter (longest run large) and
# to tight far shapes (peak coverage). Then sweep f over [15, 15+WIN).
# ----------------------------------------------------------------------------
def far_shapes(r, max_off=6):
    """Tight far-cluster offset shapes (start at 0). The coverage-maximizing far
    placements are TIGHT (adjacent / small-gap), per the resonance decay law."""
    if r == 1:
        return [(0,)]
    shapes = set()
    shapes.add(tuple(range(r)))            # adjacent block
    for g in (2, 3):
        s = tuple(i * g for i in range(r))
        if max(s) <= 7 * r:
            shapes.add(s)
    # all tight shapes within a small window (exhaustive small clusters)
    if r <= 3:
        for rest in combinations(range(1, max_off + 1), r - 1):
            shapes.add((0,) + rest)
    return sorted(shapes)


def scan_genuine_wide(k, win=70, base_dense_only=True, verbose=False):
    """Exact per (base, far-shape) sup over far position f in [15,15+win), restricted to
    genuine-wide configs. Returns max p0 per r and overall, with argmax."""
    cap = CAP[k]
    q = QVAL[k]
    best_overall = F(-1)
    arg_overall = None
    per_r = {}
    for r in range(1, k):
        b = k - r
        if b < 1 or b - 1 > 14:
            continue
        # bases: 0 + (b-1) from [1,14]. To keep it tractable & cover the binding shapes,
        # we take DENSE bases (max consecutive-diff <= 3) when base_dense_only, else all.
        best_r = F(-1)
        arg_r = None
        shapes = far_shapes(r)
        for rest in combinations(range(1, 15), b - 1):
            base = (0,) + rest
            if base_dense_only and b >= 2:
                diffs = [y - x for x, y in zip(base, base[1:])]
                if max(diffs) > 3:
                    continue
            bmax = base[-1]
            lo = max(15, bmax + 1)
            for S in shapes:
                for f in range(lo, lo + win):
                    far = tuple(f + s for s in S)
                    if min(far) <= 14:
                        continue
                    E = tuple(sorted(set(base) | set(far)))
                    if len(E) != k or not primitive(E):
                        continue
                    if not is_genuine_wide(E):
                        continue
                    pv = _p0(E)
                    if pv > best_r:
                        best_r, arg_r = pv, E
        per_r[r] = (best_r, arg_r)
        if best_r > best_overall:
            best_overall, arg_overall = best_r, arg_r
        if verbose and best_r >= 0:
            print(f"    r={r} base_size={b}: max_p0={float(best_r):.6f} "
                  f"(cap-p0={float(cap-best_r):.6f}, Q-p0={float(q-best_r):+.6f}) arg={arg_r}")
    return best_overall, arg_overall, per_r


def main():
    print("=" * 96)
    print("THREAD 2 PROBE: genuine-wide slack floor characterization k=10,11,12 (mac-mini-S7)")
    print("=" * 96)
    for k in (8, 9, 10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.6f}  Q(k-1)={float(QVAL[k]):.6f}  MARGIN={float(MARGIN[k]):.6f}")
    print()
    for k in (10, 11, 12):
        t = time.time()
        print("-" * 96)
        print(f"k={k}  (dense bases, tight far shapes, far window 70)")
        best, arg, per_r = scan_genuine_wide(k, win=70, base_dense_only=True, verbose=True)
        q = QVAL[k]
        print(f"  => genuine-wide max p0 = {float(best):.6f} = {best}")
        print(f"     Q(k-1) = {float(q):.6f},  Q - max_p0 = {float(q - best):+.6f}  "
              f"({'UNDER Q (good)' if best < q else '*** OVER Q ***'})")
        print(f"     argmax E = {arg}, base={tuple(e for e in arg if e<=14)}, far={tuple(e for e in arg if e>14)}")
        print(f"     time={time.time()-t:.1f}s")
        print()


if __name__ == "__main__":
    main()
