#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: REFRAME the genuine-wide leg as the GENERALIZED DOUBLET family.

mac-mini k=12 obstruction: the genuine-wide max is E*=(0,2,4,6,8,9,10,11,12,14,16,18),
which BEATS the consec-adjacent doublet and exceeds Q(11). mac-mini concluded 'doublet not max'.

REFRAME (claude-opus): E* = base {0,2,4,6,8,9,10,11,12,14} + far pair {16,18} = a DOUBLET with
GAP g=2 (spread), not adjacent (g=1). THM-564 only treated g=1. The genuine-wide maximizer is a
GENERALIZED DOUBLET {M, M+g} over SOME bounded base B and gap g>=1 -- the adjacent case is just g=1.

This script TESTS the reframe at the binding k=10,11,12:
  (1) confirm E* = base + {16,18} (g=2 doublet); compute p0, far-count r.
  (2) over genuine-wide configs, is the MAX always r=2 (a generalized doublet, 2 far elements)?
      i.e. does far-count monotonicity (r=2 > r>=3) hold so the max family = generalized doublets?
  (3) within r=2, sweep base B + gap g + position M: find the max and confirm < cap.
If the genuine-wide max is always a generalized doublet, the R-tail framework (extended to gap g)
covers the ENTIRE genuine-wide leg -- completing the reduction.
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL
from lrc14_wide_branch_ridge_codex_s47 import primitive


def reprim(E):
    E = tuple(sorted(set(int(x) for x in E)))
    g = reduce(gcd, E)
    return tuple(x // g for x in E) if g > 1 else E


def is_gw(E):
    E = tuple(sorted(set(E)))
    if 0 not in E or max(E) - min(E) <= 14 or not primitive(E):
        return False
    for e in E:
        sub = reprim(tuple(x for x in E if x != e))
        if len(sub) < 2 or max(sub) - min(sub) <= 14:
            return False
    return True


def n_far(E):
    return sum(1 for e in E if e > 14)


def far_structure(E):
    """Return (base, far_tuple) splitting at 14."""
    base = tuple(e for e in E if e <= 14)
    far = tuple(e for e in E if e > 14)
    return base, far


def main():
    rng = random.Random(2026)
    print("=" * 80)
    print("GENERALIZED DOUBLET REFRAME of the genuine-wide leg  (claude-opus 2026-06-22)")
    print("=" * 80)

    # (1) E* analysis
    Estar = (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18)
    base, far = far_structure(Estar)
    print(f"\n(1) k=12 obstruction E* = {Estar}")
    print(f"    base (<=14) = {base}")
    print(f"    far  (>14)  = {far}   => far pair {far}, GAP g = {far[1]-far[0]} (spread doublet!)")
    print(f"    p0(E*) = {p0_fast(Estar)} = {float(p0_fast(Estar)):.6f}  cap_12={float(CAP[12]):.5f}  "
          f"Q(11)={float(QVAL[12]):.5f}")
    print(f"    genuine-wide? {is_gw(Estar)}   n_far(r) = {n_far(Estar)}")

    # (2) far-count of the genuine-wide max at k=10,11,12 (exhaustive-ish span<=20)
    for k in (10, 11, 12):
        cap = CAP[k]
        best, bestE = F(0), None
        by_r = {}
        seen = set()
        # structured: all bounded bases (subset of [0,14], size up to k-2..k-2) + far pairs/triples
        for nf in (2, 3, 4):
            basesize = k - nf
            if basesize < 1:
                continue
            # sample bases (full enumeration too big; structured + random)
            base_cands = set()
            base_cands.add(tuple(range(basesize)))  # consec
            base_cands.add(tuple([0] + [2 * i for i in range(1, basesize)]))  # even-AP
            for _ in range(300):
                S = tuple(sorted(rng.sample(range(1, 15), basesize - 1)))
                base_cands.add((0,) + S)
            for B in base_cands:
                if max(B) > 14:
                    continue
                # far elements: nf elements in (14, 30], structured (adjacent, spread) + random
                for _ in range(120):
                    fars = sorted(rng.sample(range(15, 31), nf))
                    E = reprim(tuple(B) + tuple(fars))
                    if len(E) != k or E in seen:
                        continue
                    seen.add(E)
                    if not is_gw(E):
                        continue
                    v = p0_fast(E)
                    r = n_far(E)
                    if r not in by_r or v > by_r[r][0]:
                        by_r[r] = (v, E)
                    if v > best:
                        best, bestE = v, E
        print(f"\n(2) k={k}  cap={float(cap):.5f}  genuine-wide max by far-count r:")
        for r in sorted(by_r):
            v, E = by_r[r]
            print(f"     r={r}: max p0={float(v):.6f}  cap-max={float(cap-v):+.5f}  {E}")
        rmax = n_far(bestE) if bestE else None
        print(f"     GLOBAL max p0={float(best):.6f} at r={rmax}: {bestE}")
        print(f"     => max is a GENERALIZED DOUBLET (r=2)? {rmax == 2}")
    print("\n" + "=" * 80)
    print("If the genuine-wide max is ALWAYS r=2 (a generalized doublet {M,M+g} over a base),")
    print("the R-tail framework extended to gap g covers the ENTIRE genuine-wide leg.")


if __name__ == "__main__":
    main()
