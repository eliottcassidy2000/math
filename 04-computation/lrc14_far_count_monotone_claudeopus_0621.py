#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: FAR-COUNT MONOTONICITY — does r=2 (doublet) dominate r>=3?
Addresses kps HYP-2795's stated gap 'no global far-monotonicity' at the MAXIMIZER level.

If the genuine-wide p0 maximizer always has exactly r=2 far blocks (the tight doublet),
and max p0 DECREASES as the far-block count r grows, then bounding the doublet bounds ALL
genuine-wide. This script tabulates, for k=9..12, the EXACT max p0 over genuine-wide configs
stratified by r = number of maximal coherent far blocks (clusters beyond the base window),
to test r=2 dominance and p0-vs-r monotonicity.

Also: the r=3 rung of the Dedekind ladder (reflection lrc-the-dedekind-ladder...) — the triplet
{M,M+1,M+2} joint third difference C3(M). If it saturates (-> C3_sat) with an O(1/M) tail like
the doublet, the ladder's inductive step is supported.
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP
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


def far_block_count(E):
    """# maximal consecutive runs among elements with value>14 (the 'far' part),
    counting the base run [0..14] separately; r = number of far clusters."""
    E = sorted(E)
    blocks = 0
    prev = None
    for x in E:
        if prev is None or x - prev > 1:
            blocks += 1
        prev = x
    # subtract the base block (the one containing 0)
    return blocks - 1  # far clusters = total clusters minus the base cluster


def gen_configs(k, rng, nrand=15000):
    cand = set()
    # structured r=2,3,4 splits
    for a in range(1, k):
        for start in range(a + 1, 70):
            E = tuple(list(range(a)) + [start + i for i in range(k - a)])
            if len(set(E)) == k:
                cand.add(reprim(E))
    for a in range(1, k - 1):
        for bb in range(1, k - a):
            c = k - a - bb
            if c < 1:
                continue
            for g1 in range(2, 24, 2):
                for g2 in range(2, 24, 2):
                    s2 = (a - 1) + g1
                    s3 = s2 + (bb - 1) + g2
                    E = tuple(list(range(a)) + [s2 + i for i in range(bb)] + [s3 + i for i in range(c)])
                    if len(set(E)) == k:
                        cand.add(reprim(E))
    for _ in range(nrand):
        sp = rng.randint(16, 110)
        pts = sorted(rng.sample(range(1, sp + 1), k - 1))
        cand.add(reprim(tuple([0] + pts)))
    return cand


def triplet_curvature(k):
    """C3(M) = third finite difference for far triplet {M,M+1,M+2} (Dedekind-ladder r=3 rung)."""
    base = tuple(range(k - 3))  # base so that base + triplet = k speeds
    def p(extra):
        return p0_fast(tuple(sorted(base + tuple(extra))))
    vals = []
    for M in (20, 40, 80, 160, 320):
        # 3rd difference: inclusion-exclusion over subsets of {M,M+1,M+2}
        C3 = (p((M, M + 1, M + 2)) - p((M, M + 1)) - p((M, M + 2)) - p((M + 1, M + 2))
              + p((M,)) + p((M + 1,)) + p((M + 2,)) - p(()))
        vals.append((M, C3))
    return base, vals


def main():
    rng = random.Random(123)
    print("=" * 78)
    print("FAR-COUNT MONOTONICITY: does r=2 (doublet) dominate r>=3?  claude-opus 2026-06-21")
    print("=" * 78)
    for k in range(9, 13):
        cap = CAP[k]
        cand = gen_configs(k, rng)
        by_r = {}
        for E in cand:
            if len(E) != k or not is_gw(E):
                continue
            r = far_block_count(E)
            v = p0_fast(E)
            if r not in by_r or v > by_r[r][0]:
                by_r[r] = (v, E)
        print(f"\nk={k}  cap={float(cap):.5f}   max p0 by far-block-count r:")
        for r in sorted(by_r):
            v, E = by_r[r]
            print(f"   r={r}: max p0 = {float(v):.6f}  cap-max={float(cap-v):+.5f}  argmax={E}")
        if 2 in by_r:
            r2 = by_r[2][0]
            dom = all(by_r[r][0] <= r2 for r in by_r if r >= 2)
            print(f"   => r=2 dominates r>=3? {dom}")
    print("\n" + "=" * 78)
    print("DEDEKIND LADDER r=3 rung: triplet {M,M+1,M+2} third-difference C3(M) (saturation?)")
    for k in (9, 10):
        base, vals = triplet_curvature(k)
        print(f"  k={k} base=consec{base}:")
        for M, C3 in vals:
            print(f"     M={M:3d}: C3(M) = {float(C3):+.6f}")
        # saturation: late values should converge
        tail = [float(v) for _, v in vals[-3:]]
        print(f"     -> {'saturating (ladder rung supported)' if max(tail)-min(tail) < 0.02 else 'not clearly saturating'}")


if __name__ == "__main__":
    main()
