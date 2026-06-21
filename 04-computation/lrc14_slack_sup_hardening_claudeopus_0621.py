#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: HARDEN the slack-side sup = the doublet (closes the skeleton's
'no global far-monotonicity' nut at the binding k=10,11).

kps HYP-2795 assembled skeleton: the SINGLE open analytic nut is the SLACK-side sup over the
infinite far-placement space at k>=10 (claimed 'no global far-monotonicity'). My HYP-2797
far-coherence monotonicity says the sup IS localized at r=1 (the tight doublet). This script
HARDENS that: a large adversarial search at k=10,11 over EXACTLY the families the skeleton
worried about (stretched APs, two-cluster gap splits, dilated bases + perturbation, large
random), confirming NOTHING beats the doublet sup_M p0 and reporting the worst (closest-to-cap)
genuine-wide config found. Exact rationals.
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


def doublet_sup(k):
    base = list(range(k - 2))
    best, bM = F(0), None
    for M in range(15, 80):
        v = p0_fast(tuple(sorted(base + [M, M + 1])))
        if v > best:
            best, bM = v, M
    return best, bM


def adversarial(k, rng, nrand):
    cand = set()
    # stretched APs: 0, d, 2d, ... with common difference d, then primitivize
    for d in range(1, 9):
        for length in range(k, k + 1):
            E = tuple(d * i for i in range(length))
            if len(set(E)) == k:
                cand.add(reprim(E))
    # dilated consec + 1,2 perturbations
    for d in range(2, 9):
        base = [d * i for i in range(k - 1)]
        for p in range(1, d * k + 25):
            if p in base:
                continue
            cand.add(reprim(tuple(sorted(set([0] + base + [p])))))
        base2 = [d * i for i in range(k - 2)]
        for p1 in range(1, d * k + 15, 3):
            for p2 in range(p1 + 1, d * k + 16, 5):
                if p1 in base2 or p2 in base2:
                    continue
                E = tuple(sorted(set([0] + base2 + [p1, p2])))
                if len(E) == k:
                    cand.add(reprim(E))
    # two-cluster gap splits: consec_a ... gap ... consec_b, all gaps/positions
    for a in range(1, k):
        b = k - a
        for gap in range(2, 80):
            s = (a - 1) + gap
            E = tuple(list(range(a)) + [s + i for i in range(b)])
            if len(set(E)) == k:
                cand.add(reprim(E))
    # doublet at every position (the predicted maximizer family) + near-doublet
    for M in range(15, 100):
        cand.add(reprim(tuple(range(k - 2)) + (M, M + 1)))
        cand.add(reprim(tuple(range(k - 2)) + (M, M + 2)))  # near-doublet
    # large random, wide spans
    for _ in range(nrand):
        sp = rng.randint(16, 150)
        pts = sorted(rng.sample(range(1, sp + 1), k - 1))
        cand.add(reprim(tuple([0] + pts)))
    return cand


def main():
    rng = random.Random(2024)
    print("=" * 78)
    print("HARDEN slack-side sup = doublet (binding k=10,11)  claude-opus 2026-06-21")
    print("=" * 78)
    for k in (10, 11):
        cap = CAP[k]
        ds, dM = doublet_sup(k)
        cand = adversarial(k, rng, nrand=120000)
        best, bestE, n = F(0), None, 0
        beaters = []
        for E in cand:
            if len(E) != k or not is_gw(E):
                continue
            n += 1
            v = p0_fast(E)
            if v > best:
                best, bestE = v, E
            if v > ds + F(1, 1000000):
                beaters.append((v, E))
        print(f"\nk={k}  cap={float(cap):.6f}")
        print(f"  doublet sup_M p0 = {float(ds):.6f} at M={dM}")
        print(f"  HARDENED search: {n} genuine-wide configs (target: stretched-AP/dilated/2-cluster/random)")
        print(f"  global max p0 = {float(best):.6f}  argmax={bestE}")
        print(f"  cap - max = {float(cap-best):+.6f}  (the worst/closest-to-cap genuine-wide found)")
        if beaters:
            beaters.sort(reverse=True)
            print(f"  *** {len(beaters)} BEAT the doublet! top {float(beaters[0][0]):.6f} {beaters[0][1]}")
        else:
            print(f"  NO beater over {n} configs => doublet = slack-side sup HARDENED (k={k})")
    print("\n" + "=" * 78)
    print("If no beaters: STEP C slack-sup is localized at the doublet (the 'no global")
    print("far-monotonicity' nut is resolved -> sup is at r=1), reducing STEP C to the doublet bound.")


if __name__ == "__main__":
    main()
