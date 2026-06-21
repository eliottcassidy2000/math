#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: BROAD genuine-wide stress test at the BINDING k=10,11.

Goal: strengthen the combinatorial claim 'the genuine-wide p0 maximizer is
consec_{k-2} + tight far doublet {M,M+1} at small M' (HYP-2797), by an aggressive
adversarial search over MANY genuine-wide families and confirming none beats the
doublet's sup_M p0. If something beats it -> CRITICAL finding (report loudly).

genuine-wide: span>14 AND remove-any-one stays wide (reprimitivized). Exact rationals.

Families searched (k=10,11):
  - all 2-cluster splits [a, k-a] over block gaps and second-block positions (wide range)
  - all 3-cluster and 4-cluster splits, varied gaps
  - dilated AP d*consec + 1..2 perturbations
  - 'comb' configs (evenly spaced blocks)
  - large random genuine-wide (span up to 120), many samples
Reports the global max p0, its config, and the doublet baseline; flags any beater.
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


def doublet_baseline(k):
    base = list(range(k - 2))
    best, bestM = F(0), None
    for M in range(15, 60):
        v = p0_fast(tuple(sorted(base + [M, M + 1])))
        if v > best:
            best, bestM = v, M
    return best, bestM


def gen_all(k, rng):
    cand = set()
    # 2-cluster [a, k-a]
    for a in range(1, k):
        b = k - a
        for start in range(a + 1, 90):
            E = tuple(list(range(a)) + [start + i for i in range(b)])
            if len(set(E)) == k:
                cand.add(reprim(E))
    # 3-cluster
    for a in range(1, k - 1):
        for bb in range(1, k - a):
            c = k - a - bb
            if c < 1:
                continue
            for g1 in range(2, 30, 2):
                for g2 in range(2, 30, 2):
                    s2 = (a - 1) + g1
                    s3 = s2 + (bb - 1) + g2
                    E = tuple(list(range(a)) + [s2 + i for i in range(bb)] + [s3 + i for i in range(c)])
                    if len(set(E)) == k:
                        cand.add(reprim(E))
    # 4-cluster (coarse)
    for g in range(3, 16, 2):
        sizes = [k // 4] * 4
        for i in range(k - sum(sizes)):
            sizes[i] += 1
        E, pos = [], 0
        for s in sizes:
            E += [pos + i for i in range(s)]
            pos += s - 1 + g
        E = tuple(E)
        if len(set(E)) == k:
            cand.add(reprim(E))
    # dilated AP + perturbation
    for d in (2, 3, 4, 5):
        base = [d * i for i in range(k - 1)]
        for p in range(1, d * (k - 1) + 20):
            if p in base:
                continue
            E = tuple(sorted(set([0] + base + [p])))
            if len(E) == k:
                cand.add(reprim(E))
    # large random
    for _ in range(20000):
        sp = rng.randint(16, 120)
        pts = sorted(rng.sample(range(1, sp + 1), k - 1))
        cand.add(reprim(tuple([0] + pts)))
    return cand


def main():
    rng = random.Random(77)
    print("=" * 78)
    print("BROAD GENUINE-WIDE STRESS TEST (binding k=10,11)  claude-opus 2026-06-21")
    print("=" * 78)
    for k in (10, 11):
        cap = CAP[k]
        db, dbM = doublet_baseline(k)
        cand = gen_all(k, rng)
        best, bestE, nchk = F(0), None, 0
        beaters = []
        for E in cand:
            if len(E) != k or not is_gw(E):
                continue
            nchk += 1
            v = p0_fast(E)
            if v > best:
                best, bestE = v, E
            if v > db + F(1, 100000):  # strictly beats doublet baseline
                beaters.append((v, E))
        print(f"\nk={k}  cap={cap}={float(cap):.5f}")
        print(f"  doublet baseline sup_M p0 = {float(db):.6f} at M={dbM}  (config consec_{k-2}+{{{dbM},{dbM+1}}})")
        print(f"  GLOBAL max over {nchk} genuine-wide configs = {float(best):.6f}")
        print(f"     argmax = {bestE}")
        print(f"  cap - global_max = {float(cap-best):+.6f}")
        if beaters:
            beaters.sort(reverse=True)
            print(f"  *** {len(beaters)} configs BEAT the doublet baseline! top: {float(beaters[0][0]):.6f} {beaters[0][1]}")
        else:
            print(f"  NO config beats the doublet baseline -> maximizer = consec+doublet CONFIRMED (this search)")
    print("\n" + "=" * 78)
    print("If no beaters at both k: strong support that bounding the doublet bounds all genuine-wide.")


if __name__ == "__main__":
    main()
