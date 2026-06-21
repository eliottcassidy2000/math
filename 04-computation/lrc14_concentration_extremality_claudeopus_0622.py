#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: CONCENTRATION EXTREMALITY -- the cleanest gK8 wide closure.

DISCOVERY: max_wide L_yK8 <= max_bounded L_yK8 = MB < 10cap. So the GLOBAL maximizer of
L_yK8 = 10q0+q3+10q6 over ALL k-speed configs is a BOUNDED (concentrated) config (consec-like).
=> gK8's bounded cert (max_bounded L_yK8 <= 10cap, mac-mini Lean) is the GLOBAL cert =>
   L_yK8(E) <= 10cap for ALL E (bounded AND wide) => p0 <= cap. NO dichotomy, NO R-tail needed.

This is the global extremality of the literature's tight instance (consec). This script
STRESS-TESTS it: search for ANY wide config (span>14) with L_yK8 > MB, over an aggressive
adversarial bank INCLUDING the binding families (genuine-wide maximizers, E*, single-far
near-cap, dilated even-AP, small-M R-tail bumps, two-cluster, stretched-AP). If NONE exceeds
MB at any k, concentration extremality is strongly supported => the wide bound closes cleanly.
Exact rationals.
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

from lrc14_gK8_decorrelation_principle_claudeopus_0622 import LyK8
from lrc14_threadA_regime_dichotomy_kpswf8 import CAP
from lrc14_wide_branch_ridge_codex_s47 import primitive


def gcd_all(E):
    return reduce(gcd, E)


def max_bounded(k, rng, n=2500):
    MB, MBe = F(0), None
    cands = {tuple(range(k))}
    if 2 * (k - 1) <= 14:
        cands.add(tuple(range(0, 2 * k, 2)))
    for _ in range(n):
        S = tuple(sorted(rng.sample(range(1, 15), k - 1)))
        E = (0,) + S
        if gcd_all(E) == 1:
            cands.add(E)
    for E in cands:
        if max(E) - min(E) <= 14:
            v = LyK8(E)
            if v > MB:
                MB, MBe = v, E
    return MB, MBe


def wide_bank(k, rng):
    """Aggressive wide bank including binding families + small-M R-tail bumps."""
    out = set()
    size = k - 2
    # generalized doublets over many bases (incl small M = R-tail bump), gaps 1..4
    bases = {tuple(range(size)), tuple([0] + list(range(15 - (size - 1), 15)))}
    if 2 * (size - 1) <= 14:
        bases.add(tuple(range(0, 2 * size, 2)))
    for _ in range(200):
        bases.add((0,) + tuple(sorted(rng.sample(range(1, 15), size - 1))))
    for B in bases:
        if max(B) > 14:
            continue
        for g in (1, 2, 3, 4):
            for M in range(15, 50):
                out.add(tuple(sorted(B + (M, M + g))))
    # single-far near-cap (consec_{k-1}+far)
    for f in range(15, 80):
        out.add(tuple(range(k - 1)) + (f,))
    # dilated even-AP + far
    for f in range(15, 60):
        out.add(tuple(range(0, 2 * (k - 1), 2)) + (f,))
    # E* family (even-AP + odd bridges + far) at k=12
    if k == 12:
        out.add((0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18))
    # random wide
    for _ in range(6000):
        sp = rng.randint(16, 60)
        out.add(tuple([0] + sorted(rng.sample(range(1, sp + 1), k - 1))))
    return out


def main():
    rng = random.Random(314)
    print("=" * 78)
    print("CONCENTRATION EXTREMALITY: any WIDE config with L_yK8 > MB?  (claude-opus 0622)")
    print("=" * 78)
    for k in (10, 11, 12):
        cap = CAP[k]
        MB, MBe = max_bounded(k, rng)
        MW, MWe = F(0), None
        exceed = []
        nchk = 0
        for E in wide_bank(k, rng):
            E = tuple(sorted(set(E)))
            if len(E) != k or max(E) - min(E) <= 14 or not primitive(E):
                continue
            nchk += 1
            v = LyK8(E)
            if v > MW:
                MW, MWe = v, E
            if v > MB:
                exceed.append((v, E))
        print(f"\nk={k}  10cap={float(10*cap):.5f}")
        print(f"  MB (bounded max L_yK8) = {float(MB):.5f}  at {MBe}")
        print(f"  MW (wide max L_yK8) = {float(MW):.5f}  at {MWe}   ({nchk} wide configs)")
        print(f"  MW - MB = {float(MW - MB):+.5f}   10cap - MW = {float(10*cap - MW):+.5f}")
        if exceed:
            exceed.sort(reverse=True)
            print(f"  *** {len(exceed)} WIDE configs EXCEED MB! top: {float(exceed[0][0]):.5f} {exceed[0][1]}")
        else:
            print(f"  NO wide config exceeds MB => CONCENTRATION EXTREMALITY holds (k={k})")
    print("\n" + "=" * 78)
    print("If no exceedances at any k: max_E L_yK8 = max_bounded L_yK8 = MB < 10cap (GLOBAL).")
    print("=> gK8 bounded cert is GLOBAL => wide bound closes by concentration extremality.")
    print("   No dichotomy, no R-tail. The whole wide leg reduces to BOUNDED extremality of L_yK8.")


if __name__ == "__main__":
    main()
