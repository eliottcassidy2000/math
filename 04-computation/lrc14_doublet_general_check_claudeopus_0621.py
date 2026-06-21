#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: DOUBLET GENERAL CHECK over ALL bounded bases (the genuine-wide
analogue of THM-563's 12805-base single-far check). Closes kps-S27's named remaining task:
'genuine-wide leg reduces to all-bounded-bases doublet generalization (THM-563-style)'.

By far-count monotonicity (HYP-2797), the genuine-wide binding family is the tight far DOUBLET
{M,M+1} appended to a bounded base B subset [0,14]. E = B u {M,M+1} is always genuine-wide
(removing any one element keeps span>14) and primitive (M,M+1 consecutive => gcd 1).

This script verifies, for EVERY primitive bounded base B (0 in B, |B|=k-2, B subset [0,14]) and
all M in the window [15, F], the EXACT max p0(B u {M,M+1}) < cap_k, and reports the binding
(closest-to-cap) base + margin. The tail M>F is covered by kps's p0 <= p0_inf + J/f bound
(decaying); here F=80 is well past the M=20-21 argmax (p0 decays toward p0_inf<cap thereafter).

Output: per-k worst base, worst max-p0, margin cap-max. If all margins>0: genuine-wide doublet
is VERIFIED over ALL bounded bases (not just consec) -- the THM-563-style general check.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP
from lrc14_wide_branch_ridge_codex_s47 import primitive


def bounded_bases(size):
    """All B = {0} u S, S subset {1..14}, |B|=size, primitive."""
    for S in combinations(range(1, 15), size - 1):
        B = (0,) + S
        # primitive as a base on its own (gcd of nonzero elements)
        if reduce(gcd, B[1:]) == 1 or True:  # base gcd handled by full-config primitivity
            yield B


def main():
    F_WIN = 80
    print("=" * 80)
    print("DOUBLET GENERAL CHECK over ALL bounded bases (genuine-wide, THM-563-style)")
    print("claude-opus 2026-06-21   E = B u {M,M+1}, B subset [0,14], |B|=k-2, M in [15,80]")
    print("=" * 80)
    for k in range(8, 11):
        cap = CAP[k]
        size = k - 2
        worst = F(0)
        worst_base = None
        worst_M = None
        nbases = 0
        nfail = 0
        for B in bounded_bases(size):
            nbases += 1
            bmax = F(0)
            bM = None
            for M in range(15, F_WIN + 1):
                E = tuple(sorted(B + (M, M + 1)))
                if not primitive(E):
                    continue
                v = p0_fast(E)
                if v > bmax:
                    bmax, bM = v, M
            if bmax >= cap:
                nfail += 1
            if bmax > worst:
                worst, worst_base, worst_M = bmax, B, bM
        print(f"\nk={k}  cap={cap}={float(cap):.6f}   bases checked: {nbases}  (size {size})")
        print(f"  WORST (closest-to-cap) genuine-wide doublet: max p0 = {float(worst):.6f}")
        print(f"     at base B={worst_base}, doublet M={worst_M}")
        print(f"  margin cap - worst = {float(cap - worst):+.6f}   FAILS (>=cap): {nfail}")
        print(f"  => doublet general check {'PASSES (all bounded bases < cap)' if nfail == 0 else 'FAILS'}")
    print("\n" + "=" * 80)
    print("If all PASS: the genuine-wide doublet is VERIFIED over ALL bounded bases (window")
    print("[15,80]); tail M>80 by kps p0<=p0_inf+J/f. Combined w/ far-count monotonicity (doublet")
    print("= genuine-wide max), this is the THM-563-style general check for the genuine-wide leg.")


if __name__ == "__main__":
    main()
