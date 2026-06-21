#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: FROZEN ROOM Phi_frozen(B,g) < cap over ALL bounded bases & gaps.

The generalized-doublet reframe (HYP-2807) reduces the genuine-wide leg to:
  p0(B u {M,M+g}) = Phi_frozen(B,g) + g(M)/M,  g(M)=P+R bounded by G ~ 3.4 (uniform).
Closure for M >= M* = ceil(G/(cap-Phi)) is automatic IF Phi_frozen(B,g) < cap (the ROOM).
This script computes Phi_frozen(B,g) = lim_{M->inf} p0(B u {M,M+g}) (large-M estimate, averaged
over a window to kill the almost-periodic wobble) for ALL primitive bounded bases B (size k-2)
and gaps g=1..5, and reports the WORST (largest) Phi_frozen and cap - worst. If margin>0 for all
(B,g): the frozen room holds, and (with the uniform R-tail) the genuine-wide leg closes via a
finite window [15, M*]. Exact rationals.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP


def phi_frozen(B, g, lo=300, hi=340):
    # average p0 over a late M-window to estimate the frozen limit
    vals = [p0_fast(tuple(sorted(B + (M, M + g)))) for M in range(lo, hi + 1)]
    return sum(vals) / len(vals)


def bounded_bases(size):
    for S in combinations(range(1, 15), size - 1):
        yield (0,) + S


def main():
    print("=" * 78)
    print("FROZEN ROOM Phi_frozen(B,g) < cap over ALL bounded bases & gaps  claude-opus 0622")
    print("=" * 78)
    for k in (9, 10):
        cap = CAP[k]
        size = k - 2
        worst = {}  # per gap
        worst_overall = F(0)
        wo_at = None
        nb = 0
        for B in bounded_bases(size):
            nb += 1
            for g in (1, 2, 3, 4, 5):
                ph = phi_frozen(B, g)
                if g not in worst or ph > worst[g][0]:
                    worst[g] = (ph, B)
                if ph > worst_overall:
                    worst_overall, wo_at = ph, (B, g)
        print(f"\nk={k}  cap={float(cap):.6f}  bases={nb}  worst Phi_frozen per gap:")
        for g in sorted(worst):
            ph, B = worst[g]
            print(f"   g={g}: max Phi_frozen={float(ph):.6f}  cap-Phi={float(cap-ph):+.6f}  base={B}")
        print(f"   OVERALL worst Phi_frozen={float(worst_overall):.6f} at base={wo_at[0]} g={wo_at[1]}")
        print(f"   cap - worst = {float(cap-worst_overall):+.6f}  "
              f"=> frozen room {'HOLDS' if cap-worst_overall>0 else 'FAILS'}")
        # implied cutoff with G~3.4
        Gbound = F(34, 10)
        Hk = cap - worst_overall
        if Hk > 0:
            import math
            Mstar = math.ceil(float(Gbound) / float(Hk))
            print(f"   with G~3.4: H_K=cap-Phi={float(Hk):.4f} => M*=ceil(G/H_K)={Mstar} (finite window [15,{Mstar}])")
    print("\n" + "=" * 78)
    print("If frozen room HOLDS for all (B,g): genuine-wide closes = [frozen room] + [uniform")
    print("R-tail G~3.4] + [finite window [15,M*]]. M* small => tractable THM-563-style check.")


if __name__ == "__main__":
    main()
