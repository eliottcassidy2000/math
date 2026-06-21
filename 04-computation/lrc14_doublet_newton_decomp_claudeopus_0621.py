#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: EXACT NEWTON DECOMPOSITION of the doublet error
(explains the "almost periodic" structure; pins |M*error| to known single-far period-max).

mac-mini Thread-B Newton peel:  p0(B u {M,M+1}) = p0(B) + [p0(B u{M})-p0(B)]
  + [p0(B u{M+1})-p0(B)] + C(M,M+1),   C = joint curvature (2nd difference).
=> error(M) := p0(E_M) - Phi_2 = Delta_M + Delta_{M+1} + (C(M) - C_sat),
   Delta_w := p0(B u {w}) - Phi(B)  [single-far deviation, THM-563: w*Delta_w periodic],
   Phi(B) = single-far plateau, Phi_2 = doublet plateau, C_sat = lim C.

PREDICTION (claude-opus): M*error(M) = (M*Delta_M) + (M*Delta_{M+1}) + M*(C(M)-C_sat)
  = [periodic, THM-563] + [~periodic] + [curvature approach, decaying].
This EXPLAINS HYP-2797's 'almost periodic' (periodic single-far parts + decaying curvature),
and gives  |M*error| <= 2*period-max(B) + sup_M |M*(C(M)-C_sat)|.
If the curvature-approach term is bounded/small, the doublet error is pinned to the KNOWN
single-far period-max (THM-563), closing the genuine-wide binding case via existing machinery.

This script (exact rationals, base B = consec {0,...,k-3}):
  - verifies the Newton identity error(M) == Delta_M + Delta_{M+1} + (C(M)-C_sat) exactly
  - reports sup_M |M*Delta_M| (single-far, should ~ THM-563 period-max ~1)
  - reports sup_M |M*(C(M)-C_sat)| (the curvature-approach residual)
  - checks 2*period-max + sup|M*curv-approach| vs the measured sup|M*error| (~1.4)
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP


def base_consec(k):
    return tuple(range(k - 2))  # {0,...,k-3}, k-2 elements incl 0


def main():
    print("=" * 80)
    print("EXACT NEWTON DECOMPOSITION of the doublet error  (claude-opus 2026-06-21)")
    print("error(M) = Delta_M + Delta_{M+1} + (C(M)-C_sat);  base B = consec {0,...,k-3}")
    print("=" * 80)
    LO, HI = 15, 400
    LATE = 300
    for k in range(8, 12):
        B = base_consec(k)
        p0B = p0_fast(B)
        # single-far values and plateau
        sf = {}
        for w in range(LO, HI + 2):
            sf[w] = p0_fast(tuple(sorted(B + (w,))))
        PhiB = sum(sf[w] for w in range(LATE, HI + 1)) / (HI + 1 - LATE)  # single-far plateau
        # doublet values and plateau
        dv = {}
        for M in range(LO, HI + 1):
            dv[M] = p0_fast(tuple(sorted(B + (M, M + 1))))
        Phi2 = sum(dv[M] for M in range(LATE, HI + 1)) / (HI + 1 - LATE)
        # curvature C(M) = dv[M] - p0B - (sf[M]-p0B) - (sf[M+1]-p0B) = dv[M] - sf[M] - sf[M+1] + p0B
        C = {M: dv[M] - sf[M] - sf[M + 1] + p0B for M in range(LO, HI + 1)}
        Csat = sum(C[M] for M in range(LATE, HI + 1)) / (HI + 1 - LATE)
        # verify identity: error(M) == Delta_M + Delta_{M+1} + (C(M)-Csat)
        maxres = F(0)
        for M in range(LO, HI + 1):
            err = dv[M] - Phi2
            DeltaM = sf[M] - PhiB
            DeltaM1 = sf[M + 1] - PhiB
            recon = DeltaM + DeltaM1 + (C[M] - Csat)
            # NOTE: Phi2 vs (p0B + 2*(PhiB-p0B) + Csat) consistency check
            res = abs(err - recon)
            if res > maxres:
                maxres = res
        # the identity is exact up to the plateau-estimate consistency:
        # Phi2 should equal p0B + 2*(PhiB-p0B) + Csat
        Phi2_pred = p0B + 2 * (PhiB - p0B) + Csat
        # signed sups
        supMdelta = max(abs(M * (sf[M] - PhiB)) for M in range(LO, HI + 1))
        supMcurv = max(abs(M * (C[M] - Csat)) for M in range(LO, HI + 1))
        supMerr = max(abs(M * (dv[M] - Phi2)) for M in range(LO, HI + 1))
        print(f"\nk={k}  B=consec{B}")
        print(f"  identity max residual |err - (DM+DM1+(C-Csat))| = {float(maxres):.2e} "
              f"(0 iff plateaus consistent; Phi2 vs pred diff = {float(Phi2-Phi2_pred):.2e})")
        print(f"  Csat (saturated joint curvature) = {float(Csat):+.5f}   "
              f"(mac-mini Thread-B sup|C|~0.029)")
        print(f"  sup_M |M*Delta_M| (single-far)        = {float(supMdelta):.4f}   (THM-563 period-max ~1)")
        print(f"  sup_M |M*(C(M)-Csat)| (curv approach)  = {float(supMcurv):.4f}")
        print(f"  => bound 2*sup|M*Delta| + sup|M*curv|  = {float(2*supMdelta+supMcurv):.4f}")
        print(f"     measured sup_M |M*error|            = {float(supMerr):.4f}  "
              f"(decomposition {'consistent' if supMerr <= 2*supMdelta+supMcurv+F(1,100) else 'CHECK'})")
    print("\n" + "=" * 80)
    print("READING: the doublet |M*error| is pinned to (single-far period-max, THM-563/known)")
    print("plus a curvature-approach residual. If the latter is bounded (signed Dedekind),")
    print("the genuine-wide doublet closes by existing single-far machinery + a curvature lemma.")


if __name__ == "__main__":
    main()
