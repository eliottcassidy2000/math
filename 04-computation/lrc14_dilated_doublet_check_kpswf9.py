#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- does the DILATED-base doublet close, and via which machinery?

The k=10 genuine-wide maximizer (HYP-2805) is {0,2,4,6,8,10,12,14,15,16} = (2*consec_8) u {15,16}.
Threads claim: (B) frozen law p0_inf=91711/230496, f0=223; (C) THM-564 P/R also applies.

CHECK 1: scale-invariance. p0(2*consec_8 u {15,16}) =? p0(consec_8 u {7.5, 8})  -- but 7.5 is
NOT an integer, so this is NOT a single-far dilated config; THM-563's dilated single-far extension
does NOT directly apply (the far pair is split across the dilation, not a single runner). Confirm
the dilated doublet needs its OWN P/R / frozen analysis (THM-564 with the dilated base), not THM-563.

CHECK 2: run the dilated doublet family (2*consec_8) u {f,f+1} for f=15..40 and confirm
  (a) p0 < cap_10 for all f>=15, (b) the sup is at the tight end, (c) the margin floor.

CHECK 3: confirm the binding genuine-wide value is bounded BELOW cap with a definite margin by
scanning the dilated doublet sup AND the split-pair sup. Report the worst (the true binding margin).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast
from lrc14_wide_branch_ridge_codex_s47 import CAP

cap10 = CAP[10]


def main():
    print("=" * 92)
    print("DILATED-base DOUBLET analysis (k=10), kind-pasteur kpswf9")
    print("=" * 92)
    print(f"cap_10 = {cap10} = {float(cap10):.6f}\n")

    # CHECK 1: scale-invariance shows the far pair would map to non-integer speeds
    base_dil = tuple(2 * i for i in range(8))    # {0,2,...,14}
    print(f"dilated base = {base_dil} = 2*consec_8  (gcd 2)")
    Emax = tuple(sorted(set(base_dil) | {15, 16}))
    print(f"k=10 maximizer E = {Emax}")
    print(f"  p0(E) = {p0_fast(Emax)} = {float(p0_fast(Emax)):.6f}; cap-p0 = {float(cap10-p0_fast(Emax)):.6f}")
    print(f"  Under x->x/2 scale-invariance, far {{15,16}} -> speeds {{7.5, 8.0}} (7.5 NON-INTEGER):")
    print(f"  => the dilated doublet is NOT a dilated SINGLE-far config; THM-563's dilated")
    print(f"     single-far extension does NOT close it. It needs THM-564 (P/R) with base=2*consec_8,")
    print(f"     OR the frozen-law route (HYP-2806), applied to the DILATED base.\n")

    # CHECK 2: dilated doublet family sup
    print("CHECK 2: dilated doublet (2*consec_8) u {f,f+1}, f=15..40:")
    sup = F(-1); sarg = None
    for f in range(15, 41):
        E = tuple(sorted(set(base_dil) | {f, f + 1}))
        if reduce(gcd, [e for e in E if e]) != 1:
            continue
        pv = p0_fast(E)
        if pv > sup:
            sup, sarg = pv, f
    print(f"  sup p0 over f=15..40 = {sup} = {float(sup):.6f} at f={sarg}")
    print(f"  cap - sup = {float(cap10-sup):.6f}  (< cap: {sup < cap10})")
    tight = tuple(sorted(set(base_dil) | {15, 16}))
    print(f"  tight end f=15: p0={float(p0_fast(tight)):.6f} (the maximizer)\n")

    # CHECK 3: report the true binding margin (worst over dilated doublet + the split-pair found)
    split = (0, 2, 4, 6, 8, 10, 12, 14, 16, 29)
    print("CHECK 3: true binding genuine-wide margin at k=10:")
    print(f"  dilated tight doublet {tight}: margin {float(cap10-p0_fast(tight)):.6f}")
    print(f"  split-pair {split}: p0={float(p0_fast(split)):.6f}, margin {float(cap10-p0_fast(split)):.6f}")
    worst = min(cap10 - p0_fast(tight), cap10 - p0_fast(split))
    print(f"  WORST (binding) margin = {float(worst):.6f}")
    print(f"  => robust 0.16 {'FAILS' if worst < F(16,100) else 'holds'}; "
          f"< cap {'holds' if worst > 0 else 'FAILS'} (margin > 0)")
    print("=" * 92)


if __name__ == "__main__":
    main()
