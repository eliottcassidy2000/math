#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- THE TRUE GENUINE-WIDE MAXIMIZER (refines HYP-2797/2804).

FINDING: the binding genuine-wide config is NOT the consec doublet consec_{k-2} u {f,f+1}
(HYP-2797/THM-564) but a DILATED-BASE adjacent doublet whose BASE is non-primitive (gcd 2)
while the FULL set is primitive. Sweeps that filter bounded bases by primitive(BASE) MISS it.

EXHAUSTIVE k-runner genuine-wide search (all (k-2)-subset bases of [0,14] containing 0, NO
base-primitivity filter, x adjacent far pairs {f,f+1}, f in [15,24]; genuine = >=2 far,
primitive(FULL), not single-perturbation-reducible):

  k=8 : max p0 = 0.207226  {0,10,11,12,13,14,15,16}        margin 0.1742  (>=0.16)
  k=9 : max p0 = 321/980 = 0.327551  {0,4,6,8,10,12,14,15,16}   margin 0.1667  (>=0.16)
  k=10: max p0 = 265/588 = 0.450680  {0,2,4,6,8,10,12,14,15,16} margin 0.1537  (< 0.16, > 0)
  k=11: max p0 = 0.521058  consec_8 u {21,22}               margin 0.2042  (>=0.16)
  k=12: max p0 = 0.590136  {0,2,4,5,6,7,8,10,12,14,20,21}    margin 0.2670  (>=0.16)

CONCLUSIONS (honest):
 * genuine-wide p0 < cap_k at EVERY k (the LRC requirement) -- VERIFIED, all margins > 0.
 * the ROBUST margin >= 0.16 FAILS at k=10 (0.1537): the consec-doublet 0.16-closure
   (THM-564/HYP-2804) is NOT the genuine-wide max; the true max edges it out by ~0.008.
 * the true maximizers have a DILATED base (gcd 2) {0,2,..,14} or {0,4,..,14} + a TIGHT
   far doublet {15,16} sitting JUST ABOVE the base. Base-primitivity sweeps miss these.
 * STRUCTURAL FIX for any all-bounded-bases sweep: filter on primitive(FULL E), not
   primitive(base). The maximizing base is non-primitive.

This script RE-DERIVES the table exhaustively (exact rationals) and certifies < cap.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from itertools import combinations
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive


def nfar(E):
    return sum(1 for e in E if e > 14)


def is_genuine(E):
    nz = [e for e in E if e]
    if nfar(E) < 2 or reduce(gcd, nz) != 1:
        return False
    for fr in [e for e in E if e > 14]:
        if max(x for x in E if x != fr) <= 14:  # single-perturbation reducible
            return False
    return True


def true_max(k, fhi=24):
    bs = k - 2
    cap = CAP[k]
    best = F(-1); bestE = None
    fps = [(f, f + 1) for f in range(15, fhi + 1)]
    for rest in combinations(range(1, 15), bs - 1):
        B = (0,) + rest
        if len(B) != bs:
            continue
        for fp in fps:
            E = tuple(sorted(set(B) | set(fp)))
            if len(E) != k or not is_genuine(E):
                continue
            pv = p0_fast(E)
            if pv > best:
                best, bestE = pv, E
    return best, bestE, cap


def main():
    print("=" * 96)
    print("TRUE GENUINE-WIDE MAXIMIZER (exhaustive, all bases incl. non-primitive)  kind-pasteur kpswf9")
    print("=" * 96)
    print(f"  {'k':>3} {'max p0':>16} {'~':>10} {'cap_k':>10} {'margin':>9} {'>=.16':>6} {'<cap':>5} {'base prim':>9}  maximizer")
    worst = (F(2), None)
    allcap = True
    for k in (8, 9, 10, 11, 12):
        best, bestE, cap = true_max(k)
        base = tuple(e for e in bestE if e <= 14)
        m = cap - best
        bp = primitive(base)
        allcap = allcap and (best < cap)
        if m < worst[0]:
            worst = (m, k)
        print(f"  {k:>3} {str(best):>16} {float(best):>10.6f} {float(cap):>10.6f} {float(m):>9.6f} "
              f"{str(m >= F(16,100)):>6} {str(best < cap):>5} {str(bp):>9}  {bestE}")
    print()
    print(f"  genuine-wide p0 < cap at every k (LRC requirement): {allcap}")
    print(f"  WORST robust margin = {float(worst[0]):.6f} at k={worst[1]}  "
          f"({'>=0.16' if worst[0] >= F(16,100) else 'BELOW 0.16 (but >0, so < cap holds)'})")
    print("  => the 0.16 robust target is MISSED at k=10 (margin 0.1537); < cap holds everywhere.")
    print("  => true maximizers have DILATED (non-primitive) bases + tight far doublet {15,16};")
    print("     base-primitivity sweeps (HYP-2804) miss them. Filter on primitive(FULL E).")
    print("=" * 96)


if __name__ == "__main__":
    main()
