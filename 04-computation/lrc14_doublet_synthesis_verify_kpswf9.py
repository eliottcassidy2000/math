#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- SYNTHESIS verification of the genuine-wide DOUBLET leg.

Independent cross-check of the three threads' load-bearing numbers:
  (1) p0(consec_8 u {f,f+1}) at k=10, tight end f=15..21 -> consec-doublet sup = 1301/2940.
  (2) p0({0,2,4,6,8,10,12,14,15,16}) = 265/588 (the DILATED-base k=10 maximizer, threads B/C),
      and confirm it BEATS the consec doublet => robust margin 0.16 FAILS at k=10 (0.15372).
  (3) BROADER k=10 genuine-wide search: non-adjacent far pairs, far range [15,30], many base
      shapes -- does ANYTHING beat 265/588? does < cap ever fail? (the real LRC requirement).
  (4) Recompute p0 by a 2ND independent integrator (direct rational sweep of the step function
      x|->#sectors hit) to certify p0_fast is not the source of the numbers.

All EXACT rationals.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive

ALL_INNER = 0b1111110


def p0_independent(E):
    """2nd integrator. measure{ x in [0,1) : {floor(7 e x) mod 7 : e in E} covers sectors 1..6 }.

    Sector of element e at point x is floor(7 e x) mod 7 (= which of the 7 arcs e*x lands in,
    scaled). E 'covers' iff the 6 nonzero sectors 1..6 all appear. Breakpoints of the piecewise-
    constant cover indicator are at x = j/(7 e) for all e, j. Integrate exactly over each cell.
    """
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    # common denominator for breakpoints x=j/(7e): D = 7*lcm(e)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    D = 7 * l
    bps = {0, D}
    for e in nz:
        # breakpoints at x = j/(7e) -> in units of 1/D that's j*(l//e)
        step = l // e
        x = 0
        while x <= D:
            bps.add(x)
            x += step
    bps = sorted(bps)
    num = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        # midpoint in units of 1/D ; x_mid = (lo+hi)/(2D)
        # sector of e = floor(7 e x_mid) mod 7 = floor(7 e (lo+hi)/(2D)) mod 7
        # 7 e (lo+hi)/(2D) = 7 e (lo+hi)/(2*7*l) = e (lo+hi)/(2 l)
        mask = 0
        for e in nz:
            mask |= 1 << ((e * (lo + hi) // (2 * l)) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num += hi - lo
    return F(num, D)


def nfar(E):
    return sum(1 for e in E if e > 14)


def is_genuine(E):
    nz = [e for e in E if e]
    if nfar(E) < 2 or reduce(gcd, nz) != 1:
        return False
    for fr in [e for e in E if e > 14]:
        if max(x for x in E if x != fr) <= 14:
            return False
    return True


def main():
    print("=" * 96)
    print("DOUBLET SYNTHESIS VERIFICATION (kind-pasteur kpswf9)")
    print("=" * 96)

    cap10 = CAP[10]
    print(f"\ncap_10 = {cap10} = {float(cap10):.6f}\n")

    # (1) consec doublet at k=10, tight end
    print("(1) CONSEC doublet consec_8 u {f,f+1}, k=10, tight end f=15..21:")
    base8 = tuple(range(8))
    consec_sup = F(-1); consec_arg = None
    for f in range(15, 22):
        E = tuple(sorted(set(base8) | {f, f + 1}))
        pv = p0_fast(E)
        pv2 = p0_independent(E)
        flag = "  <-- p0_fast != p0_indep!!" if pv != pv2 else ""
        if pv > consec_sup:
            consec_sup, consec_arg = pv, f
        print(f"    f={f:>2}: p0={str(pv):>12}={float(pv):.6f}  (indep {float(pv2):.6f}){flag}")
    print(f"  consec-doublet sup over f>=15 = {consec_sup} = {float(consec_sup):.6f} at f={consec_arg}")
    print(f"  cap - consec_sup = {cap10-consec_sup} = {float(cap10-consec_sup):.6f}  "
          f"(THREAD A's 'binding' margin 0.16188)\n")

    # (2) dilated-base maximizer
    print("(2) DILATED-base maximizer {0,2,4,6,8,10,12,14,15,16}:")
    Edil = (0, 2, 4, 6, 8, 10, 12, 14, 15, 16)
    pdil = p0_fast(Edil); pdil2 = p0_independent(Edil)
    print(f"    p0_fast      = {pdil} = {float(pdil):.6f}")
    print(f"    p0_independent = {pdil2} = {float(pdil2):.6f}   match={pdil==pdil2}")
    print(f"    is_genuine={is_genuine(Edil)}  primitive(full)={primitive(Edil)}  "
          f"primitive(base)={primitive(tuple(e for e in Edil if e<=14))}")
    print(f"    cap - p0 = {cap10-pdil} = {float(cap10-pdil):.6f}")
    print(f"    BEATS consec doublet by {float(pdil-consec_sup):.6f}; "
          f"robust 0.16 {'FAILS' if cap10-pdil < F(16,100) else 'holds'} "
          f"(margin {float(cap10-pdil):.5f}); < cap {'holds' if pdil<cap10 else 'FAILS'}\n")

    # (3) BROADER k=10 genuine-wide search: non-adjacent pairs, far range [15,30]
    print("(3) BROADER k=10 genuine-wide search (non-adjacent far pairs, far range [15,30],")
    print("    all 8-element bases from [0,14] with 0, gcd-1 full set):")
    best = F(-1); bestE = None; capfail = []
    farpairs = [(a, b) for a in range(15, 31) for b in range(a + 1, 31)]
    nbase = 0
    for rest in combinations(range(1, 15), 7):
        B = (0,) + rest
        nbase += 1
        for (fa, fb) in farpairs:
            E = tuple(sorted(set(B) | {fa, fb}))
            if len(E) != 10:
                continue
            if not is_genuine(E):
                continue
            pv = p0_fast(E)
            if pv >= cap10:
                capfail.append((E, pv))
            if pv > best:
                best, bestE = pv, E
    print(f"    searched {nbase} bases x {len(farpairs)} far pairs")
    print(f"    BROADER max p0 = {best} = {float(best):.6f}  at {bestE}")
    print(f"    cap - max = {float(cap10-best):.6f}")
    print(f"    beats 265/588? {best > pdil}   ( 265/588 = {float(pdil):.6f} )")
    print(f"    any p0 >= cap (LRC FAIL)? {len(capfail)>0}  ({len(capfail)} configs)")
    if capfail:
        for E, pv in capfail[:5]:
            print(f"       FAIL: {E}  p0={float(pv):.6f}")
    print()

    # (4) verdict on the 0.16 vs <cap question
    print("=" * 96)
    print("VERDICT:")
    print(f"  consec-doublet sup (THM-564 family) k=10 : {float(consec_sup):.6f}  margin {float(cap10-consec_sup):.5f} (>=0.16)")
    print(f"  TRUE genuine-wide max k=10 (this search)  : {float(best):.6f}  margin {float(cap10-best):.5f}")
    rob = (cap10 - best) >= F(16, 100)
    print(f"  robust margin >= 0.16 at k=10 ? {rob}   (binding config {bestE})")
    print(f"  p0 < cap at k=10 (the LRC requirement)   ? {best < cap10}")
    print("=" * 96)


if __name__ == "__main__":
    main()
