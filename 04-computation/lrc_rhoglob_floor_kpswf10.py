#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_rhoglob_floor_kpswf10.py  (kind-pasteur 2026-06-21, THM-527 Thread A)

THE CORRECTED OBJECT: rho*_glob = meas{ x in G_P AND maxgap{frac(e_i x)} > 1/7 }
                                = the GLOBAL-WITNESS density (gap > 1/7),
NOT rho*_crit (gap > 2/7, the via-Vmax criterion density, which we proved hits 0
on admissible covering sets).

THM-527 part A's "good period" uses gap>2/7 because it certifies via the Vmax
CRITERION (remove Vmax, one Vmax-period hosts a full-S-safe subarc).  But LRC(14)
only needs M(S)>=1/14, i.e. ONE global witness tau* safe for ALL of S -- which
(slow-fast, THM-526 lines 306-310) needs only gap > 1/7 (w_theta>0), NOT >2/7.
So the RIGHT density to floor is rho*_glob.

This script tests the DECISIVE question for the corrected program:

   Is  inf rho*_glob(P,E) > 0  over the shape space?
   In particular, does rho*_glob ever = 0 on ADMISSIBLE shapes?

If rho*_glob > 0 uniformly (and =0 only on inadmissible shapes), the corrected
compactness floor closes the GLOBAL-WITNESS existence -- which is what LRC needs.
We:
  (1) exhaustively scan k=9 shapes (spread<=14) x worst-P for rho*_glob=0;
  (2) compare the rho*_crit=0 census (38 cases) -- do they all have rho*_glob>0?
  (3) push spread-growing families to see the rho*_glob floor;
  (4) report the min rho*_glob and whether any admissible shape gives 0.

The pigeonhole (THM-526): with m=|L|=k cluster points, maxgap >= 1/k, so gap>1/7
is AUTOMATIC when k <= 6 (1/k >= 1/6 > 1/7).  The hard regime is k >= 7, where
gap>1/7 is NOT automatic -- exactly where rho*_glob could in principle hit 0.
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import gcd
from functools import reduce


def circ_maxgap_at(E, x):
    phases = sorted(set((Fr(e) * x) % 1 for e in E))
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def gp_breaks(P):
    bps = set()
    for p in P:
        if p == 0:
            continue
        for m in range(0, p):
            for r in (1, 13):
                v = Fr(14 * m + r, 14 * p)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def good_breaks(E):
    bps = set()
    diffs = set()
    El = list(E)
    for i in range(len(El)):
        for j in range(i + 1, len(El)):
            d = abs(El[i] - El[j])
            if d != 0:
                diffs.add(d)
    for d in diffs:
        for t in range(1, d):
            bps.add(Fr(t, d))
        for m in range(0, 7 * d + 1):
            for s in (2, -2):
                v = Fr(7 * m + s, 7 * d)
                if 0 < v < 1:
                    bps.add(v)
    return bps


def in_GP(P, x, thr=Fr(1, 14)):
    for p in P:
        f = (Fr(p) * x) % 1
        dd = f if f <= Fr(1, 2) else 1 - f
        if dd < thr:
            return False
    return True


def density(P, E, gapthr):
    bps = sorted({Fr(0), Fr(1)} | gp_breaks(P) | good_breaks(E))
    tot = Fr(0)
    for x0, x1 in zip(bps, bps[1:]):
        mid = (x0 + x1) / 2
        if in_GP(P, mid) and circ_maxgap_at(E, mid) > gapthr:
            tot += (x1 - x0)
    return tot


def main():
    print("=" * 78)
    print("THM-527 Thread A: the CORRECTED floor -- inf rho*_glob (gap>1/7) > 0 ?")
    print("=" * 78)
    sys.stdout.flush()
    g27 = Fr(2, 7)
    g17 = Fr(1, 7)

    # (1) exhaustive k=9 scan: rho*_crit=0 vs rho*_glob=0
    print("\n--- (1) k=9, spread<=14, 3 worst P: census of rho*_crit=0 vs rho*_glob=0 ---")
    worstP = [[1, 2, 3, 12], [1, 6, 7, 13], [1, 2, 3, 13]]
    crit0 = 0
    glob0 = 0
    tot = 0
    min_glob = None
    min_glob_arg = None
    glob0_examples = []
    for tail in itertools.combinations(range(1, 15), 8):
        E = [0] + list(tail)
        for P in worstP:
            tot += 1
            rc = density(P, E, g27)
            if rc == 0:
                crit0 += 1
                rg = density(P, E, g17)
                if rg == 0:
                    glob0 += 1
                    if len(glob0_examples) < 8:
                        glob0_examples.append((E, P))
                if min_glob is None or rg < min_glob:
                    min_glob = rg
                    min_glob_arg = (E, P)
    print(f"    rho*_crit=0 cases: {crit0} / {tot}")
    print(f"    of those, rho*_glob=0 cases: {glob0}")
    print(f"    min rho*_glob over the rho*_crit=0 cases = {min_glob} = {float(min_glob):.6f}")
    print(f"      at E={min_glob_arg[0]} P={min_glob_arg[1]}")
    if glob0_examples:
        print(f"    rho*_glob=0 examples (CHECK ADMISSIBILITY):")
        for E, P in glob0_examples:
            print(f"      E={E} P={P}")
    sys.stdout.flush()

    # (2) FULL k=9 scan for rho*_glob=0 (not just within crit=0)
    print("\n--- (2) FULL k=9 (spread<=14, 3 worst P): ANY rho*_glob=0? ---")
    g0 = 0
    g0ex = []
    minall = None
    minall_arg = None
    for tail in itertools.combinations(range(1, 15), 8):
        E = [0] + list(tail)
        for P in worstP:
            rg = density(P, E, g17)
            if rg == 0:
                g0 += 1
                if len(g0ex) < 8:
                    g0ex.append((E, P))
            if minall is None or rg < minall:
                minall = rg
                minall_arg = (E, P)
    print(f"    rho*_glob=0 cases: {g0}")
    print(f"    min rho*_glob (all) = {minall} = {float(minall):.6f} at "
          f"E={minall_arg[0]} P={minall_arg[1]}")
    for E, P in g0ex:
        print(f"      glob0: E={E} P={P}")
    sys.stdout.flush()

    # (3) is the min rho*_glob attained on an ADMISSIBLE shape?
    print("\n--- (3) admissibility of the min-rho*_glob shape ---")
    E, P = minall_arg
    sp = max(E)
    cov = 0
    Mbad = 0
    for Vmax in range(14 + sp, 14 + sp + 200):
        L = [Vmax - e for e in E]
        if len(set(L)) != len(L) or min(L) <= 13:
            continue
        S = sorted(set(P) | set(L))
        if len(S) != 13 or reduce(gcd, S) != 1:
            continue
        if all(any(v % q == 0 for v in S) for q in range(2, 15)):
            cov += 1
    print(f"    min-rho*_glob shape E={E} P={P} (rho*_glob={float(minall):.6f}):")
    print(f"    admissible covering+primitive S over Vmax sweep: {cov}")
    print(f"    => the min global-witness density is {'>0' if minall>0 else '=0'} "
          f"and {'IS' if cov>0 else 'is NOT'} admissible.")
    sys.stdout.flush()

    # (4) spread-growing: does rho*_glob floor positive?
    print("\n--- (4) rho*_glob along spread-growing families (floor positive?) ---")
    Pf = [1, 2, 3, 12]
    print("    F1 {0..7,M}:")
    for M in [9, 13, 23, 47, 101, 401]:
        E = list(range(8)) + [M]
        print(f"      M={M:4d}: rho*_glob = {float(density(Pf, E, g17)):.6f}")
    print("    AP {0,s,..,8s}:")
    for s in [11, 23, 43, 101]:
        E = [s * j for j in range(9)]
        print(f"      s={s:4d}: rho*_glob = {float(density(Pf, E, g17)):.6f}")

    print("\n" + "=" * 78)
    print("VERDICT (corrected program):")
    if g0 == 0:
        print("  rho*_glob (gap>1/7, the GLOBAL-WITNESS density) is POSITIVE on")
        print("  EVERY k=9 shape scanned -- including all rho*_crit=0 cases.")
        print("  => The CORRECT object to floor is rho*_glob, NOT rho*_crit.")
        print("     The via-Vmax CRITERION (>2/7) was the wrong target; the")
        print("     global-witness density (>1/7) is what LRC needs and it stays >0.")
    else:
        print(f"  rho*_glob=0 occurs ({g0} cases) -- must check admissibility;")
        print("  if any admissible, the global-witness route also has gaps.")
    print("=" * 78)
    print("\nDONE.")


if __name__ == "__main__":
    main()
