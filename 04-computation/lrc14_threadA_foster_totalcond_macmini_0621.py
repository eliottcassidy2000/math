#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2 PAYOFF: the FOSTER TOTAL-CONDUCTANCE
extremality. We found:
  * sumW is EXACTLY a function of the conductance vector c_r = sum_{e mod7=r, e!=0} 1/|e|.
  * consec UNIQUELY MAXIMIZES the TOTAL conductance T(E) = sum_r c_r over the
    full-residue stratum (0 shapes beat it).
This script tests whether the clean, PROVABLE statement
   "consec maximizes T(E) = sum_{e in E, e!=0} 1/|e|"
(a) holds at all strata k=8 / spans, and (b) COINCIDES with sumW-max, giving an
electrical proof route: T-max is a trivial harmonic-sum fact (consec uses the
SMALLEST possible magnitudes => largest 1/|e|), and IF sumW-max <=> T-max, the
wall reduces to this trivial fact.

NOTE the subtlety: T(E) = sum_{e!=0} 1/|e| does NOT depend on residues at all --
it is just the harmonic sum of the magnitudes. consec = {0,1,...,7} obviously
maximizes sum 1/|e| among any 8 distinct integers containing 0 (you cannot do
better than {1,...,7} for the nonzero part). So T-max is TRIVIALLY consec, over
ANY stratum (not just full-residue!). The CONTENT is whether sumW-max = T-max.

TESTS:
 (T1) T(E) = sum_{e!=0} 1/|e|; consec is the unique T-max over ALL k=8 shapes
      with 0 (trivial). Confirm.
 (T2) Is sumW-max == T-max across the WHOLE k=8 stratum (not just full-residue)?
      If sumW is maximized by the same shape that maximizes T, that is the
      electrical reduction. Test exhaustively span<=W.
 (T3) Is sumW a MONOTONE function of T alone (rank correlation)? Or does sumW
      need the full c-vector? (We know c-vector determines it; does T suffice
      for the ARGMAX even if not for the whole function?)
 (T4) The DECISIVE refinement: restrict to shapes with the SAME residue partition
      as consec (each residue once, plus residue 0 doubled). Within that exact
      "consec-type" topology, is sumW monotone in T = harmonic sum? (Rayleigh on
      a fixed network topology -- the clean electrical monotonicity.)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    secs = set()
    for e in E:
        v = e * xm; v = v - (v.numerator // v.denominator)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7_arcs(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    arcs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: arcs.append((lo, hi))
    return arcs

def W_a(E, a):
    arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14)
    w = F(0)
    for lo, hi in arcs:
        l = max(lo, lo_b); h = min(hi, hi_b)
        if h > l: w += h - l
    return w

def sumW(E): return sum(W_a(E, a) for a in range(1, 7))
def measS7(E): return sum((hi-lo) for lo, hi in measS7_arcs(E))
def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))
def Tcond(E): return sum(F(1, abs(e)) for e in E if e != 0)
def res_partition(E):
    """multiset of residues (how many generators per residue)."""
    c = defaultdict(int)
    for e in E:
        if e != 0: c[e % 7] += 1
    c[0] += 1  # the 0 generator sits in residue 0 (sector 0, frozen)
    return tuple(c[r] for r in range(7))

if __name__ == "__main__":
    print("#" * 78)
    print("# FOSTER TOTAL-CONDUCTANCE: is sumW-max == T-max (harmonic-sum)?")
    print("#" * 78)
    C = consec(8)
    print(f"\nconsec={C}  T(consec)=sum 1/|e| = {Tcond(C)} = {float(Tcond(C)):.6f}")
    print(f"   measS7(consec)={measS7(C)}={float(measS7(C)):.6f}, sumW={sumW(C)} "
          f"(== measS7: {sumW(C)==measS7(C)})")

    for W in [12, 14, 16]:
        bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
        allsh = [list(E) for E in bank if primitive(E)]
        full = [E for E in allsh if is_full_residue(E)]
        print(f"\n=== span<= {W}: {len(allsh)} primitive shapes ({len(full)} full-residue) ===")

        # (T1) T-max over ALL shapes is consec (trivial)
        Tmax_shape = max(allsh, key=Tcond)
        print(f"  (T1) argmax T over ALL k=8 shapes: {Tmax_shape} (==consec: {Tmax_shape==C})")

        # (T2) sumW-max over the FULL-residue stratum and over ALL shapes
        sw_full_max = max(full, key=sumW)
        print(f"  (T2) argmax sumW over FULL-residue: {sw_full_max} (==consec: {sw_full_max==C})")
        # over ALL shapes (sumW only meaningful when full cover possible; non-full
        # have sumW could be lower). check measS7 argmax over all:
        ms_all_max = max(allsh, key=measS7)
        print(f"       argmax measS7 over ALL shapes: {ms_all_max} (==consec: {ms_all_max==C})")

        # (T3) rank correlation of sumW vs T over full-residue
        full_s = sorted(full, key=Tcond)
        Ts = [float(Tcond(E)) for E in full_s]
        Ss = [float(sumW(E)) for E in full_s]
        # Spearman-ish: count concordant/discordant pairs sign
        conc = disc = 0
        N = len(full_s)
        # sample to keep it fast at large N
        idx = list(range(N))
        if N > 300:
            import random; random.seed(0); idx = random.sample(idx, 300)
        for ii in range(len(idx)):
            for jj in range(ii+1, len(idx)):
                a, b = idx[ii], idx[jj]
                dt = Ts[a]-Ts[b]; ds = Ss[a]-Ss[b]
                if dt*ds > 0: conc += 1
                elif dt*ds < 0: disc += 1
        tot = conc+disc
        print(f"  (T3) sumW vs T concordance over full-res: {conc}/{tot} concordant "
              f"= {conc/tot:.3f}  (1.0 => sumW monotone in T)")

    # (T4) DECISIVE: restrict to consec-type residue partition; is sumW monotone in T?
    print("\n(T4) CONSEC-TYPE TOPOLOGY (residue partition = consec's): sumW monotone in T?")
    W = 18
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    pc = res_partition(C)
    print(f"   consec residue partition = {pc} (residue 0 has 2 gens incl e=0)")
    same_top = [E for E in full if res_partition(E) == pc]
    print(f"   shapes with consec's residue partition: {len(same_top)}")
    same_top.sort(key=Tcond)
    # monotone test within this topology
    sw_st = {tuple(E): sumW(E) for E in same_top}
    T_st = {tuple(E): Tcond(E) for E in same_top}
    viol = 0; tested = 0
    for i in range(len(same_top)):
        for j in range(len(same_top)):
            if i == j: continue
            E1, E2 = same_top[i], same_top[j]
            if T_st[tuple(E1)] > T_st[tuple(E2)]:
                tested += 1
                if sw_st[tuple(E1)] < sw_st[tuple(E2)] - F(1, 10**12):
                    viol += 1
    print(f"   within consec-topology: T-larger but sumW-smaller pairs: {viol}/{tested} "
          f"(sumW monotone in T on fixed topology: {viol==0})")
    # also: is consec the T-max within its topology? (trivially yes -> sumW-max)
    print(f"   consec is T-max within its topology: {max(same_top, key=Tcond)==C}")
