#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2, THE CRUX: sumW is a FUNCTION of the
per-residue CONDUCTANCE vector c_r = sum_{e: e mod7=r} 1/|e| (verified 0/2260).
This is the electrical invariant the leg-profile missed. Now PIN the extremality.

THE FOSTER/CONDUCTANCE EXTREMALITY HYPOTHESIS (falsifiable):

  H-COND-MONO: sumW(E) = Phi(c_1,...,c_6 ; c_0) is MONOTONE-INCREASING in the
    full conductance vector. Increasing any c_r (adding a faster-or-extra clock
    to residue r, i.e. a SMALLER |e|) can only INCREASE survival. (Electrical:
    adding conductance shrinks resistance; here it lengthens the survival window.)

  H-COND-MAX: consec uniquely MAXIMIZES sumW because it achieves the
    LEXICOGRAPHICALLY / DOMINANT conductance profile achievable with k=8 distinct
    integers including 0 covering all residues: c = (1, 1/2, 1/3, 1/4, 1/5, 1/6)
    for residues 1..6 and c_0=1/7. Any other full-residue shape has, for each
    residue, its single smallest |e| >= the consec value => c_r <= consec's
    => by monotonicity sumW <= sumW(consec). (Caveat: extra clocks in a residue
    ADD conductance, so domination must account for the FULL multiset.)

TESTS:
 (T1) MONOTONICITY of Phi in the conductance vector: build many shapes, sort by
      c-vector domination, and check sumW respects it. Report violations + the
      MINIMAL violating pair (the exact obstruction, if any).
 (T2) The DOMINATION CLAIM: does consec's conductance vector DOMINATE every
      other full-residue shape's c-vector componentwise? (If yes + T1 mono => DONE.)
      Account for parallel clocks: a competitor could ADD a 2nd clock to some
      residue and beat consec's c_r there while paying elsewhere.
 (T3) If neither pure domination nor pure monotonicity holds, test the SCHUR /
      total-conductance budget: is sumW Schur-concave or Schur-convex in c, and
      does consec extremize sum_r c_r? Foster says a TOTAL is conserved; here
      maybe sum_r c_r is NOT conserved but consec maximizes it.
 (T4) THE SHARP FORM: regress sumW against simple symmetric functions of c to
      find the exact Phi (or its leading behavior) -- if Phi is (e.g.) increasing
      and consec dominates, the proof is an electrical Rayleigh-monotonicity.
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
def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def cond_vec(E):
    """c_r = sum_{e in E, e!=0, e mod7=r} 1/|e|. (e=0 contributes 0: never drifts.)"""
    byres = defaultdict(F)
    for e in E:
        if e != 0:
            byres[e % 7] += F(1, abs(e))
    return tuple(byres[r] for r in range(7))

if __name__ == "__main__":
    print("#" * 78)
    print("# CONDUCTANCE EXTREMALITY: is sumW monotone & consec conductance-dominant?")
    print("#" * 78)
    C = consec(8)
    cC = cond_vec(C)
    sC = sumW(C)
    print(f"\nconsec = {C}")
    print(f"  conductance c_r (r=0..6) = {[str(x) for x in cC]}")
    print(f"  sumW(consec) = {sC} = {float(sC):.6f}")

    W = 18
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    print(f"\nfull-residue stratum span<= {W}: {len(full)} shapes")
    cv = {tuple(E): cond_vec(E) for E in full}
    sw = {tuple(E): sumW(E) for E in full}

    # (T1) MONOTONICITY of Phi in the conductance vector
    print("\n(T1) Phi monotone in conductance? (c(E1)>=c(E2) componentwise => sumW(E1)>=sumW(E2))")
    def dom(a, b):  # a >= b componentwise, a != b
        return all(x >= y for x, y in zip(a, b)) and a != b
    pool = full[:600]
    viol = []
    for E1 in pool:
        c1 = cv[tuple(E1)]
        for E2 in pool:
            if E1 is E2: continue
            c2 = cv[tuple(E2)]
            if dom(c1, c2):
                if sw[tuple(E1)] < sw[tuple(E2)] - F(1, 10**12):
                    viol.append((E1, E2, sw[tuple(E1)], sw[tuple(E2)]))
    print(f"   conductance-domination pairs tested in pool of {len(pool)}; "
          f"monotonicity violations: {len(viol)}")
    for E1, E2, v1, v2 in viol[:5]:
        print(f"      c-dom {E1} (c={[str(x) for x in cv[tuple(E1)]]}) sumW {float(v1):.6f}")
        print(f"           >= {E2} (c={[str(x) for x in cv[tuple(E2)]]}) sumW {float(v2):.6f}  VIOL")

    # (T2) does consec's conductance dominate every other full-residue shape?
    print("\n(T2) consec conductance-dominant over all full-residue shapes?")
    nondom = [E for E in full if any(cv[tuple(E)][r] > cC[r] for r in range(7))]
    print(f"   shapes exceeding consec's c_r in some residue: {len(nondom)} "
          f"(consec dominant: {len(nondom)==0})")
    # show a few that exceed, and whether they still have lower sumW
    for E in nondom[:6]:
        higher = [r for r in range(7) if cv[tuple(E)][r] > cC[r]]
        print(f"      {E}: higher c in residues {higher}; sumW {float(sw[tuple(E)]):.6f} "
              f"(< consec {float(sC):.6f}: {sw[tuple(E)] < sC})")

    # (T3) total conductance budget: does consec maximize sum_r c_r?
    print("\n(T3) TOTAL conductance T(E) = sum_r c_r: does consec maximize it?")
    def T(E): return sum(cv[tuple(E)])
    TC = T(C)
    bigger = [E for E in full if T(E) > TC]
    print(f"   T(consec) = {TC} = {float(TC):.6f}; shapes with HIGHER T: {len(bigger)}")
    for E in sorted(bigger, key=lambda E: -T(E))[:5]:
        print(f"      {E}: T={float(T(E)):.6f} > consec, but sumW={float(sw[tuple(E)]):.6f} "
              f"(< consec sumW {float(sC):.6f}: {sw[tuple(E)] < sC})")

    # (T4) Schur direction: is sumW Schur-convex or -concave in c at fixed T?
    # Compare consec (c spread out: 1,1/2,...,1/7) vs shapes with similar T but
    # more even c. Test: among shapes, does MORE SPREAD c => HIGHER sumW (Schur-convex)?
    print("\n(T4) Schur direction: sumW vs conductance spread (variance of c) at comparable T:")
    # bucket by T rounded, compare sumW vs Var(c)
    import statistics
    def var_c(E):
        cs = [float(x) for x in cv[tuple(E)]]
        return statistics.pvariance(cs)
    # correlation sign of (var_c, sumW) within the stratum
    xs = [var_c(E) for E in full]; ys = [float(sw[tuple(E)]) for E in full]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    cov = sum((x-mx)*(y-my) for x, y in zip(xs, ys))/len(xs)
    print(f"   cov(Var(c), sumW) = {cov:.6f}  -> sign {'+' if cov>0 else '-'} "
          f"({'spread helps (Schur-convex)' if cov>0 else 'spread hurts (Schur-concave)'})")
    print(f"   consec Var(c) = {var_c(C):.6f} (extreme spread: 1 down to 1/7)")
