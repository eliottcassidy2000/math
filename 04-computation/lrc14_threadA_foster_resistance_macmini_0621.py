#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2: EFFECTIVE-RESISTANCE / FOSTER MODEL
of the LRC(14) wall. Goal: express sum_a W_a as an electrical-network quantity
on the cycle graph C_7 and test whether Foster's theorem / a network-extremal
principle forces the consecutive (AP) assignment.

----------------------------------------------------------------------------
THE ELECTRICAL DICTIONARY (what we are testing)
----------------------------------------------------------------------------
g(t) = 2 t (7-t) is, up to scale, the effective resistance R_eff(0,t) on the
cycle graph C_7: on C_N with unit resistors, R_eff(i,j) = d(N-d)/N where
d = |i-j|. So R_eff(0,t) on C_7 = t(7-t)/7, and g(t) = 2*7*R_eff(0,t).

The per-residue SURVIVAL LEG at resonance a:
  Near x=a/7, clock e sits at sector (e*a mod 7) and drifts at SPEED 7e.
  Residue r in Z/7 is "covered" while at least one e with sector r is present.
  The slowest representative (min |e| in the class that maps to r) holds r
  longest. Survival half-width contributed ~ proportional to 1/(7 * speed) of
  the *fastest* clock that has to be relied upon, i.e. governed by min|e| legs.

HYPOTHESIS H-FOSTER (falsifiable):
  Define a network N_a on vertex set Z/7 (the 7 sectors). For the AP/consec
  assignment, the survival width W_a is an effective-resistance-type quantity,
  and the consec assignment is the unique RESISTOR PLACEMENT that, summed over
  all 7 (or 6) resonances, maximizes a total-resistance / total-coverage
  functional subject to a Foster-type CONSERVATION law sum = const.

CONCRETE TESTS:
  (T1) Is g(t)/7 EXACTLY R_eff(0,t) on C_7?  (sanity of the dictionary)
  (T2) FOSTER SUM RULE: on C_7, sum over the 7 cycle edges of R_eff(edge) = N-1
       = 6. Does any analogous CONSERVATION hold for sum_a W_a or for the legs?
  (T3) Express W_a (consec) via the C_7 resistances and see if 7*W_a equals a
       combination of R_eff values / Foster-type partial sums.
  (T4) The KEY extremal test: model sum_a W_a as a function of a "conductance
       placement" c_r (per residue), and ask whether consec's placement is the
       extremizer of a concave/Schur functional. Foster => the SUM of a resistance
       functional is conserved; if W_a-sum is a SCHUR-CONCAVE functional of the
       residue legs that has consec at the boundary of the conservation simplex.
----------------------------------------------------------------------------
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ------------------------------------------------------------------ measS7 / W_a
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

def measS7(E):
    return sum((hi-lo) for lo, hi in measS7_arcs(E))

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

# ------------------------------------------------------------------ C_7 resistance
def R_eff_cycle(i, j, N=7):
    """Effective resistance between i,j on cycle C_N (unit resistors)."""
    d = abs(i - j) % N
    d = min(d, N - d)
    return F(d * (N - d), N)

def g(t, N=7): return 2 * t * (N - t)

# ------------------------------------------------------------------ MAIN
if __name__ == "__main__":
    print("#" * 78)
    print("# FOSTER / EFFECTIVE-RESISTANCE MODEL OF THE LRC WALL (ANGLE 2)")
    print("#" * 78)

    N = 7
    C = consec(8)

    # (T1) dictionary check: g(t)/7 == R_eff(0,t) on C_7
    print("\n(T1) g(t)=2t(7-t) vs C_7 effective resistance R_eff(0,t)=t(7-t)/7:")
    for t in range(1, 7):
        print(f"   t={t}: g(t)={g(t)}, g(t)/14={F(g(t),14)}, R_eff(0,t)={R_eff_cycle(0,t)}, "
              f"g(t)/14 == R_eff? {F(g(t),14)==R_eff_cycle(0,t)}")

    # (T2) Foster sum rule on C_7
    print("\n(T2) FOSTER sum rule on C_7: sum over 7 cycle EDGES of R_eff(edge) = N-1 = 6?")
    fsum = sum(R_eff_cycle(i, (i+1) % N) for i in range(N))
    print(f"   sum_edges R_eff = {fsum} (expect {N-1})  Foster holds: {fsum == N-1}")
    # also the all-pairs sum (Kirchhoff index)
    Kf = sum(R_eff_cycle(i, j) for i in range(N) for j in range(i+1, N))
    print(f"   Kirchhoff index Kf(C_7) = sum_{{i<j}} R_eff = {Kf} = {float(Kf):.4f} "
          f"(closed form N(N^2-1)/12 = {F(N*(N*N-1),12)})")

    # (T3) express 7*W_a (consec) and compare to C_7 resistance combinations
    print("\n(T3) consec 7*W_a per resonance vs C_7 resistance values:")
    for a in range(1, 7):
        wa = W_a(C, a)
        print(f"   a={a}: 7*W_a={7*wa}={float(7*wa):.5f}")
    tot = sum(W_a(C, a) for a in range(1, 7))
    print(f"   sum_a W_a (consec) = {tot} = {float(tot):.6f}")
    print(f"   measS7(consec)     = {measS7(C)} = {float(measS7(C)):.6f}  (== sum_a? {tot==measS7(C)})")

    # (T4) THE CORE TEST: is sum_a W_a a Schur-CONCAVE functional of the residue
    # LEGS (min|e| per residue), with a Foster-type CONSERVATION constraint that
    # pins consec at the extremum? We model: legs L = (L_0,...,L_6), and survival
    # contribution of residue r ~ a function of L_r. Foster's law would say a
    # total flux is conserved. Test the conjecture:
    #   F1: sum_a W_a is a function only of the MULTISET of legs (sorted profile).
    #   F2: among full-residue shapes with the SAME sorted leg profile, sum_a W_a
    #       is constant (testing if legs alone determine the electrical quantity).
    #   F3: the consec leg-profile (0,1,2,3,4,5,6) is the UNIQUE Foster-minimizer
    #       of sum_r R_eff-type leg weights.
    print("\n(T4) Schur/Foster test on the residue LEG profile:")
    def leg_profile(E):
        byres = defaultdict(list)
        for e in E: byres[e % 7].append(abs(e))
        if not all(r in byres for r in range(7)): return None
        return tuple(sorted(min(byres[r]) for r in range(7)))

    W = 16
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    print(f"   full-residue stratum span<= {W}: {len(full)} shapes")

    byprof = defaultdict(list)
    for E in full:
        p = leg_profile(E)
        byprof[p].append((E, sum(W_a(E, a) for a in range(1, 7))))
    nonconst = 0
    for p, lst in byprof.items():
        vals = set(v for _, v in lst)
        if len(vals) > 1:
            nonconst += 1
    print(f"   F2: profiles where sum_a W_a is NOT constant: {nonconst}/{len(byprof)} "
          f"-> legs alone determine the sum: {nonconst==0}")

    # Foster-minimizer test: among profiles, consec=(0,1,2,3,4,5,6). Build a
    # candidate Foster weight W_foster(profile) = sum_r R_eff(0, leg_r mod 7)-like
    # and see whether MAXIMIZING sum_a W_a corresponds to MINIMIZING a leg-energy.
    pC = leg_profile(C)
    print(f"   consec leg profile = {pC}")
    # leg-energy candidate: sum_r g(leg_r mod 7)? or sum_r leg_r? Try several.
    def E_sumleg(p):  return sum(p)
    def E_sumg(p):    return sum(g(t % 7) for t in p)
    def E_sumsq(p):   return sum(t*t for t in p)
    best_by_sum = max(full, key=lambda E: sum(W_a(E, a) for a in range(1, 7)))
    print(f"   shape maximizing sum_a W_a = {best_by_sum} "
          f"(== consec? {best_by_sum==C})")
    # rank profiles by sum_a W_a (max representative) and show their energies
    prof_rep = {}
    for p, lst in byprof.items():
        prof_rep[p] = max(v for _, v in lst)
    ranked = sorted(prof_rep.items(), key=lambda kv: -kv[1])[:8]
    print("\n   Top-8 leg profiles by max sum_a W_a, with leg-energies:")
    print("   profile                        | sum_aW_a  | sumleg | sum g(t) | sum t^2")
    for p, v in ranked:
        print(f"   {str(p):30s} | {float(v):.6f} | {E_sumleg(p):6d} | {E_sumg(p):8d} | {E_sumsq(p):7d}")
