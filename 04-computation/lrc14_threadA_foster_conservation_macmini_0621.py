#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2 FOLLOW-UP: THE FOSTER CONSERVATION LAW.

DISCOVERY (lrc14_threadA_foster_resistance): for EVERY full-residue shape, the
residue legs (one min|e| per residue r in Z/7) satisfy
   sum_{r in Z/7} g(leg_r mod 7) = sum_{t=0}^{6} 2 t (7-t) = 112   (CONSTANT).
This is a FOSTER-TYPE CONSERVATION: the legs form a transversal of Z/7, so the
sum of cycle-graph resistances g(residue) is fixed at 112 = 14 * Kirchhoff-ish.
It is the EXACT analogue of Foster's theorem (sum of edge resistances = const).

So the residue structure gives a CONSERVED total resistance. The real question:
given this conservation, why does consec MAXIMIZE sum_a W_a? The answer must be:
sum_a W_a is a CONVEX/SCHUR functional of the actual MAGNITUDES (not just residues)
under the conserved residue-resistance budget, and consec puts ALL magnitudes at
their MINIMUM (the contracted network).

This script tests a sharper, FALSIFIABLE electrical hypothesis:

H-FOSTER-2:  W_a is, EXACTLY, governed by the SLOWEST clock per sector, and the
  total survival sum_a W_a is monotone-DECREASING in EACH magnitude |e| (holding
  the residue structure / network topology fixed). I.e. the AP minimizes every
  leg magnitude => maximizes survival. The obstruction to a pure leg-profile proof
  (F2 failed) is that W_a also depends on the SECOND-slowest clock per sector
  (the network has parallel resistors), so we need the FULL conductance vector,
  not just the min. Test:
   (A) monotone-coordinatewise: is sum_a W_a monotone-decreasing in each |e_i|
       when all other generators fixed? (the electrical-contraction principle)
   (B) the CONSERVATION 112 holds universally (verify on big stratum).
   (C) PARALLEL-RESISTOR refinement: define per-residue conductance
       c_r = sum_{e: e mod 7 == r} 1/(7|e|)  (parallel clocks add conductance),
       survival ~ 1/c_r-type. Does sum_a W_a track a function of the FULL c_r
       vector (not just min)? Test if (c_r) determines sum_a W_a (sharper than legs).
   (D) consec as the MAXIMAL-CONDUCTANCE placement: consec has the smallest
       magnitudes => largest conductances c_r => is it the unique maximizer of a
       concave conductance functional under the residue-budget?
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
def g(t, N=7): return 2 * t * (N - t)

if __name__ == "__main__":
    print("#" * 78)
    print("# FOSTER CONSERVATION LAW + CONTRACTION/CONDUCTANCE TEST (ANGLE 2)")
    print("#" * 78)
    C = consec(8)

    # ---- (B) Conservation: sum_r g(leg_r mod7) = 112 universal ----
    print("\n(B) FOSTER CONSERVATION: sum_{r in Z/7} g(leg_r mod 7) = 112 for ALL full-res?")
    W = 18
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    def legs(E):
        byres = defaultdict(list)
        for e in E: byres[e % 7].append(abs(e))
        return {r: min(byres[r]) for r in range(7)}
    consts = set()
    for E in full:
        L = legs(E)
        consts.add(sum(g(L[r] % 7) for r in range(7)))
    # NOTE: leg_r mod 7 == r always (since leg is in residue class r), so g(leg_r mod7)=g(r)
    print(f"   distinct values of sum_r g(leg_r mod 7) over {len(full)} shapes: {consts}")
    print(f"   sum_{{t=0..6}} g(t) = {sum(g(t) for t in range(7))} "
          f"(this is the conserved Foster total; leg_r mod 7 = r always => trivially const)")
    print("   -> So the RESIDUE-resistance is conserved by construction (transversal).")
    print("      The DEGREE OF FREEDOM is the actual magnitude per residue, not residue.")

    # ---- (A) CONTRACTION: is sumW monotone-decreasing in each |e|? ----
    # Take consec and increase one generator e_i by +7 (stays full-residue, same
    # residue r) and check sumW decreases. Then a broader random test.
    print("\n(A) CONTRACTION PRINCIPLE: increasing one magnitude (same residue) lowers sumW?")
    base = sumW(C)
    print(f"   consec sumW = {float(base):.6f}")
    viol = 0; tested = 0
    # systematic: for each non-zero generator, bump by +7 (preserve residue/full)
    for i in range(1, 8):
        E2 = C.copy()
        E2[i] = E2[i] + 7  # increase magnitude, same residue
        if not is_full_residue(E2) or not primitive(E2): continue
        v = sumW(sorted(E2)); tested += 1
        dec = v < base
        if not dec: viol += 1
        print(f"   bump e={C[i]}->{E2[i]} (res {C[i]%7}): sumW {float(v):.6f}  "
              f"decreased? {dec}")
    print(f"   single-bump from consec: {tested} tests, {viol} violations of contraction")

    # broader: among full-res shapes, test pairwise contraction: if E2 dominates E1
    # in magnitudes (same residue assignment, each |e| >=, at least one >), is
    # sumW(E2) <= sumW(E1)? This is the electrical monotonicity (Rayleigh).
    print("\n   Rayleigh-monotonicity over full-res stratum (residue-matched domination):")
    # group by residue->position assignment is complex; instead test: for shapes
    # with identical residue partition pattern, magnitude-domination => sumW-anti.
    def res_sorted_mags(E):
        byres = defaultdict(list)
        for e in E: byres[e % 7].append(abs(e))
        return tuple(tuple(sorted(byres[r])) for r in range(7))
    pool = full[:400]  # cap for time
    ms = {tuple(E): sumW(E) for E in pool}
    rm = {tuple(E): res_sorted_mags(E) for E in pool}
    def dom(a, b):
        # a dominates b: same shape (count per residue), each entry a>=b
        if any(len(x) != len(y) for x, y in zip(a, b)): return False
        return all(all(xi >= yi for xi, yi in zip(x, y)) for x, y in zip(a, b)) and a != b
    rviol = 0; rtest = 0
    for E1 in pool:
        for E2 in pool:
            if E1 is E2: continue
            if dom(rm[tuple(E1)], rm[tuple(E2)]):
                rtest += 1
                # E1 magnitudes >= E2 => expect sumW(E1) <= sumW(E2)
                if ms[tuple(E1)] > ms[tuple(E2)] + F(1, 10**12):
                    rviol += 1
    print(f"   residue-matched magnitude-domination pairs: {rtest}, "
          f"Rayleigh violations: {rviol}  (monotone contraction: {rviol==0})")

    # ---- (C) CONDUCTANCE refinement: does the full c_r vector determine sumW? ----
    print("\n(C) PARALLEL-CONDUCTANCE: c_r = sum_{e: e mod7=r} 1/|e|; does (c_r) -> sumW?")
    def cond_vec(E):
        byres = defaultdict(F)
        for e in E:
            if e == 0:
                byres[0] += F(0)  # 0 has infinite period; treat as residue-0 base
            else:
                byres[e % 7] += F(1, abs(e))
        return tuple(byres[r] for r in range(7))
    bycond = defaultdict(set)
    for E in full:
        bycond[cond_vec(E)].add(sumW(E))
    noncond = sum(1 for v in bycond.values() if len(v) > 1)
    print(f"   conductance vectors with >1 distinct sumW: {noncond}/{len(bycond)} "
          f"-> (c_r) determines sumW: {noncond==0}")

    # ---- (D) consec = max-conductance: largest c_r per residue ----
    print("\n(D) consec maximizes per-residue conductance? (smallest magnitudes)")
    cC = cond_vec(C)
    print(f"   consec conductance vector c_r = {[str(x) for x in cC]}")
    # consec should dominate: every other shape has c_r <= consec's for each r
    notdom = 0
    for E in full:
        cE = cond_vec(E)
        if any(cE[r] > cC[r] for r in range(1, 7)):  # residue 0 special
            notdom += 1
    print(f"   shapes with HIGHER conductance than consec in some residue: {notdom} "
          f"(consec is conductance-maximal: {notdom==0})")
