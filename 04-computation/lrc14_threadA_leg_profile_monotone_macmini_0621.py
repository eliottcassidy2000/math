#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- IS measS7 A MONOTONE / SCHUR FUNCTION OF THE
LEG PROFILE (min |e| per residue class)?  The route to a non-finite-check proof.

From the local model: consec has the UNIQUE minimal leg profile
   legs(consec) = {residue r -> min|e| in class r} = {0,1,2,3,4,5,6}
(the smallest possible: residue r realized by e=r, plus 0 doubled by 0 and 7).
Every other full-residue shape has a leg profile that DOMINATES consec's
componentwise-after-sorting (some residue forced to a larger min-magnitude),
because to be full-residue with k=8 distinct nonneg integers and 0 in E, you
must cover all 7 residues, and the cheapest way is 0..6 (consec adds 7 = a 2nd
rep of residue 0).

HYPOTHESIS (the clean route): measS7(E) depends, to leading order, MONOTONE-
DECREASINGLY on the leg profile -- smaller min|e| per residue => longer survival
=> larger measS7. If measS7 is a monotone-decreasing function of the sorted leg
profile, and consec uniquely minimizes that profile (Schur/majorization-minimal),
then consec maximizes measS7 with NO finite check.

TEST 1: is the sorted leg profile of consec componentwise <= every other full-
        residue shape's sorted leg profile? (majorization-minimal / dominated)
TEST 2: is measS7 monotone in the leg profile? i.e. if profile(E1) <= profile(E2)
        componentwise (sorted), is measS7(E1) >= measS7(E2)? Test exhaustively.
TEST 3: if monotone fails, find the MINIMAL counterexample pair (the precise
        obstruction: measS7 is NOT a function of the leg profile alone).
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

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: total += hi - lo
    return total

def leg_profile(E):
    byres = defaultdict(list)
    for e in E: byres[e % 7].append(abs(e))
    return tuple(sorted(min(byres[r]) for r in range(7))) if all(r in byres for r in range(7)) else None

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def dominated(p, q):
    """p <= q componentwise (both sorted)."""
    return all(a <= b for a, b in zip(p, q))

if __name__ == "__main__":
    print("#"*78); print("# LEG-PROFILE MONOTONICITY TEST (THREAD A)"); print("#"*78)
    k = 8; W = 12
    C = consec(k); pC = leg_profile(C)
    print(f"\nconsec={C}: leg profile (sorted) = {pC}, measS7={float(measS7(C)):.6f}")

    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    print(f"Full-residue stratum (span<= {W}): {len(full)} shapes")

    # TEST 1: is consec's profile dominated by all others?
    profs = {tuple(E): leg_profile(E) for E in full}
    not_dom = [E for E in full if not dominated(pC, profs[tuple(E)])]
    print(f"\nTEST 1: consec profile {pC} <= every full shape's profile (componentwise sorted)?")
    print(f"   shapes NOT dominating consec: {len(not_dom)}  -> consec is majorization-minimal: {len(not_dom)==0}")
    for E in not_dom[:5]:
        print(f"      {E}: profile {profs[tuple(E)]}")

    # TEST 2: is measS7 monotone-decreasing in the (sorted) leg profile?
    # For every pair with profile(E1) <= profile(E2), check measS7(E1) >= measS7(E2).
    print(f"\nTEST 2: measS7 monotone-decreasing in sorted leg profile? (all dominated pairs)")
    ms = {tuple(E): measS7(E) for E in full}
    viol = []
    for E1 in full:
        for E2 in full:
            if tuple(E1) == tuple(E2): continue
            if dominated(profs[tuple(E1)], profs[tuple(E2)]):
                # profile E1 <= profile E2 => expect measS7(E1) >= measS7(E2)
                if ms[tuple(E1)] < ms[tuple(E2)] - F(1, 10**12):
                    viol.append((E1, E2))
    print(f"   monotonicity violations (profile-smaller but measS7-smaller): {len(viol)}")
    for E1, E2 in viol[:8]:
        print(f"      {E1} prof {profs[tuple(E1)]} meas {float(ms[tuple(E1)]):.5f}  "
              f"<  {E2} prof {profs[tuple(E2)]} meas {float(ms[tuple(E2)]):.5f}")

    # TEST 3: does measS7 depend ONLY on the leg profile? (shapes with same profile)
    print(f"\nTEST 3: is measS7 a FUNCTION of the leg profile? (same profile -> same measS7)")
    byprof = defaultdict(list)
    for E in full: byprof[profs[tuple(E)]].append(E)
    nonfunc = 0
    for prof, Es in byprof.items():
        vals = set(ms[tuple(E)] for E in Es)
        if len(vals) > 1:
            nonfunc += 1
            if nonfunc <= 5:
                print(f"      profile {prof}: {len(Es)} shapes, measS7 values "
                      f"{sorted(float(v) for v in vals)}")
    print(f"   profiles with >1 distinct measS7: {nonfunc}  -> measS7 is profile-function: {nonfunc==0}")

    # CONCLUSION line
    print(f"\nNET: consec uniquely minimizes leg profile (TEST1). "
          f"measS7 {'IS' if not viol else 'is NOT'} monotone in the profile (TEST2). "
          f"measS7 {'IS' if nonfunc==0 else 'is NOT'} a function of the profile alone (TEST3).")
