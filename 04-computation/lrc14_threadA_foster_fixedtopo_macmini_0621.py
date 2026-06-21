#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2, DECISIVE FIXED-TOPOLOGY TEST.

Global concordance of sumW vs T(harmonic sum) is only ~0.5 (T-max coincides with
sumW-max at consec, but the FUNCTIONS don't track). So a naive "T-max => sumW-max"
is NOT a proof. The electrical-Rayleigh idea survives only if monotonicity holds
WITHIN A FIXED NETWORK TOPOLOGY (fixed residue partition = fixed parallel/series
structure). This is the honest Rayleigh statement: decreasing a single resistance
(= increasing 1/|e|) on a FIXED graph can only increase total conductance/flow.

This script isolates the consec-type topology (residue partition (2,1,1,1,1,1,1):
residue 0 doubled by e=0 & one more, residues 1..6 each once) and asks:
  Within this fixed topology, is sumW monotone-INCREASING in each 1/|e|?
  Equivalently: if we decrease one magnitude (raise its conductance), holding the
  others fixed, does sumW go up? (Pure Rayleigh.)

If YES, then: consec is the unique element of its topology with ALL magnitudes
minimal => unique max of sumW WITHIN the topology. The remaining gap (proving the
consec topology beats OTHER topologies) is the residue-budget/Foster part.

TESTS:
 (T1) consec-type topology, exhaustive coordinatewise Rayleigh: for every shape in
      the topology and every coordinate, decreasing that magnitude by the minimal
      step (to the next-smaller integer of the same residue, if it keeps the shape
      valid) increases sumW. Report violations.
 (T2) FULL Rayleigh within topology: for all pairs E1 (all mags <= E2's,
      residue-matched), is sumW(E1) >= sumW(E2)? (the global monotone form)
 (T3) Do OTHER topologies ever beat consec topology's max? (the budget gap)
      List the top sumW per topology.
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
def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def res_partition(E):
    c = defaultdict(int)
    for e in E:
        if e != 0: c[e % 7] += 1
    c[0] += 1
    return tuple(c[r] for r in range(7))

def res_mags(E):
    """sorted magnitudes per residue (residue 0 includes the e=0 -> treat e=0 mag as 0)."""
    byres = defaultdict(list)
    for e in E:
        byres[e % 7].append(abs(e))
    return tuple(tuple(sorted(byres[r])) for r in range(7))

if __name__ == "__main__":
    print("#" * 78)
    print("# FIXED-TOPOLOGY RAYLEIGH MONOTONICITY (ANGLE 2 DECISIVE TEST)")
    print("#" * 78)
    C = consec(8)
    pc = res_partition(C)
    print(f"\nconsec={C}, residue partition (topology) = {pc}")
    print(f"   sumW(consec) = {float(sumW(C)):.6f}")

    W = 20
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    same = [E for E in full if res_partition(E) == pc]
    print(f"\nfull-residue span<= {W}: {len(full)} shapes; consec-topology: {len(same)} shapes")

    sw = {tuple(E): sumW(E) for E in same}
    rm = {tuple(E): res_mags(E) for E in same}

    # (T2) FULL Rayleigh within topology: E1 mags <= E2 mags (residue-matched) => sumW(E1)>=sumW(E2)
    def le(a, b):  # a <= b componentwise across the per-residue sorted tuples
        return all(all(xi <= yi for xi, yi in zip(x, y)) for x, y in zip(a, b))
    viol = []
    for E1 in same:
        for E2 in same:
            if E1 is E2: continue
            if le(rm[tuple(E1)], rm[tuple(E2)]) and rm[tuple(E1)] != rm[tuple(E2)]:
                if sw[tuple(E1)] < sw[tuple(E2)] - F(1, 10**12):
                    viol.append((E1, E2))
    print(f"\n(T2) FULL Rayleigh in consec-topology: "
          f"mag-domination pairs with sumW INVERTED: {len(viol)}")
    for E1, E2 in viol[:6]:
        print(f"   {E1} (sumW {float(sw[tuple(E1)]):.6f}) mags<= "
              f"{E2} (sumW {float(sw[tuple(E2)]):.6f})  VIOL")
    print(f"   => sumW Rayleigh-monotone within consec topology: {len(viol)==0}")

    # (T1) coordinatewise: decreasing one magnitude (next valid smaller, same residue)
    # raises sumW. Restrict to checking consec is the in-topology max:
    in_topo_max = max(same, key=sumW)
    print(f"\n(T1) in-topology argmax sumW = {in_topo_max} (==consec: {in_topo_max==C})")

    # (T3) per-topology max: does any other topology beat consec topology?
    print("\n(T3) Top sumW per residue topology (does any beat consec topology?):")
    bytopo = defaultdict(list)
    for E in full: bytopo[res_partition(E)].append(E)
    rows = []
    for topo, Es in bytopo.items():
        best = max(Es, key=sumW)
        rows.append((topo, best, sumW(best)))
    rows.sort(key=lambda r: -r[2])
    for topo, best, v in rows[:10]:
        mark = "  <== CONSEC TOPOLOGY" if topo == pc else ""
        print(f"   topo {topo}: max sumW {float(v):.6f} by {best}{mark}")
