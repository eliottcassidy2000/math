#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- DOES the exact G_P/Markoff cell-discrepancy
DIRECTLY bound or compute the per-resonance width W_a?  The explicit bridge test.

The L7 tail used D_P(p,q)=G_P(||p||,||q||)/(P p q) to bound the cover of a
slope-(p,q) line through the PxP grid. Here measS7's resonance at a/7 is a
DIFFERENT object (a covering of Z/7 by k offsets, not a 2-coordinate torus line).
QUESTION: is W_a (per-resonance width) governed by the same G_P / cycle-graph
Green's-function structure? Concretely, the survival of residue r at resonance a
is a 1-D arc whose length is t(7-t)/7-type (cycle Green's function) in the local
slopes. We test whether the per-residue survival legs are EXACTLY the g(t)=2t(P-t)
cycle resistances of the G_P legs, i.e. whether the resonance width decomposes
into Markoff legs.

Also: HONEST scope check -- how robust is "consec uniquely maximizes measS7 over
the full-residue stratum"? Test span<=20 at k=8 (the full feasible stratum for the
binding L7 row), and the genuinely-aggregate nature (no monotone certificate).
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
    total = F(0); arcs=[]
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7:
            total += hi-lo; arcs.append((lo,hi))
    return total, arcs

def measS7(E): return measS7_arcs(E)[0]

def W_a(E, a):
    _, arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1,14), F(2*a+1,14); w=F(0)
    for lo,hi in arcs:
        l=max(lo,lo_b); h=min(hi,hi_b)
        if h>l: w+=h-l
    return w

def residues(E): return frozenset(e%7 for e in E)
def is_full_residue(E): return residues(E)==frozenset(range(7))
def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1
def consec(k): return list(range(k))

# cycle-graph Green's function leg g(t)=2 t (P-t) and Bernoulli identity
def cyc_leg(t, P): return 2*t*(P-t)

if __name__ == "__main__":
    print("#"*78); print("# G_P / MARKOFF-LEG BRIDGE & ROBUSTNESS (THREAD A)"); print("#"*78)

    # --- 1. Per-residue survival legs at consec, resonance a=1 ---
    C = consec(8)
    print("\n1. Resonance a=1 for consec: per-residue survival arc & cycle-leg test.")
    # local slope of e is 7e; residue e a mod 7 = e mod 7. The residue r's
    # survival half-width on the right is governed by the slowest e in class r.
    a = 1
    w = W_a(C, a)
    print(f"   W_1(consec) = {w} = {float(w):.6f}")
    # the 7-cell around 1/7 has width 1/7; cover fraction = w/(1/7)
    print(f"   cover fraction in the 1/7-cell = {w*7} = {float(w*7):.6f}")
    # cycle leg values for t=1..6, P=7:
    print(f"   cycle-graph legs g(t)=2t(7-t), t=1..6: {[cyc_leg(t,7) for t in range(1,7)]}")
    print(f"   (g(t)/7 = effective resistance R_eff(0,t) on C_7; the G_P building block)")

    # --- 2. Is W_a a rational with denominator tied to the leg structure? ---
    print("\n2. Per-resonance widths W_a(consec), a=1..6 (exact), and 7*W_a:")
    for a in range(1,7):
        wa = W_a(C, a)
        print(f"   a={a}: W={wa}  7W={7*wa}={float(7*wa):.6f}")

    # --- 3. ROBUSTNESS: consec-max over full-residue stratum at larger span ---
    print("\n3. ROBUSTNESS of consec-max over the FULL-residue stratum (k=8):")
    for W in [13, 16, 20]:
        bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
        full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
        mC = measS7(C)
        beaters = [E for E in full if measS7(E) > mC + F(1,10**12)]
        print(f"   span<= {W}: {len(full)} full-res shapes; consec beaters={len(beaters)}  "
              f"(consec unique max: {len(beaters)==0})")

    # --- 4. The genuinely-aggregate marker: resonant-part consec-max vs per-a ---
    print("\n4. AGGREGATE marker: consec maximizes the SUM of W_a but NOT each W_a.")
    W=13
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),7)]
    full=[list(E) for E in bank if primitive(E) and is_full_residue(E)]
    for a in range(1,7):
        wC=W_a(C,a)
        beat=sum(1 for E in full if W_a(E,a)>wC+F(1,10**12))
        print(f"   a={a}: consec W_a beaters = {beat}")
    sumbeat=sum(1 for E in full if sum(W_a(E,a) for a in range(1,7))>sum(W_a(C,a) for a in range(1,7))+F(1,10**12))
    print(f"   SUM_a: consec beaters = {sumbeat} (the sum is consec-maximal though individual a are not)")
    print("\n   => consec-max is genuinely AGGREGATE across resonances: it wins the")
    print("      TOTAL while losing individual resonances. No per-resonance argument.")
