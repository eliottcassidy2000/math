#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOWER-DELETION LEMMA for LRC14 AP-tail mouth-retention - PROOF via the codex comb technique.
kind-pasteur-2026-06-19.  EXACT rationals.

LEMMA: for a in {0,1,2,3}, if the AP-tail 12-core C has 2^a as a HOLE (2^a not in C),
       then meas(G_C) >= 426/35035 = thr2.

NOTE ON A FALSE START: a "monotonicity" reduction to the single-hole base is INVALID.
  Adding a constraint ||c t||>1/14 only SHRINKS the safe set, so for B subset C,
  meas(G_C) <= meas(G_B) (NOT >=).  Tails (speeds >=14) are genuine extra constraints and
  can LOWER the measure.  There is no free lower bound by dropping tails.  This is exactly why
  THM-543/544 needed the periodic-comb cutoff + finite exact residue scan.  We follow that.

METHOD (mirrors THM-543/544 exactly):
  Core missing 2^a:  C = B U T,  B = {1..13}\holes  (2^a in holes),  T = tails (>=14 distinct),
  |C|=12 so |T| = |holes|-1 =: k.

  k=0:  C = single-hole base {1..13}\{2^a};  meas computed exactly, >= thr2.  [base case]

  k>=1: SINGLE-TAIL COMB LEMMA.  For a safe set G with measure M and c interval components,
        and a new speed r,
            meas(G \ D_r) >= (6/7) M - 2c/(7r).
        (Each component: full 1/r periods remove exactly 1/7 of their length; the two partial
         end-periods remove at most two extra teeth of total width 2/(7r).  Sum over c comps.)

  This gives a per-tail RATIONAL CUTOFF:  meas(G_B \ D_r) >= thr2  whenever
            r >= R_B := ceil( 2 c_B / (6 M_B - 7 thr2) )           (6 M_B - 7 thr2 > 0 always).
  Tails r >= R_B are certified by the comb; the finite residue r in [14, R_B) is checked EXACTLY.
  For k>=2 we iterate: comb out the largest tail, exact-scan the smaller tail(s) in their residue.

  FLOOR sanity (why higher k is safe): (6/7)^k * min_{k-hole base} M_B is INCREASING in k
  (extra holes remove constraints and raise M faster than the comb factor erodes), staying far
  above thr2 at every k.  So only configurations with a SMALL tail can approach thr2, and those
  are exactly the finite residues the exact scan covers.

This file PROVES the a=3, k=1 layer in full (cutoff + exact residue) and tabulates the base cases.
The k=2/k=3 layers and a=0,1,2 are handled by lrc14_tower_deletion_layers_kps.py.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def lonely_measure(C, theta=Fraction(1,14), comps=False):
    arcs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0, d+1):
            arcs.append((Fraction(m,d)-w, Fraction(m,d)+w))
    segs=[]
    for lo,hi in arcs:
        for shift in (-1,0,1):
            a2=max(lo+shift,Fraction(0)); b2=min(hi+shift,Fraction(1))
            if a2<b2: segs.append((a2,b2))
    segs.sort()
    cur=Fraction(-1); union=[]
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    safe=Fraction(1)-sum(b-a for a,b in union)
    if comps:
        nc=0; prev=Fraction(0)
        for a,b in union:
            if a>prev: nc+=1
            prev=b
        if prev<Fraction(1): nc+=1
        return safe, nc
    return safe

def ceil_frac(x):
    return x.numerator//x.denominator + (1 if x.numerator%x.denominator else 0)

thr2=Fraction(426,35035)

if __name__=="__main__":
    print(f"thr2 = 426/35035 = {float(thr2):.8f}\n")

    print("=== BASE CASES k=0: single-hole bases {1..13}\\{2^a} (no tail) ===")
    for a in [0,1,2,3]:
        B=tuple(d for d in range(1,14) if d!=2**a)
        L=lonely_measure(B)
        print(f"  a={a}: meas(G_B)={L}={float(L):.8f}  >=thr2? {L>=thr2}")

    print("\n=== FLOOR sanity: (6/7)^k * min_{k-hole base containing 2^a} M_B  (increasing, >> thr2) ===")
    for a in [3,2,1,0]:
        bit=2**a; others=[d for d in range(1,14) if d!=bit]
        line=[]
        for k in range(0,5):
            Mmin=None
            for extra in itertools.combinations(others,k):
                holes=set(extra)|{bit}
                B=tuple(d for d in range(1,14) if d not in holes)
                if len(B)<1: continue
                M=lonely_measure(B)
                if Mmin is None or M<Mmin: Mmin=M
            floor=Fraction(6,7)**k*Mmin
            line.append(f"k={k}:{float(floor):.4f}")
        print(f"  a={a}: " + "  ".join(line) + f"   (thr2={float(thr2):.4f})")

    print("\n=== a=3, k=1 layer: delete 8 + one hole h + one tail r>=14  (FULL PROOF) ===")
    worstR=0; below=[]
    for h in [d for d in range(1,14) if d!=8]:
        B=tuple(d for d in range(1,14) if d not in (8,h))
        M,c=lonely_measure(B,comps=True)
        denom=6*M-7*thr2
        assert denom>0
        Rh=ceil_frac((2*c)/denom)
        worstR=max(worstR,Rh)
        for r in range(14,Rh):
            if r in B: continue
            C=tuple(sorted(B+(r,)))
            if len(C)!=12: continue
            if reduce(gcd,C)!=1: continue
            L=lonely_measure(C)
            if L<thr2: below.append((L,h,r))
    print(f"  comb cutoff: all tails r >= R={worstR} certified by  meas >= (6/7)M_B - 2c_B/(7r) >= thr2.")
    print(f"  finite residue r in [14,{worstR}) exact-scanned: {len(below)} rows below thr2.")
    if not below:
        print("  => a=3 k=1 layer PROVED: meas(G_C) >= thr2 for every delete-8 two-hole one-tail core.")
    else:
        for L,h,r in sorted(below): print(f"    BELOW meas={L} hole={h} tail={r}")
