#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE A part 5: exact sector-symmetry of consec's avoid measures, and a
last coupling probe -- the reflection x -> 1-x and the sector pairing j<->7-j.

Observed (k=8): per-sector avoid J({j}) for j=1..6 was
  [799/2940,179/735,117/490,179/735,799/2940,82/245]
which is PALINDROMIC in j=1..5 (J1=J5, J2=J4) with J6 special, J3 center.
Conjecture: J({j},consec)=J({7-j},consec) for j=1..6 (so J1=J6? no -> check the
right involution). Actually under x->1-x, frac(e(1-x))=frac(-ex), sector j -> 6-j
(since the circle reflects). So J({j}) = J({6-j}) for the reflection that fixes
sector 0. Test the EXACT identity J({j})=J({6-j}) for consec, j=1..6 (6-j in 0..5,
sector 0 excluded -> j and 6-j both in 1..5 pair up; j=6 -> 6-j=0 excluded).

Then: is p_0(E) <= p_0(consec) provable by a 2-element MOON-TYPE switching /
the three-gap (Steinhaus) theorem for the AP orbit? We at least pin the symmetry.
kind-pasteur-2026-06-19 ANGLE-A.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def frac(v): return v-(v.numerator//v.denominator)
def sect(v): return (v.numerator*7)//v.denominator

def avoid_sector_meas(E,j):
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        if all(sect(frac(e*mid))!=j for e in E): tot+=(hi-lo)
    return tot

def avoid_pair_meas(E,A):
    """meas{x: orbit avoids EVERY sector in set A}."""
    E=sorted(set(E)); A=set(A)
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        secs=set(sect(frac(e*mid)) for e in E)
        if not (secs & A): tot+=(hi-lo)
    return tot

if __name__=="__main__":
    for k in [8,9]:
        print(f"\n{'='*70}\nk={k}\n{'='*70}")
        C=list(range(k))
        J=[None]+[avoid_sector_meas(C,j) for j in range(1,7)]  # J[1..6]
        print(f"single-sector avoid J[j], j=1..6:")
        for j in range(1,7): print(f"   J[{j}] = {J[j]} = {float(J[j]):.6f}")
        print(f"reflection identity J[j]==J[6-j] (j=1..5; 6-j in 1..5):")
        for j in range(1,6):
            jj=6-j
            if 1<=jj<=6:
                print(f"   J[{j}] vs J[{jj}]: {'EQUAL' if J[j]==J[jj] else 'DIFFER'}")
        # also test the full alternating inclusion-exclusion p0 = sum (-1)^|A| J(A)
        # over the 6 inner sectors, exact, as a sanity reconstruction:
        from itertools import combinations
        p0=Fraction(0)
        for r in range(7):
            for A in combinations(range(1,7),r):
                p0 += (-1)**r * avoid_pair_meas(C,A)
        print(f"reconstructed p0 (incl-excl) = {p0} = {float(p0):.6f}")
