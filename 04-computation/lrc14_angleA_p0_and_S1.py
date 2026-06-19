#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE A part 3: the two clean sub-lemmas.

(L1) consec maximizes p_0 = meas(S7) exactly (up to dilation) -- strict?
(L2) consec MINIMIZES S_1 = E[N] = sum_{j=1..6} P(sector j empty).
     S_1 = sum_{j=1}^6 meas{x: no e in E has frac(e x) in [j/7,(j+1)/7)}.
     By symmetry of the 6 inner sectors? NO -- sector 0 has e=0 fixed.
     S_1 = sum_{A: |A|=1} J(A,E). Each J({j}) = meas{x: orbit avoids sector j}.

Key structural idea to test for a PROOF of L2 (E[N] minimized by consec):
  E[N] = 6 - E[#nonzero sectors hit]. Equivalently consec MAXIMIZES expected
  number of distinct sectors covered. AP is most equidistributed => covers most.
  Try: J({j},E) = meas{x : for all e, frac(ex) not in sector j}.
  For a SINGLE sector this is a covering-by-arcs measure; test if consec min.

Also test the THREE-DISTANCE / equidistribution intuition: for AP, the orbit
{frac(e x)} = arithmetic progression on circle => three-gap theorem => most uniform.

kind-pasteur-2026-06-19 ANGLE-A.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

def frac(v): return v - (v.numerator//v.denominator)
def sect(v): return (v.numerator*7)//v.denominator

def avoid_sector_meas(E, j):
    """meas{x in [0,1): for all e in E, frac(e x) NOT in sector j} (exact)."""
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1):
            bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        ok=all(sect(frac(e*mid))!=j for e in E)
        if ok: tot+=(hi-lo)
    return tot

def dist_p(E):
    E=sorted(set(E))
    bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1):
            bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2
        hit=set(sect(frac(e*mid)) for e in E)
        t=sum(1 for jj in range(1,7) if jj not in hit)
        p[t]+=(hi-lo)
    return p

def gen_competitors(k, maxspread):
    out=[]
    for combo in itertools.combinations(range(1,maxspread+1),k-1):
        E=[0]+list(combo)
        if reduce(gcd,E)!=1: continue
        out.append(E)
    return out

def is_dilation(E,C):
    E=sorted(E)
    if E[0]!=0: return False
    g0=reduce(gcd,E)
    return [e//g0 for e in E]==C

if __name__=="__main__":
    for k in [8,9]:
        print(f"\n{'='*70}\nk={k}\n{'='*70}")
        C=list(range(k))
        pc=dist_p(C)
        p0c=pc[0]
        S1c=sum(t*pc[t] for t in range(7))
        # per-sector avoid measures for consec
        avc=[avoid_sector_meas(C,j) for j in range(1,7)]
        print(f"consec p0={p0c}={float(p0c):.6f}  E[N]=S1={S1c}={float(S1c):.6f}")
        print(f"consec per-sector avoid J({{j}}) j=1..6: {[str(a) for a in avc]}")
        print(f"  ~ {[f'{float(a):.5f}' for a in avc]}  (sum={float(sum(avc)):.5f}=S1)")

        comps=gen_competitors(k,{8:14,9:13}[k])
        p0_strictbeat=0; p0_nondil_tie=0
        S1_strictbeat=0; S1_nondil_tie=0   # competitor S1 < consec S1 (consec should be MIN)
        for E in comps:
            if E==C: continue
            p=dist_p(E)
            if p[0]>p0c: p0_strictbeat+=1
            elif p[0]==p0c and not is_dilation(E,C): p0_nondil_tie+=1
            S1=sum(t*p[t] for t in range(7))
            if S1<S1c: S1_strictbeat+=1
            elif S1==S1c and not is_dilation(E,C): S1_nondil_tie+=1
        print(f"  p0: competitors beating consec p0 = {p0_strictbeat}; non-dilation ties = {p0_nondil_tie}")
        print(f"  S1=E[N]: competitors with SMALLER E[N] than consec = {S1_strictbeat}; non-dil ties = {S1_nondil_tie}")
        print(f"  => consec {'MAXIMIZES p0 (strict up to dilation)' if p0_strictbeat==0 and p0_nondil_tie==0 else 'does NOT uniquely max p0'}")
        print(f"  => consec {'MINIMIZES E[N] (strict up to dilation)' if S1_strictbeat==0 and S1_nondil_tie==0 else 'does NOT uniquely min E[N]'}")
