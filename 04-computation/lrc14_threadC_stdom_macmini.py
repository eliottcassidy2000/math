#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- THE FINDING: stochastic dominance of the empty-count N.

The single-crossing run proved (0 violations, k=8,9,10, span up to 17) that
   F_consec(t) := P(N_consec <= t)  >=  P(N_E <= t) =: F_E(t)   for ALL t=0..7.
Equivalently  N_consec  <=_st  N_E  in the usual STOCHASTIC ORDER: consec's
empty-count is stochastically the SMALLEST over all competitors E.

This is a *clean one-directional certificate that STRICTLY DOMINATES the crux*:
   measS7(E) = P(N_E=0) = F_E(0) <= F_consec(0) = measS7(consec).
The crux "consec maximizes measS7" is the t=0 slice of a uniform CDF dominance.

CAUTION / earlier discrepancy: lrc14_threadC_Nlaw_orders reported "(A) consec
stochastically SMALLEST N: fails on ALL".  That used SURVIVAL G_b=P(N>=b) and an
incorrect joint test.  st-order N_consec<=_st N_E is EXACTLY
   P(N_consec>=b) <= P(N_E>=b) for all b   <=>   G_consec(b) <= G_E(b) all b
   <=>  F_consec(t) >= F_E(t) all t.
The old run tested  GC[b]<=G[b]  which IS the right inequality -- so it should NOT
fail.  This script RE-RESOLVES the discrepancy with a careful, exact recomputation
and reports the TRUE st-order status, plus the EXACT worst-case margin
   min over E, t of (F_consec(t) - F_E(t)).
If the st-order truly holds with 0 violations, that is the THREAD-C deliverable:
a NEW, PROVED-target probabilistic statement (consec st-minimizes N) that contains
the crux as its t=0 face and is amenable to a COUPLING proof.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce

def occupancy_pi(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi

def N_law(E):
    pi=occupancy_pi(E); nl=[F(0)]*8
    for h in range(8): nl[7-h]+=pi[h]
    return nl

def cdf(nl):
    F_=[F(0)]*8; acc=F(0)
    for t in range(8):
        acc+=nl[t]; F_[t]=acc
    return F_

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def run(k,W):
    C=consec(k); FC=cdf(N_law(C))
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E) and list(E)!=C]
    violations=0; worst=None; worstE=None; worstT=None
    # tight at t where margin minimal but still >=0
    minmargin=None; mm_E=None; mm_t=None
    for E in bank:
        FE=cdf(N_law(list(E)))
        for t in range(8):
            d=FC[t]-FE[t]
            if d<0:
                violations+=1
                if worst is None or d<worst: worst=d; worstE=list(E); worstT=t
            else:
                if (minmargin is None or d<minmargin) and not (t==7):  # t=7 trivially both=1
                    minmargin=d; mm_E=list(E); mm_t=t
    return dict(k=k,W=W,n=len(bank),viol=violations,worst=worst,worstE=worstE,worstT=worstT,
                minmargin=minmargin,mm_E=mm_E,mm_t=mm_t,FC=FC)

if __name__=="__main__":
    print("="*78)
    print("THREAD C FINDING: consec stochastically MINIMIZES the empty-count N")
    print("   N_consec <=_st N_E   <=>   F_consec(t) >= F_E(t)  for all t, all E")
    print("   t=0 slice = the LRC crux  measS7(consec)>=measS7(E).")
    print("="*78)
    allok=True
    for k,W in [(8,13),(9,13),(10,12),(8,15),(9,15),(10,15),(8,17),(9,17),(10,17)]:
        r=run(k,W)
        st="ST-DOMINANCE HOLDS" if r['viol']==0 else f"FAILS ({r['viol']} viol)"
        print(f"\n--- k={k}, span<={W}, competitors={r['n']} ---  {st}")
        print(f"  F_consec(t), t=0..7 = {[str(x) for x in r['FC']]}")
        if r['viol']==0:
            print(f"  tightest margin (min F_consec(t)-F_E(t), t<7) = {r['minmargin']}"
                  f" = {float(r['minmargin']):.3e}  at t={r['mm_t']}, E={r['mm_E']}")
        else:
            allok=False
            print(f"  WORST violation = {r['worst']} at t={r['worstT']}, E={r['worstE']}")
    print("\n"+"="*78)
    print(f"OVERALL: {'consec st-MINIMIZES N on ALL tested boxes (0 violations)' if allok else 'st-dominance FAILS somewhere'}")
    print("="*78)
    print("DONE.")
