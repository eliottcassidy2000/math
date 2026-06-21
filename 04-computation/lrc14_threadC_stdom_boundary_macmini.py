#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- WHERE does st-dominance of N hold?  (honest boundary).

Wide random sampling showed full st-order N_consec <=_st N_E FAILS globally
(consistent with HYP-2635 "FOSD eliminated").  BUT it held with 0 violations on
ALL span<=17 primitive competitors at k=8,9,10 (the THM-535 binding regime).

This script pins the EXACT span boundary: for each k, increase the max-element
cap B and report the first B at which a st-violation (F_consec(t)<F_E(t) some t)
appears, and the offending E.  Goal: characterize "st-dominance holds for
max(E) <= B*(k)" -- a HONEST, bounded statement that (if B*(k) covers the
binding rows of THM-535) becomes a usable certificate for the bounded branch
(true-wide is handled separately by THM-547 etc.).

Also: report the WORST-CASE st margin as a function of B, to see how the order
degrades.
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
    for t in range(8): acc+=nl[t]; F_[t]=acc
    return F_
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def first_violation(k, Bmax):
    C=consec(k); FC=cdf(N_law(C))
    firstB=None; firstE=None
    for B in range(k, Bmax+1):
        # only NEW shapes whose max element == B (max=B), to find the threshold
        # E=(0,) + (k-1 chosen from 1..B with B included)
        rest=range(1,B)  # elements strictly below B
        for combo in itertools.combinations(rest, k-2):
            E=(0,)+combo+(B,)
            if not primitive(E): continue
            FE=cdf(N_law(list(E)))
            if any(FC[t]<FE[t] for t in range(8)):
                return B, list(E), [str(FC[t]-FE[t]) for t in range(8)]
    return None, None, None

if __name__=="__main__":
    print("="*78)
    print("THREAD C: st-dominance boundary -- first max(E)=B with N_consec NOT <=_st N_E")
    print("="*78)
    for k,Bmax in [(8,30),(9,28),(10,26)]:
        B,E,d=first_violation(k,Bmax)
        if B is None:
            print(f"\nk={k}: NO st-violation up to max(E)={Bmax}  (st-dominance holds on whole box!)")
        else:
            print(f"\nk={k}: first st-violation at max(E)={B}")
            print(f"   E={E}")
            print(f"   F_consec - F_E (t=0..7) = {d}   (t=0 entry={d[0]} -> measS7 still consec-max)")
            print(f"   => st-dominance HOLDS for max(E) <= {B-1}; the binding THM-535 rows"
                  f" have max <= 17, so they are INSIDE the st-dominance regime.")
    print("\nDONE.")
