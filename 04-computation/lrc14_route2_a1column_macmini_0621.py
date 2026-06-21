#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ROUTE 2 -- the a=1 column. At a=1, sector s native occupant is residue s (since a=1,
e*a==e). Going right, sector s covered by clock e=s native, leaves at 1/(7*minmag(s)),
refilled by clock of residue s-1 entering at 1/(7*e_{s-1}) shifted... At a=1 the geometry
is the cleanest. Test: is consec the UNIQUE max of (T+(1)+T-(1)) over full-residue? And
is the a=1 binding speed EXACTLY max_{r in 1..6} minmag(r) (=6 for consec)?"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)
def sector_of(e,a,y):
    pos=F(e*a)+F(7*e)*y; return (pos.numerator//pos.denominator)%7
def covered(E,a,y): return len({sector_of(e,a,y) for e in E})==7
def breakpoints(E,a):
    bps={F(0),-HALF,HALF}
    for e in E:
        if e==0: continue
        lo_val=F(7*e)*(-HALF)+F(e*a); hi_val=F(7*e)*(HALF)+F(e*a)
        lo_i=min(lo_val,hi_val);hi_i=max(lo_val,hi_val); m=lo_i.numerator//lo_i.denominator
        while m<=hi_i.numerator//hi_i.denominator+1:
            y=F(m-e*a,7*e)
            if -HALF<=y<=HALF: bps.add(y)
            m+=1
    return sorted(bps)
def window_TpTm(E,a):
    bps=breakpoints(E,a); ivals=list(zip(bps,bps[1:])); Tp=F(0)
    for lo,hi in ivals:
        if hi<=0: continue
        lo2=max(lo,F(0))
        if covered(E,a,(lo2+hi)/2): Tp=hi
        else:
            if lo2==F(0): Tp=F(0)
            break
    Tp=min(Tp,HALF); Tm=F(0)
    for lo,hi in reversed(ivals):
        if lo>=0: continue
        hi2=min(hi,F(0))
        if covered(E,a,(lo+hi2)/2): Tm=-lo
        else:
            if hi2==F(0): Tm=F(0)
            break
    Tm=min(Tm,HALF); return Tp,Tm
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))
def minmag(E):
    d={}
    for e in E:
        r=e%7
        if r not in d or abs(e)<d[r]: d[r]=abs(e)
    return d

# Per-resonance a, is consec the unique max of T+(a)+T-(a)?
for k,span in [(8,14),(9,16),(10,18)]:
    C=consec(k)
    bank=[list(c) for c in itertools.combinations(range(0,span+1),k) if 0 in c and is_full_residue(c)]
    print(f"\nk={k}: per-resonance consec-max of (T+(a)+T-(a))?")
    for a in range(1,7):
        Cval=sum(window_TpTm(C,a))
        beat=[E for E in bank if sum(window_TpTm(E,a))>Cval+F(1,10**12)]
        tie=[E for E in bank if sum(window_TpTm(E,a))==Cval and E!=C]
        print(f"  a={a}: consec (T++T-)={float(Cval):.6f}; beaters={len(beat)} ties={len(tie)}")
