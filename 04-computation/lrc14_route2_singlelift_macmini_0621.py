#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ROUTE 2 -- SINGLE LIFT from consec. Replace one element e by e+7m (m>=1),
keeping the multiset full-residue & distinct. Measure delta WIN. The cleanest
exchange step: e -> e+7 (smallest lift). Does it ALWAYS strictly decrease WIN?
And which residue/element lift hurts most?"""
import sys, itertools
from fractions import Fraction as F
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
def WIN(E): return sum(sum(window_TpTm(E,a)) for a in range(1,7))
def is_full_residue(E): return frozenset(e%7 for e in E)==frozenset(range(7))
def consec(k): return list(range(k))

for k in (8,9,10):
    C=consec(k); Cw=WIN(C)
    print(f"\nk={k} consec={C} WIN={float(Cw):.6f}={Cw}")
    print("  single lift e->e+7 (must stay distinct & full-residue):")
    for idx,e in enumerate(C):
        E2=sorted(C[:idx]+[e+7]+C[idx+1:])
        if len(set(E2))!=k: 
            print(f"   lift {e}->{e+7}: collides, skip"); continue
        if not is_full_residue(E2):
            print(f"   lift {e}->{e+7}: breaks full-residue, skip"); continue
        w2=WIN(E2); d=w2-Cw
        print(f"   lift {e}->{e+7}: E={E2} WIN={float(w2):.6f} dWIN={float(d):+.6f} {'DOWN' if d<0 else 'UP!!'}")
