#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 2 -- the binding-speed MULTISET and the exact functional WIN = sum 1/(7 b).
We have 12 values T_i (a=1..6, side +/-). Define binding speed b_i = 1/(7 T_i) for
T_i not capped/empty. WIN = sum_i T_i = (1/7) sum_i 1/b_i = harmonic sum.

Goal: characterize consec's binding-speed multiset and test the exchange direction.
Print consec's 12 b's at k=8,9,10 and the failing-shape's b's side by side.
"""
import sys, itertools
from fractions import Fraction as F
from collections import Counter
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
def bspeeds(E):
    out=[]
    for a in range(1,7):
        Tp,Tm=window_TpTm(E,a)
        for T in (Tp,Tm):
            if T==0: out.append(('EMPTY',F(0)))
            elif T==HALF: out.append(('CAP',F(7,2)))  # b=1/(7*1/14)=2
            else: out.append((str(1/(7*T)),T))
    return out
def consec(k): return list(range(k))
def WIN(E): return sum(sum(window_TpTm(E,a)) for a in range(1,7))

for k in (8,9,10):
    C=consec(k)
    bs=bspeeds(C)
    bvals=Counter(b for b,T in bs)
    print(f"\nk={k} consec WIN={WIN(C)}={float(WIN(C)):.6f}")
    print(f"  binding-speed multiset (b=1/(7T)): {dict(bvals)}")
    print(f"  Tvalues sorted: {sorted((float(T) for b,T in bs), reverse=True)}")

print("\n--- failing shape k=10 [0..8,12] ---")
E=[0,1,2,3,4,5,6,7,8,12]
bs=bspeeds(E); print(f"  WIN={WIN(E)}={float(WIN(E)):.6f}")
print(f"  binding-speed multiset: {dict(Counter(b for b,T in bs))}")
