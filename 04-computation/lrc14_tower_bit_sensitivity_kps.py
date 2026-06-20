#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which tower bit is 'cheapest' to break?  For each single tower bit b in {1,2,4,8}, find the
minimum lonely measure over primitive 12-cores in [1,cap] that contain the OTHER three tower
bits but NOT b.  Reveals the per-bit 'cost' of leaving the mouth.  Writes to file. EXACT.
kind-pasteur-2026-06-19.
"""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd
THETA=Fraction(1,14); THR2=Fraction(426,35035); THR1=Fraction(7,858); TOWER=(1,2,4,8)
def meas(C):
    segs=[]
    for d in C:
        w=THETA/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for s in (-1,0,1):
                a=max(lo+s,Fraction(0)); b=min(hi+s,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); cur=Fraction(-1); U=[]
    for a,b in segs:
        if a>cur: U.append([a,b]); cur=b
        elif b>cur: U[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in U)
out=open("05-knowledge/results/lrc14_tower_bit_sensitivity_kps.out","w",encoding="utf-8")
def P(*a): print(*a,file=out); out.flush()
cap=20
P(f"THR1(mouth)={float(THR1):.6f}  THR2(2nd)={float(THR2):.6f}  cap={cap}")
for kill in TOWER:
    keep=[b for b in TOWER if b!=kill]
    best=None; n=0
    for C in itertools.combinations(range(1,cap+1),12):
        if kill in C: continue
        if not all(k in C for k in keep): continue
        if reduce(gcd,C)!=1: continue
        n+=1
        L=meas(C)
        if best is None or L<best[0]: best=(L,C)
    P(f"  break ONLY bit {kill} (keep {keep}): cores={n} min_meas={best[0]}={float(best[0]):.6f} C={best[1]} >=THR2? {best[0]>=THR2}")
out.close()
