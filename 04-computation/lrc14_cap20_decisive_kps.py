#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decisive cap=20 test: tower-broken min + all sub-THR2 cores. Writes result to a file."""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd
THETA=Fraction(1,14); THR2=Fraction(426,35035); TOWER=(1,2,4,8)
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
out=open("05-knowledge/results/lrc14_cap20_decisive_kps.out","w",encoding="utf-8")
def P(*a): print(*a,file=out); out.flush()
for cap in [20,22]:
    best=None; n=0
    for C in itertools.combinations(range(1,cap+1),12):
        if all(b in C for b in TOWER): continue
        if reduce(gcd,C)!=1: continue
        n+=1
        L=meas(C)
        if best is None or L<best[0]: best=(L,C)
    P(f"cap={cap} tower-BROKEN primitive cores: {n}")
    P(f"  min meas={best[0]}={float(best[0]):.6f} C={best[1]}")
    P(f"  THR2={float(THR2):.6f} min>=THR2? {best[0]>=THR2}")
    below=[]
    for C in itertools.combinations(range(1,cap+1),12):
        if reduce(gcd,C)!=1: continue
        L=meas(C)
        if L<THR2: below.append((float(L),C,all(b in C for b in TOWER)))
    below.sort()
    P(f"  sub-THR2 cores in [1,{cap}]: {len(below)}")
    for L,C,tw in below: P(f"     meas={L:.6f} C={C} tower_full={tw}")
    P("")
out.close()
