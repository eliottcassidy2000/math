#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ROUTE 2 -- EXACT T+(a), T-(a) formula for consec at general k.
Compute T+/T- for consec at k=8..16 and look for the closed pattern in (a,k)."""
import sys
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
HALF=F(1,14)
def cov_right(E,a):
    cov=defaultdict(list)
    for e in E:
        if e==0: cov[0].append((F(0),F(10))); continue
        j=0
        while F(j,7*e)<HALF:
            cov[(e*a+j)%7].append((F(j,7*e),F(j+1,7*e))); j+=1
    return cov
def cov_left(E,a):
    cov=defaultdict(list)
    for e in E:
        if e==0: cov[0].append((F(0),F(10))); continue
        j=0
        while F(j,7*e)<HALF:
            cov[(e*a-1-j)%7].append((F(j,7*e),F(j+1,7*e))); j+=1
    return cov
def fg(iv):
    iv=sorted(iv); cur=F(0)
    for lo,hi in iv:
        if lo>cur: return cur
        cur=max(cur,hi)
    return min(cur,HALF)
def Tp(E,a): return min(fg(cov_right(E,a)[s]) for s in range(7))
def Tm(E,a): return min(fg(cov_left(E,a)[s]) for s in range(7))
def consec(k): return list(range(k))

print("consec T+(a), T-(a) for a=1..6 across k:")
for k in range(8,17):
    C=consec(k)
    tps=[Tp(C,a) for a in range(1,7)]
    tms=[Tm(C,a) for a in range(1,7)]
    win=sum(tps)+sum(tms)
    print(f" k={k:2d}: T+ = {[str(t) for t in tps]}")
    print(f"        T- = {[str(t) for t in tms]}  WIN={win}={float(win):.6f}")

# look at b+ = 1/(7 T+) and b- = 1/(7 T-)
print("\nbinding speeds b+=1/(7T+):")
for k in range(8,17):
    C=consec(k)
    bp=[(1/(7*Tp(C,a))) for a in range(1,7)]
    bm=[(1/(7*Tm(C,a))) for a in range(1,7)]
    print(f" k={k:2d}: b+ = {[str(b) for b in bp]}")
    print(f"        b- = {[str(b) for b in bm]}")
