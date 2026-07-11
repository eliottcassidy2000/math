#!/usr/bin/env python3
"""Slice 2: (A) sharpened far-element constant at the Phi level (one IBP on the
miss-count indicators, not 127 atoms): C_sharp = V*||H|| <= (14 Sum e')*(3/7) = 6 Sum e'.
(B) COMPRESSION TEST: does every elementary compression (move an offset one step toward
the packed interval, keeping primitivity/distinctness) weakly INCREASE Phi = p0 + p1/3?
If yes over exhaustive small + random moderate cores, the consec-extremality lemma
reduces to local-move induction -- a finite, provable skeleton."""
from fractions import Fraction as F
import random, itertools
def sector_p(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(p for p in pts if 0<=p<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
def Phi(E):
    p=sector_p(E); return p[0]+p[1]/3
print("(B) COMPRESSION MONOTONICITY TEST: e -> e-1 (toward packing), Phi should RISE")
viol=0; tests=0
random.seed(11)
# exhaustive small: all 4-subsets of {1..10}
for E in itertools.combinations(range(1,11),4):
    E=list(E)
    for i,e in enumerate(E):
        if e-1>=1 and e-1 not in E:
            E2=sorted(E[:i]+[e-1]+E[i+1:]); tests+=1
            d=Phi(E2)-Phi(E)
            if d<0: viol+=1
            if d<0 and viol<=5: print(f"  VIOLATION: {E} -> {E2}: dPhi={float(d):+.5f}")
print(f"  exhaustive 4-subsets of [1,10]: {tests} moves, {viol} violations")
v2=0; t2=0
for _ in range(60):
    k=random.randint(5,7); E=sorted(random.sample(range(1,40),k))
    for i,e in enumerate(E):
        if e-1>=1 and e-1 not in E:
            E2=sorted(E[:i]+[e-1]+E[i+1:]); t2+=1
            if Phi(E2)-Phi(E)<0: v2+=1
print(f"  random k=5..7 in [1,40): {t2} moves, {v2} violations")
