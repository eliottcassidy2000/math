#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Does the f=2 tightening gap GROW with u_max? (=> v_max(U) effectively bounded => finite confinement check).
opus-S62. min over random 11-runner even parts U' (by u_max band) of gap(U')=min_{odd w1,w2}(M(2U' u {w1,w2})-1/14).
gap=0 <=> tight q*=28 exists. If MIN gap stays bounded away from 0 and grows with u_max, large even parts are
robustly infeasible -> u_max is effectively bounded (supports mac-mini's finite per-U check).
"""
import sys, itertools, random
import numpy as np
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
random.seed(2)
def M_float(S, G=8009):
    t=(np.arange(G)+0.5)/G; F=np.full(G,1.0)
    for v in S:
        fr=(v*t)%1.0; F=np.minimum(F,np.minimum(fr,1.0-fr))
    return float(F.max())
def prim(S): return reduce(gcd,S)==1
TARGET=1.0/14
def gap(Up, wmax):
    E=[2*u for u in Up]; odds=[w for w in range(1,wmax+1,2) if w not in E]
    best=1.0
    for w1,w2 in itertools.combinations(odds,2):
        S=E+[w1,w2]
        if not prim(S): continue
        M=M_float(S)
        if M>=TARGET-1e-4:
            g=M-TARGET
            if g<best: best=g
    return best
print("MIN f=2 tightening gap over random 11-runner even parts U', by u_max band (G=8009 float)")
print("  gap=0 => tight q*=28 family exists; gap>0 with margin => robustly infeasible at that u_max")
print("  %12s %9s %10s %11s" % ('u_max band','#sampled','MIN gap','median gap'))
for lo,hi in [(11,13),(14,18),(19,25),(26,34),(35,48)]:
    gaps=[]; tries=0
    while len(gaps)<18 and tries<3000:
        tries+=1
        umax=random.randint(lo,hi)
        Up=sorted(set(random.sample(range(1,umax),10)+[umax]))
        if len(Up)!=11: continue
        gaps.append(gap(Up, min(41,6*umax)))
    gaps.sort()
    print("  %12s %9d %10.5f %11.5f" % (f"{lo}-{hi}", len(gaps), min(gaps), gaps[len(gaps)//2]))
print("READING: MIN gap rising (or bounded away from 0) with u_max => no tight q*=28 at large u_max =>")
print("  u_max effectively bounded; the near-feasible threats sit at SMALL u_max (AP-like, near the")
print("  imprimitive even block). So bounding u_max ~ ruling out small AP-like even parts = the rigidity.")
print("DONE.")
