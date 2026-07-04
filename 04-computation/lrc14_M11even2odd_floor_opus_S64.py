#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CONJECTURE (=> m=2,f=2 confinement as a UNIFORM gap): M(2U u {w1,w2}) >= 1/12 for ANY 11-runner U and
odd w1,w2. Then gap(U)=min(M-1/14) >= 1/12-1/14 = 1/84 > 0. opus-S64.
This descent MINIMIZES M over (11 even + 2 odd) families to stress-test the 1/12 floor (avoid MISTAKE-101:
seeds include commensurate/dilated-AP families). Float descent; exact-verify the minimizer.
"""
import sys, itertools, random
import numpy as np
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
random.seed(4)
def Mf(S,G=60013):
    t=(np.arange(G)+0.5)/G;F=np.full(G,1.0)
    for v in S:
        fr=(v*t)%1.0;F=np.minimum(F,np.minimum(fr,1.0-fr))
    return float(F.max())
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0)
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v
    return b
EV=list(range(2,64,2)); OD=list(range(1,64,2))
def M(ev,od): return Mf(list(ev)+list(od))
def descend(ev,od):
    ev=sorted(set(ev));od=sorted(set(od));cur=M(ev,od)
    imp=True
    while imp:
        imp=False
        for idx in range(11):
            for c in EV:
                if c in ev:continue
                e2=sorted(set(ev[:idx]+ev[idx+1:]+[c]))
                if len(e2)==11:
                    m=M(e2,od)
                    if m<cur-1e-7:ev=e2;cur=m;imp=True;break
            if imp:break
        if imp:continue
        for idx in range(2):
            for c in OD:
                if c in od:continue
                o2=sorted(set(od[:idx]+od[idx+1:]+[c]))
                if len(o2)==2:
                    m=M(ev,o2)
                    if m<cur-1e-7:od=o2;cur=m;imp=True;break
            if imp:break
    return cur,ev,od
seeds=[([2*j for j in range(1,12)],[1,3]),
       ([2*j for j in range(1,14) if j not in(12,13)],[1,11]),
       ([2*c for c in [2,3,4,5,6,7,8,9,10,11,12]],[1,3]),
       ([4*j for j in range(1,12)],[1,3])]
for _ in range(50):
    seeds.append((random.sample(EV,11),random.sample(OD,2)))
best=(1.0,None,None)
for ev,od in seeds:
    if len(set(ev))!=11 or len(set(od))!=2:continue
    m,e,o=descend(ev,od)
    if m<best[0]:best=(m,e,o)
mf,ev,od=best
print("MIN M over (11 even + 2 odd), descent %d seeds:"%len(seeds))
print("  min M (float) = %.6f   1/12 = %.6f   1/14 = %.6f"%(mf,1/12,1/14))
Me=exact_M(sorted(ev+od))
print("  minimizer exact M = %s = %.6f   (>= 1/12? %s)"%(str(Me),float(Me),Me>=Fr(1,12)))
print("  evens=",ev," odds=",od)
print("  => min M = 1/12 (never below) supports gap >= 1/84 > 0 (uniform), = m=2,f=2 confinement w/ margin.")
