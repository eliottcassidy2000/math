#!/usr/bin/env python3
"""Is max L_yK8 over size-k configs MONOTONE DECREASING in far-count r=#{e>14}?
If max_{r far} L_yK8 is decreasing in r, then r=0 (bounded) is the global max => gK8 concentration,
and each step r->r+1 is a single-far-type comparison (periodicity domain).
Adversarial search of max L_yK8 per far-count."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(7)
def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if t<=6: q[t]+=x1-x0
    return q
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
for k in (8,9,10):
    print(f"--- k={k} (cap*10={float(10*caps[k]):.4f}) ---")
    best={}
    # r=0: bounded (search bounded configs)
    bestL={0:F(0),1:F(0),2:F(0),3:F(0)}; bestE={0:None,1:None,2:None,3:None}
    # bounded r=0
    for _ in range(400):
        E=tuple(sorted([0]+random.sample(range(1,15),k-1)))
        L=Lyk8(missdist(E))
        if L>bestL[0]: bestL[0]=L; bestE[0]=E
    # include consec explicitly
    L=Lyk8(missdist(tuple(range(k))))
    if L>bestL[0]: bestL[0]=L; bestE[0]=tuple(range(k))
    # r=1,2,3 far: base bounded + r far elements
    for r in (1,2,3):
        for _ in range(500):
            nb=k-r
            if nb<1: continue
            base=sorted([0]+random.sample(range(1,15),nb-1))
            fars=sorted(random.sample(range(15,80),r))
            E=tuple(sorted(base+fars))
            L=Lyk8(missdist(E))
            if L>bestL[r]: bestL[r]=L; bestE[r]=E
    for r in (0,1,2,3):
        print(f"  r={r} far: max L_yK8={float(bestL[r]):.4f}  @ {str(bestE[r])[:40]}")
    mono = bestL[0]>=bestL[1]>=bestL[2]>=bestL[3]
    print(f"  monotone decreasing in r? {mono}  (=> bounded r=0 is global max => concentration)")
