#!/usr/bin/env python3
"""mac-mini-S99b: map the GENUINE residual = covering families in NO tile: NOT near-AP (< 10 speeds in
{1..14}, so kps THM-733/734/738 j<=3 don't apply) AND no k<=13 shadow (mac-mini/klein). Record M and
diameter -- is the residual bounded-diameter (opus bounded-W finite check) and/or low-M (binding)?"""
from fractions import Fraction as F
from math import gcd
import numpy as np, random
c114=F(1,14)
def shadow_interval(S,a,k):
    lo=F(0); hi=F(1,2)
    for c in S:
        r=(c*a)%k
        if r==0: lo=max(lo,F(1,14*c)); hi=min(hi,F(13,14*c))
        else:
            s=r if r<=k//2 else r-k; base=F(abs(s),k)
            hi=min(hi,(F(13,14)-base)/c if s>0 else (base-c114)/c)
        if lo>=hi: return None
    return (lo,hi)
def lonely_at(S,t): return all(min((c*t)%1,1-(c*t)%1)>=c114 for c in S)
def has_shadow(C):
    S=[1]+list(C)
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(S,a,k)
            if iv and lonely_at(S,F(a,k)+(iv[0]+iv[1])/2): return True
    return False
def Mnum(S,N=800_000):
    t=(np.arange(1,N)/N); mn=None
    for v in S:
        r=(v*t)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    return mn.max()
def covering(C): return all(any(v%q==0 for v in C) for q in range(2,15))
def n_in_114(S): return sum(1 for v in S if 1<=v<=14)  # incl runner 1
rng=random.Random(31415)
resid=[]; checked=0
for _ in range(80000):
    # multi-killer: few small + several outliers >=15 (target <10 in {1..14})
    small=rng.sample(range(2,15), rng.randint(4,8))
    C=set(small)
    for q in range(2,15):
        if any(v%q==0 for v in C): continue
        C.add(q*rng.randint(2,15))
    while len(C)<12: C.add(rng.randint(15,220))
    C=sorted(list(C)[:12])
    if 1 in C or len(set(C))!=12 or not covering(C): continue
    S=[1]+C
    if n_in_114(S) >= 10: continue           # near-AP => kps covers, not the residual
    checked+=1
    if not has_shadow(C):                      # not shadow => genuine residual
        M=Mnum(S); resid.append((M, max(C)/min(C), sorted(C)))
    if checked>=2500: break
print(f"multi-killer (<10 in {{1..14}}) covering families checked: {checked}")
print(f"GENUINE RESIDUAL (not near-AP, not k<=13 shadow): {len(resid)}")
if resid:
    Ms=[r[0] for r in resid]; diams=[r[1] for r in resid]
    print(f"  M: min={min(Ms):.4f} median={sorted(Ms)[len(Ms)//2]:.4f} max={max(Ms):.4f}  (1/14={float(c114):.4f}, 2/23={2/23:.4f})")
    print(f"  diameter (max/min of cluster): min={min(diams):.1f} median={sorted(diams)[len(diams)//2]:.1f} max={max(diams):.1f}")
    lowM=[r for r in resid if r[0]<0.10]
    print(f"  LOW-M residual (M<0.10, near covering-min): {len(lowM)}")
    for M,d,C in sorted(lowM)[:6]: print(f"     M={M:.4f} diam={d:.1f}: {C}")
print()
print("READ: if the residual is BOUNDED-diameter and/or opus's floor (W>W0) covers the high-diameter part,")
print("then residual = a finite bounded-W check (opus's tile). If low-M UNBOUNDED-diameter residual exists")
print("outside all tiles, that is the genuine open crux. This maps exactly what is NOT yet tiled.")
