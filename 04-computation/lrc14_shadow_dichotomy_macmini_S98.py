#!/usr/bin/env python3
"""mac-mini-S98b: HONEST correction + the dichotomy test. My S97 'k<=13 shadow closes ALL covering' is
FALSE (high-ratio spread clusters escape). Test the SALVAGE: is every k<=13-shadow-ESCAPEE loose (large
M)? If min-M over escapees >> 1/14, then covering => (k<=13 shadow, low-M) OR (M large, loose) -- both
give M>=1/14, a clean dichotomy. Find the LOWEST-M escapee (the danger: a low-M escapee = genuine gap)."""
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
def lonely_at(S,t):
    for c in S:
        r=(c*t)%1
        if min(r,1-r)<c114: return False
    return True
def has_shadow(C):
    S=[1]+list(C)
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(S,a,k)
            if iv and lonely_at(S,F(a,k)+(iv[0]+iv[1])/2): return True
    return False
def Mnum(S,N=1_500_000):
    t=(np.arange(1,N)/N); mn=None
    for v in S:
        r=(v*t)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    return mn.max()
def covering(C): return all(any(v%q==0 for v in C) for q in range(2,15))
rng=random.Random(101)
escapee_Ms=[]; low_escapees=[]; nshad=0; ntot=0
for _ in range(50000):
    C=set()
    for q in range(2,15):
        if any(v%q==0 for v in C): continue
        C.add(q*rng.randint(1,30))
    while len(C)<12: C.add(rng.randint(2,400))
    C=sorted(list(C)[:12])
    if 1 in C or not covering(C): continue
    ntot+=1
    if not has_shadow(C):
        nshad+=1; M=Mnum([1]+C); escapee_Ms.append(M)
        if M < 0.12: low_escapees.append((round(M,4),C))
    if ntot>=3000: break
print(f"covering clusters: {ntot}; k<=13-shadow ESCAPEES: {nshad} ({100*nshad/ntot:.0f}%)")
if escapee_Ms:
    print(f"escapee M: min={min(escapee_Ms):.4f}, median={sorted(escapee_Ms)[len(escapee_Ms)//2]:.4f}, max={max(escapee_Ms):.4f}")
    print(f"  (1/14 = {1/14:.4f}; covering-min 14/183 = {14/183:.4f})")
    print(f"escapees with M < 0.12 (would threaten the dichotomy): {len(low_escapees)}")
    for M,C in sorted(low_escapees)[:6]: print(f"    M={M}: {C}")
print()
print("DICHOTOMY: if min escapee-M >> 1/14, then covering => (k<=13 shadow witness, the LOW-M/binding")
print("families incl single-killer S97, packed ratio-[6,13] klein-S299) OR (M large, loose) -- a clean")
print("shadow-or-loose split. A LOW-M escapee would be a genuine gap (shadow misses a binding case).")
