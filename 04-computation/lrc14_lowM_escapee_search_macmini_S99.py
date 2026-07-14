#!/usr/bin/env python3
"""mac-mini-S99: the MAKE-OR-BREAK test for the S98 dichotomy. Search adversarially for a LOW-M escapee:
a covering 13-set with NO k<=13 shadow witness AND M near the covering-min (< 0.15). If none exists,
'escapee => loose' holds (shadow closes all binding families); if one exists, the shadow route misses a
binding case (genuine gap). Also confirm all LOW-M covering families HAVE a k<=13 shadow."""
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
    return all(min((c*t)%1,1-(c*t)%1)>=c114 for c in S)
def has_shadow(C):
    S=[1]+list(C)
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(S,a,k)
            if iv and lonely_at(S,F(a,k)+(iv[0]+iv[1])/2): return True
    return False
def Mnum(S,N=1_000_000):
    t=(np.arange(1,N)/N); mn=None
    for v in S:
        r=(v*t)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    return mn.max()
def covering(C): return all(any(v%q==0 for v in C) for q in range(2,15))
rng=random.Random(2718)
# adversarial: bias toward LOW M -- perturbations of covering-min extremals + near-AP + tight-ish spread
seeds=[list(range(2,13))+[182], [*range(2,12),13,84], [v for v in range(2,15) if v!=6],
       [*range(2,13),182],[2,3,4,5,6,7,8,9,10,13,22,84]]
min_escapee_M=(9,None); low_shadowless=[]; checked=0; nesc=0
def process(C):
    global min_escapee_M,nesc,checked
    if 1 in C or len(set(C))!=len(C) or not covering(C): return
    checked+=1
    if not has_shadow(C):
        nesc+=1; M=Mnum([1]+C)
        if M<min_escapee_M[0]: min_escapee_M=(M,sorted(C))
        if M<0.15: low_shadowless.append((round(M,4),sorted(C)))
# (1) perturb the extremals (small edits -> stay low-M)
for _ in range(4000):
    base=list(rng.choice(seeds)); C=set(base)
    for _ in range(rng.randint(1,3)):
        if C: C.discard(rng.choice(list(C)))
        C.add(rng.choice([rng.randint(2,60), rng.choice([13,26,14,28,84,182,169])]))
    while len(C)<12: C.add(rng.randint(2,60))
    process(sorted(list(C)[:12]))
# (2) near-AP + one/two far
for _ in range(3000):
    drop=rng.sample(range(2,15), rng.randint(1,3)); C=[v for v in range(2,15) if v not in drop]
    for _ in range(len(drop)): C.append(rng.choice([13,26,84,182,14,28,169,90,91]))
    process(sorted(set(C))[:12] if len(set(C))>=12 else None) if len(set(C))>=12 else None
# (3) broad random (context)
for _ in range(4000):
    C=set()
    for q in range(2,15):
        if any(v%q==0 for v in C): continue
        C.add(q*rng.randint(1,20))
    while len(C)<12: C.add(rng.randint(2,250))
    process(sorted(list(C)[:12]))
print(f"covering families checked: {checked}; k<=13-shadow ESCAPEES: {nesc}")
print(f"MINIMUM M among escapees: {min_escapee_M[0]:.4f}  at {min_escapee_M[1]}")
print(f"  (1/14={float(c114):.4f}; covering-min 14/183={14/183:.4f}; 1/13={1/13:.4f})")
print(f"escapees with M<0.15 (would BREAK the dichotomy): {len(low_shadowless)}")
for M,C in sorted(low_shadowless)[:8]: print(f"    M={M}: {C}")
if not low_shadowless:
    print("=> NO low-M escapee found: 'escapee => M large' holds under adversarial low-M search.")
    print("   Dichotomy ROBUST: shadow closes the BINDING low-M region; escapees are loose (opus mid-band).")
