#!/usr/bin/env python3
"""mac-mini-S101d: census the THM-751 PEELING closure. Peel the largest speed w from a covering family;
at the core's tight point t0, classify w: ALIGNED (w*t0 in Z, tooth-narrowing => M rises to M(core)),
NON-ALIGNED-SAFE (||w*t0||>=M(core) => M(family)=M(core), recurse on core), or NON-ALIGNED-UNSAFE
(||w*t0||<M(core) => tight point shifts, the residual gap). How often does each occur?"""
from fractions import Fraction as F
from math import gcd
import numpy as np, random
def M_t0(S,dens=1_500_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def classify_peel(S):
    w=max(S); core=[v for v in S if v!=w]
    mu0,t0=M_t0(core)
    # ||w t0||
    r=(w*t0)%1.0; wn=min(r,1-r)
    if wn < 1.5/max(w,1):  # ~ aligned (tooth: w*t0 near integer)
        return "aligned"
    if wn >= mu0 - 2e-4:
        return "nonaligned-safe"
    return "nonaligned-UNSAFE"
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
rng=random.Random(99)
cnt={"aligned":0,"nonaligned-safe":0,"nonaligned-UNSAFE":0}; unsafe_ex=[]; tot=0
# structured extremals + near-AP-with-far + random covering
fams=[[v for v in range(1,14) if v!=6]+[182], [v for v in range(1,14) if v!=6]+[169],
      list(range(1,13))+[182], list(range(1,12))+[13,84], list(range(1,11))+[13,14,45]]
for _ in range(2500):
    S=set()
    for q in range(2,15):
        if any(v%q==0 for v in S): continue
        S.add(q*rng.randint(1,18))
    while len(S)<13: S.add(rng.randint(1,300))
    S=sorted(set(list(S)[:13])); 
    if len(S)==13 and 1 in S: fams.append(S)
for S in fams:
    if 1 not in S or len(set(S))!=13 or not covering(S): continue
    tot+=1; c=classify_peel(S); cnt[c]+=1
    if c=="nonaligned-UNSAFE" and len(unsafe_ex)<8: unsafe_ex.append(sorted(S))
print(f"covering 13-sets peeled (largest speed): {tot}")
for k,v in cnt.items(): print(f"  {k}: {v} ({100*v/tot:.1f}%)")
print(f"\nNON-ALIGNED-UNSAFE (the residual gap after one peel): {cnt['nonaligned-UNSAFE']}")
for S in unsafe_ex[:6]:
    Mm,_=M_t0(S); print(f"   M={Mm:.4f}: {S}")
print()
print("THM-751 two branches (aligned tooth + nonaligned-safe->core) close the aligned & safe peels;")
print("recursion on the core terminates at a NON-COVERING core (M>=1/14 by t=1/q sieve, THM-366) or")
print("bounded outliers (finite check). The residual = repeated NON-ALIGNED-UNSAFE peels = opus's")
print("spread stratum. Frequency of unsafe measures how much is left for the density floor.")
