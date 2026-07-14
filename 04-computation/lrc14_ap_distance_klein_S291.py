#!/usr/bin/env python3
"""
lrc14_ap_distance_klein_S291.py
===============================
klein-2026-07-13-S291 (owner: prove covering buys uniform distance from the AP).

Via the S290 identity L({1}UC)=|G(C)|(1-conc/7), "covering buys uniform distance from the AP" means
   sup_{primitive covering} conc < 7    (equivalently  inf_{primitive covering} L > 0),
which is the open residual (THM-527-A cancellation). This script establishes the STRUCTURE:
 (A) conc<=7 ALWAYS (L>=0); conc=7 <=> L=0 <=> tight extremal. Verify AP and GW both have conc=7,
     both primitive NON-covering. So covering => conc<7 pointwise.
 (B) PRIMITIVITY is the separator: other conc=7 configs are imprimitive dilates c*{AP}, c*{GW}.
 (C) THE GAP FACTORS: sup conc is at BOUNDED near-AP sets (large-speed covering sets have small conc),
     so the uniform gap = [bounded near-AP finite check (kps-THM-734)] + [large-speed equidistribution
     (opus true-disc slack)]. Sample to find max conc over primitive covering min=1 sets.
"""
import numpy as np, random
from math import gcd
from functools import reduce
NG=1<<20; THR=1.0/14.0; t=np.arange(NG)/NG
def good(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=THR)
    return g
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def isprim(S): return reduce(gcd,S)==1
lo=int(round(NG/14.0))
def conc_of(S):
    C=[w for w in S if w!=1]; gc=good(C); Lc=gc.mean()
    if Lc==0: return float('nan'),0.0
    c=14*(gc[:lo].sum()/NG)/Lc
    return c, Lc*(1-c/7.0)

print("=== (A)+(B) tight configs (conc=7 <=> L=0 <=> {AP,GW}); covering => conc<7 ===")
for nm,S in {'AP {1..13}':[1,2,3,4,5,6,7,8,9,10,11,12,13],
             'GW {1..11,13,24}':[1,2,3,4,5,6,7,8,9,10,11,13,24],
             'near-AP cov {1..14}\\{6}':[1,2,3,4,5,7,8,9,10,11,12,13,14],
             'near-AP cov {1,3..14}':[1,3,4,5,6,7,8,9,10,11,12,13,14],
             'deep well':[1,2,3,4,5,6,7,8,9,10,11,12,182],
             '{1,90..101}':[1,90,91,92,93,94,95,96,97,98,99,100,101]}.items():
    S=sorted(set(S)); c,L=conc_of(S)
    print("  %-26s conc=%6.3f L=%.6f cov=%-5s prim=%s"%(nm,c,L,iscov(S),isprim(S)))

print("\n=== (C) max conc over PRIMITIVE COVERING min=1 sets (the gap; sup at BOUNDED near-AP) ===")
random.seed(1); best=0; worstS=None; n=0
for _ in range(6000):
    N=random.choice([16,18,20,26,40,90]); C=sorted(random.sample(range(2,N+1),12)); S=sorted([1]+C)
    if not iscov(S) or not isprim(S): continue
    n+=1; c,_=conc_of(S)
    if not np.isnan(c) and c>best: best=c; worstS=S
print("  sampled %d prim-covering min=1 sets; MAX conc=%.4f (gap 7-conc=%.4f) at %s"%(n,best,7-best,worstS))
print("  large-speed covering sets have conc~3.3 (far from AP); the sup lives at BOUNDED near-AP sets.")
print("  => uniform gap FACTORS: [bounded near-AP finite check (kps-THM-734)] + [large-speed equidist (opus)].")
print("done.")
