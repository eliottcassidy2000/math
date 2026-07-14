#!/usr/bin/env python3
"""
lrc14_tworesonance_klein_S299.py
================================
klein-2026-07-13-S299 (owner: close the ratio-[6,13] residual with a two-resonance witness).

FINDING: TWO resonances are NOT enough, but the FULL low-height family k in {2..13} closes it empirically.
 - k in {2..6} resonance-shadow scan covers 88-96% of covering ratio-[6,13] clusters.
 - k in {2..13} covers 100% (sampled). The uncovered-by-{2..6} clusters have their middle good points
   near t=1/k for k in {9..14} -- the LEFT-EDGE resonances (near 1/14), not t=1/2.
KEY: loneliness at a/k depends only on residues mod k, so the witness is a BOUNDED-HEIGHT rational
(k<=13) -- supports the THM-527/663 bounded-denominator realization. But the high-k shadow gaps are
MULTI-SPEED (all 12 residues must jointly leave a gap) = the equidistribution, discretized to the
resonance grid; so NO clean two-resonance closed-form. Clean provable part stays THM-744 (6 / 6-13 parity).
"""
import numpy as np, random
from math import gcd
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def is_good(S,tt): return all(min((c*tt)%1.0,1-((c*tt)%1.0))>=1.0/14.0-1e-12 for c in S)
def reson_hit(C,ks):
    S=[1]+list(C); cmax=max(C); half=1.0/cmax
    for k in ks:
        for a in range(1,k):
            if gcd(a,k)!=1: continue
            p=a/k
            if not(1.0/14.0<p<13.0/14.0): continue
            for dd in np.linspace(-half,half,400):
                tt=p+dd
                if 1.0/14.0<=tt<=13.0/14.0 and is_good(S,tt): return True
    return False
random.seed(31)
print('%-8s %4s %10s %11s'%('ratio','n','k in{2..6}','k in{2..13}'))
for band in [(6,8),(8,10),(10,13)]:
    n=k6=k13=0
    for _ in range(40000):
        cmin=random.choice([15,20,30,45,90,150]); cmax=int(cmin*random.uniform(*band))
        if cmax-cmin<11: continue
        C=sorted(random.sample(range(cmin,cmax+1),12)); S=sorted([1]+C)
        if not iscov(S): continue
        r=max(C)/min(C)
        if r<band[0] or r>band[1]: continue
        n+=1
        if reson_hit(C,[2,3,4,5,6]): k6+=1
        if reson_hit(C,list(range(2,14))): k13+=1
        if n>=60: break
    if n: print('%-8s %4d %9.0f%% %10.0f%%'%('%d-%d'%band,n,100*k6/n,100*k13/n))
print('=> two resonances insufficient; k<=13 family closes empirically; witnesses are low-height rationals.')
print('done.')
