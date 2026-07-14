#!/usr/bin/env python3
"""
lrc14_ap_stability_klein_S296.py
================================
klein-2026-07-13-S296 (owner: attempt the AP-stability theorem directly).

AP-stability = "only the AP confines its good set to the ends; covering forces G(C) into the middle."
Formulation: L({1}UC)>0 <=> G(C) reaches [1/14,13/14] <=> the bad sets U D_c FAIL to cover the middle
(AP-cluster {2..13} exactly TILES it, L=0).

MECHANISM (verified): the middle good-arcs open at the RESONANCES j/k of a "released" speed k.
CLEAN ARITHMETIC (proven): t=j/k (gcd(j,k)=1, k<=13) is good for C  <=>  no c in C is a multiple of k
  [||c*j/k|| >= 1/k > 1/14 unless k|cj; gcd(j,k)=1 => k|cj iff k|c].
So the clean middle-witness t=1/k EXISTS <=> C misses a mult of k <=> {1}UC NON-COVERING at q=k = THM-523.
CENTRAL FINDING: covering = hits every q<=14 = NO clean t=1/q witness => reaches the middle only via the
subtle resonance-SHADOW gaps (a wide bad arc j/k replaced by narrower multiples) = the LRC extremal
rigidity. So the AP-stability, attempted directly, IS the covering rigidity -- characterized, not proven.
"""
import numpy as np
from math import gcd
NG=1<<23; t=np.arange(NG)/NG
def profile(C):
    m=np.full(NG,1.0)
    for c in C:
        fr=(c*t)%1.0; d=np.minimum(fr,1.0-fr); m=np.minimum(m,d)
    return m
def midgood(C):
    C=sorted(set(C)); m=profile(C); G=m>=1.0/14.0-1e-9
    mid=(t>=1.0/14.0)&(t<=0.5+1e-9); Gm=G&mid
    idx=np.where(Gm)[0]; loc=[]
    if len(idx):
        sp=np.where(np.diff(idx)>1)[0]
        ss=np.concatenate([[idx[0]],idx[sp+1]]); ee=np.concatenate([idx[sp],[idx[-1]]])
        for s,e in zip(ss,ee): loc.append(round(t[(s+e)//2],4))
    return (G&mid).sum()/NG*2, loc[:6]
def is_good_pt(C,tt): return all(min((c*tt)%1.0,1-((c*tt)%1.0))>=1.0/14.0-1e-12 for c in C)
def misses_mult(C,k): return not any(c%k==0 for c in C)

print("|G(C)∩mid| and good-arc centers in [1/14,1/2] (AP-cluster is rigid: 0):")
for C,nm in [([2,3,4,5,6,7,8,9,10,11,12,13],'{2..13} AP-cluster'),
             ([2,3,4,5,7,8,9,10,11,12,13,14],'{2..14}\\{6} (covering)'),
             ([3,4,5,6,7,8,9,10,11,12,13,14],'{3..14} (covering)'),
             ([2,3,4,5,6,7,8,9,10,12,13,14],'{2..10,12,13,14}')]:
    mg,loc=midgood(C); print("  %-22s |G∩mid|=%.5f  arcs at %s"%(nm,mg,loc))
print("\nverify t=1/k good <=> C misses mult of k (=THM-523 non-covering witness):")
for nm,C in {'{2..12,14}(noncov,no13)':[2,3,4,5,6,7,8,9,10,11,12,14],
             '{2..14}\\{6}(covering)':[2,3,4,5,7,8,9,10,11,12,13,14],
             '{2..13}AP':[2,3,4,5,6,7,8,9,10,11,12,13]}.items():
    miss=[k for k in range(2,14) if misses_mult(C,k)]
    for k in range(2,14): assert is_good_pt(C,1.0/k)==misses_mult(C,k)
    print("  %-24s clean t=1/k witness for k in %s"%(nm,miss))
print("\n=> covering hits every q<=14 => NO clean witness => shadow-gap rigidity = the covering core. done.")
