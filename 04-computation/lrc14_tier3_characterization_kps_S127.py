# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont39: CHARACTERIZING the localized near-unit tier-3 of the bounded-window covering.
# Tier-3 = non-detuned, prime-rich families with NO in-window ruler at q<=14, needing a near-unit q in [15,43].
# FINDINGS:
#  * REDUCTION: extending tier-1 from the CLEAN lemma (maxBand<=5) to the EXACT B5=liveCount-penalty (THM-707)
#    at q<=14 absorbs the #even=6 families (bandCount=6 with liveCount>penalty), shrinking tier-3 ~6% -> ~1.7%.
#  * The tier-3 residue is LOOSE (M~0.19, tight=1/14=0.071) -- genuinely lonely, just hard to DETECT with q<=14.
#  * STRUCTURE: some divisor hits EXACTLY 6 of the 13 speeds -- the BOUNDARY between tier-1 (<=5, clean) and
#    tier-2 (>=7, detuned). Not near-AP (longest-AP ~3.7). Detected by ~3 near-unit moduli in [15,43].
from math import comb, gcd
from fractions import Fraction as F
import random
from collections import Counter
def B5(v,q):
    S=[0]*6
    for p in range(1,q):
        b=sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
        for d in range(6): S[d]+=comb(b,d)
    return sum((-1)**d*S[d] for d in range(6))
def tier1_exact(v): return any(B5(v,q)>0 for q in range(8,15))
def heavy(v): return any(sum(1 for vi in v if vi%g==0)>=7 for g in [2,3,5,7])
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def M(v):
    best=F(-1)
    for i in range(13):
        for j in range(i+1,13):
            q=v[i]+v[j]
            for p in range(1,q):
                if gcd(p,q)==1:
                    m=min(norm(F(vi*p,q)) for vi in v)
                    if m>best: best=m
    return float(best)

def main():
    random.seed(11); n=8000; te=t2=t3=0
    for _ in range(n):
        v=[random.randrange(1,80) for _ in range(13)]
        if tier1_exact(v): te+=1
        elif heavy(v): t2+=1
        elif any(B5(v,q)>0 for q in range(15,44)): t3+=1
    print(f"COVERING with EXACT B5 tier-1 ({n} families):")
    print(f"  tier1 EXACT B5 at q<=14:  {100*te/n:.1f}%   tier2 detuned: {100*t2/n:.1f}%   tier3 near-unit: {100*t3/n:.2f}%")
    random.seed(4); t3f=[]
    while len(t3f)<40:
        v=sorted(random.sample(range(1,70),13))
        if not tier1_exact(v) and not heavy(v) and any(B5(v,q)>0 for q in range(15,44)): t3f.append(v)
    Ms=[M(v) for v in t3f]
    md=[max(sum(1 for x in v if x%g==0) for g in [2,3,5,7,11,13]) for v in t3f]
    print(f"TIER-3 characterization ({len(t3f)} families):")
    print(f"  loose: mean M={sum(Ms)/len(Ms):.3f} min {min(Ms):.3f} (tight 0.071); max-divisor-count values {sorted(set(md))} (boundary 6)")
    print(f"  => tier-3 = LOOSE, prime-rich, non-detuned, some-divisor-hits-exactly-6 families; near-unit-detected.")

if __name__=='__main__': main()
