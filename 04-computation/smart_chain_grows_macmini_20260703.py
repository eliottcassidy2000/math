#!/usr/bin/env python3
"""
SMART CHAIN GROWS (mac-mini-2026-07-03-S27): retraction of the S27 'chain is census-able' claim.
A geometric chain with lcm-PRODUCT runners at distinct scales (each blocking a GROUP of band moduli) has
witness q GROWING ~log(topN) -- q=43,59,89,97 at topN=10^7..10^15 -- SAME as the one-scale cluster. My earlier
'chain census-able (q<=46)' used a weak construction (one band modulus per runner). The free-q depends on the
lcm-GROUPS blocked (CRT capacity: 13 runners<=M block {15..Q}, Q~log M), NOT the runners' scale structure.
So the crux is the general lcm-blocker at any scale, q~3.6 ln M (= HYP-4040), not narrowed by scale.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd,xs)
def lcm(a,b): return a*b//gcd(a,b)
def nd(x): x=x%1; return min(x,1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def min_witness_q(sp,qmax=400):
    for q in range(2,qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q)>=1/14 for v in sp): return q
    return None

def smart_chain(topN, rng):
    mods=list(range(15,80)); rng.shuffle(mods)
    runners=[]; i=0; scale=topN
    while i<len(mods) and len(runners)<11 and scale>=15:
        grp=[]; L=1
        while i<len(mods):
            L2=lcm(L,mods[i])
            if L2<=scale and len(grp)<6: L=L2; grp.append(mods[i]); i+=1
            else: break
        if L>1:
            runners.append(L*max(1, round(scale/L))); scale=scale//8   # distinct scale, ratio-step 8
        else: i+=1
    for q in [8,9,5,7,11,13,2,3,4,6]:
        if len(runners)>=13: break
        if not any(s%q==0 for s in runners): runners.append(q)
    while len(runners)<13: runners.append(rng.randint(2,22))
    return sorted(set(runners))[:13]

if __name__ == "__main__":
    rng=random.Random(555)
    print("SMART chain (lcm-product runners, distinct scales): witness q vs topN")
    for topN in [10**7, 10**9, 10**12, 10**15]:
        worst=0; nc=0
        for _ in range(3000):
            sp=smart_chain(topN, rng)
            if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
            nc+=1
            q=min_witness_q(sp,400)
            if q and q>worst: worst=q
        print(f"  topN={topN:>16}: #cov={nc:>5} MAX witness q={worst}")
    print("=> grows ~log(topN); the chain is NOT uniformly census-able; free-q = CRT capacity, scale-independent")
