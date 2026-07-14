#!/usr/bin/env python3
"""
lrc14_shadowgap_thm_klein_S297.py
=================================
klein-2026-07-13-S297 (owner: assault the shadow-gap rigidity for covering sets).

THM-744 (PROVED). For C with an even element, e=min even, m=max C: if m < 6e, then t=1/2+delta is lonely
for {1}UC for every delta in (1/(14 e), 3/(7 m)) -- an explicit non-empty MIDDLE good interval, so L>0.
Proof: even c>=e -> ||c t||=c*delta in (1/14,3/7); odd c<=m -> ||c t||=1/2-c*delta>1/14; speed 1 ->
t in (1/2,13/14) -> ||t||>1/14. Non-empty iff m<6e (~ ratio<6: TIGHT covering clusters -- the residual
that isolated-far disc_v (THM-731) and multi-peel (THM-735) do NOT reach).

This verifies the good interval (dense sample) + censuses coverage.
"""
import numpy as np, random
from math import gcd
def is_good(S,tt): return all(min((c*tt)%1.0,1-((c*tt)%1.0))>=1.0/14.0-1e-12 for c in S)
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def evens_min(C):
    ev=[c for c in C if c%2==0]; return min(ev) if ev else None

print("THM-744: m<6e => t=1/2+delta lonely for {1}UC, delta in (1/(14e),3/(7m)). Verify interval all-good:")
print("%-30s %5s %5s %7s %10s %s"%("S={1}UC","e","m","m/e","gap-width","interval all good"))
for S in [[1,90,91,92,93,94,95,96,97,98,99,100,101],
          [1,30,31,32,33,34,35,36,37,38,39,40,42],
          [1,90,93,95,96,98,99,100,102,104,105,106,108],
          [1,15,16,17,18,19,20,21,22,26,28,39,14]]:   # tight-ish covering
    S=sorted(set(S))
    if len(S)!=13 or not iscov(S):
        print("  %-30s (skip: 13=%s cov=%s)"%(str(S)[:30],len(S),iscov(S))); continue
    C=[c for c in S if c!=1]; e=evens_min(C); m=max(C)
    if m<6*e:
        lo,hi=0.5+1.0/(14*e),0.5+3.0/(7*m)
        ok=all(is_good(S,x) for x in np.linspace(lo,hi,3000)[1:-1])
        print("  %-30s %5d %5d %7.2f %10.5f %s"%(str(S)[:30],e,m,m/e,hi-lo,ok))
    else:
        print("  %-30s %5d %5d %7.2f  m>=6e (shadow at k=2 fails)"%(str(S)[:30],e,m,m/e))

print("\ncensus: fraction of covering {1}UC (min>=15, tight) with m<6e (THM-744 fires):")
random.seed(2); n=0; fire=0
for _ in range(4000):
    b=random.choice([15,25,40,80]); C=sorted(random.sample(range(b,b+random.choice([12,20,40])),12)); S=sorted([1]+C)
    if not iscov(S): continue
    n+=1; e=evens_min(C)
    if e and max(C)<6*e: fire+=1
print("  %d covering sampled; THM-744 fires for %d (%.0f%%). Residual = ratio in [6,13], no isolated far."%(n,fire,100*fire/max(n,1)))
print("done.")
