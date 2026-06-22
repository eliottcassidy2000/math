#!/usr/bin/env python3
"""Kuratowski/Robertson-Seymour angle on HYP-2864: is there a FINITE set of
denominators D such that every primitive covering 13-set has a witness a/D for some D
in the set? (= a finite 'certificate basis', the RS finite-obstruction analogue.)
kind-mendel-2026-06-22-S5."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from collections import Counter
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def witnesses_for_D(S, D):
    "set of multipliers a (coprime to D, 1<=a<=D-1) with ||s a/D||>=1/14 for all s"
    out=[]
    for a in range(1,D):
        if gcd(a,D)!=1: continue
        if all(nrm(F(s*a,D))>=F(1,14) for s in S): out.append(a)
    return out
def min_witness_D(S, Dmax=200):
    for D in range(2,Dmax+1):
        if witnesses_for_D(S,D): return D
    return None

random.seed(0)
# collect a big bank of primitive covering 13-sets (varied scales)
bank=[]
for hi in [60,200,1000,10000]:
    tot=0; tr=0
    while tot<150 and tr<200000:
        tr+=1; S=sorted(random.sample(range(1,hi),13))
        if gall(S)==1 and is_covering(S): bank.append(S); tot+=1
# add the loosest known + structured
bank.append([1,2,3,4,5,6,7,8,9,10,11,13,84])
bank.append(list(range(2,15)))
print(f"bank: {len(bank)} primitive covering 13-sets")

# 1) which D is each set's MINIMAL witness denominator?
minDs=Counter()
for S in bank:
    d=min_witness_D(S, Dmax=200)
    minDs[d]+=1
print("\nminimal witness-D distribution (D : #sets):")
for d,c in sorted(minDs.items(), key=lambda x:(x[0] is None,x[0])): print(f"  D={d}: {c}")

# 2) GREEDY: find a small FIXED set of D's covering all sets (each set has a witness for some chosen D)
print("\n=== Greedy minimal denominator-basis (finite certificate set) ===")
remaining=list(range(len(bank)))
chosen=[]
candidateDs=list(range(2,100))
while remaining:
    # pick D covering the most remaining sets
    best=None
    for D in candidateDs:
        cov=[i for i in remaining if witnesses_for_D(bank[i],D)]
        if best is None or len(cov)>len(best[1]): best=(D,cov)
    D,cov=best
    if not cov:
        print(f"  STUCK: {len(remaining)} sets uncovered by any D<100"); break
    chosen.append(D); remaining=[i for i in remaining if i not in set(cov)]
    print(f"  +D={D} covers {len(cov)} more; remaining={len(remaining)}")
print(f"FINITE CERTIFICATE BASIS: {chosen}  (size {len(chosen)})")
print("Interpretation: every covering set in the bank is lonely via a/D for one of these few D.")
