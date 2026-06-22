#!/usr/bin/env python3
"""Proof push: test the OVER-DETERMINATION crux (mac-mini HYP-2878).
For covering 13-sets, count primes p (15<=p<=P) where S 'covers Z/p' (N(S,p)=0, no witness).
If robustly bounded << #primes, a witness prime exists => prime-certification.
kind-mendel-2026-06-22-S6."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def primes_upto(N):
    sieve=[True]*(N+1); sieve[0]=sieve[1]=False
    for i in range(2,int(N**.5)+1):
        if sieve[i]:
            for j in range(i*i,N+1,i): sieve[j]=False
    return [i for i in range(2,N+1) if sieve[i]]
def has_witness_modp(S,p):
    "True if exists a coprime to p with ||s a/p||>=1/14 for all s"
    for a in range(1,p):
        if gcd(a,p)!=1: continue
        if all(nrm(F(s*a,p))>=F(1,14) for s in S): return True
    return False

P=120
primes=[p for p in primes_upto(P) if p>=15]
print(f"primes in [15,{P}]: {len(primes)} -> {primes}")
random.seed(2)
print("\n=== #primes COVERED (no witness) per primitive covering 13-set ===")
def covered_count(S):
    return sum(1 for p in primes if not has_witness_modp(S,p))
# random covering sets at various scales + the loosest + the (non-covering) AP for contrast
samples=[]
for hi in [60,500,5000]:
    t=0
    while t<15:
        S=sorted(random.sample(range(1,hi),13))
        if gall(S)==1 and is_covering(S): samples.append((f'cov<{hi}',S)); t+=1
samples.append(('loosest{1..11,13,84}',[1,2,3,4,5,6,7,8,9,10,11,13,84]))
mx=0
for tag,S in samples:
    c=covered_count(S); mx=max(mx,c)
    if tag.startswith('loosest') or c>=3:
        print(f"  {tag:20s} covered {c}/{len(primes)} primes  S={S}")
print(f"  ... max covered over {len(samples)} covering sets = {mx}/{len(primes)}")
# the AP {1..13} (non-covering, tight): should cover ~all primes>=15 (witness only at 14)
S=list(range(1,14))
print(f"\n  AP {{1..13}} (non-covering, M=1/14): covered {covered_count(S)}/{len(primes)} primes "
      f"(expected ~all: witness only at D=14, not these primes)")
print(f"\nINTERPRETATION: covering sets cover FEW primes (max {mx}); the AP covers nearly all.")
print("=> over-determination: a covering set is prime-certified (witness prime exists) IF max-covered < #primes.")
print("   BUT NOTE: for p > max speed, residues=speeds, so 'covers Z/p' for large p ~ continuous non-loneliness (=Node-3).")
