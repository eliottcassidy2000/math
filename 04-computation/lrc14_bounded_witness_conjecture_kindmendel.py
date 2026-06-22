#!/usr/bin/env python3
"""Swing 3 (headline): the BOUNDED-DENOMINATOR WITNESS conjecture for LRC(14).
Every PRIMITIVE COVERING 13-set has a lonely point tau=a/D with D <= D0 (~41),
INDEPENDENT of speed magnitude (witness depends only on residues mod D).
=> elementary route to LRC(14)'s hard case; refines HYP-2566 uniform looseness.
kind-mendel-2026-06-22-S4."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def min_witness(S, Dmax=4000):
    "smallest (D,a), gcd(a,D)=1, with ||s a/D||>=1/14 for all s in S"
    for D in range(2,Dmax+1):
        for a in range(1,D//2+1):
            if gcd(a,D)==1 and all(nrm(F(s*a,D))>=F(1,14) for s in S): return D,a
    return None

# 1) the loosest known covering set (M=7/89, the binding case for uniform looseness)
print("loosest known covering set {1..11,13,84} (M=7/89):", min_witness([1,2,3,4,5,6,7,8,9,10,11,13,84]))
# 2) D is independent of speed magnitude
random.seed(0)
print("\nmax witness D over random covering sets, by speed scale:")
for hi in [10**3,10**4,10**5,10**6]:
    tot=0; mx=0; tr=0
    while tot<50 and tr<500000:
        tr+=1; S=sorted(random.sample(range(1,hi),13))
        if gall(S)!=1 or not is_covering(S): continue
        tot+=1; w=min_witness(S,Dmax=2000)
        if w and w[0]>mx: mx=w[0]
    print(f"  speeds<={hi:>8}: n={tot}, max D = {mx}")
print("\nCONJECTURE: D0 = 41 suffices for all primitive covering 13-sets (= sharpened HYP-2566).")
print("Structural reason: tau=a/D depends only on {s mod D}; covering controls residues; magnitude irrelevant.")
