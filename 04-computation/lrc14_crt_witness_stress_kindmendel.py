#!/usr/bin/env python3
"""Swing 3 deep: is the witness denominator D bounded over ALL covering 13-sets?
(= a constructive/elementary route to LRC(14), refining HYP-2566 uniform looseness.)
kind-mendel-2026-06-22-S4."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def min_witness_D(S, Dmax=3000):
    for D in range(2,Dmax+1):
        for a in range(1,D//2+1):
            if gcd(a,D)!=1: continue
            if all(nrm(F(s*a,D))>=F(1,14) for s in S): return D
    return None
def M_exact(S):
    "exact max-min reach via arrangement vertices tau in (0,1/2]"
    cands=set()
    for i,si in enumerate(S):
        for k in range(1,si): cands.add(F(2*k-1,2*si))   # (2k-1)/(2si)
        for sj in S:
            for d in (si+sj, abs(si-sj)):
                if d: 
                    for k in range(1,d): cands.add(F(k,d))
    best=F(0)
    for t in cands:
        if 0<t<=F(1,2):
            m=min(nrm(s*t) for s in S)
            if m>best: best=m
    return best

random.seed(11)
print("=== max witness D and min M over adversarial covering 13-sets ===")
for label, gen in [
   ("speeds<=300", lambda: sorted(random.sample(range(1,300),13))),
   ("speeds<=3000", lambda: sorted(random.sample(range(1,3000),13))),
   ("speeds<=30000", lambda: sorted(random.sample(range(1,30000),13))),
   ("near-loosest {1..11,13,84*m}+perturb", None),
]:
    if gen is None:
        # adversarial: perturbations of the loosest family
        sets=[]
        for m in [1,2,3,5,10]:
            base=[1,2,3,4,5,6,7,8,9,10,11,13,84*m]
            sets.append(sorted(set(base)))
        for extra in range(20):
            S=[1,2,3,4,5,6,7,8,9,10,11,13, 84*random.randint(1,50)]
            S=sorted(set(S)); 
            if len(S)==13: sets.append(S)
        maxD=0; minM=F(1); worstD=None
        for S in sets:
            if len(S)!=13 or gall(S)!=1 or not is_covering(S): continue
            D=min_witness_D(S, Dmax=3000); 
            if D and D>maxD: maxD=D; worstD=S
            M=M_exact(S); 
            if M<minM: minM=M
        print(f"{label:38s}: max D={maxD} (at {worstD}); min M={float(minM):.5f} (7/89={7/89:.5f})")
        continue
    tot=0; maxD=0; minM=F(1); trials=0; worstD=None; worstM=None; nofind=0
    while tot<50 and trials<300000:
        trials+=1; S=gen()
        if gall(S)!=1 or not is_covering(S): continue
        tot+=1
        D=min_witness_D(S, Dmax=3000)
        if D is None: nofind+=1; continue
        if D>maxD: maxD=D; worstD=S
        M=M_exact(S)
        if M<minM: minM=M; worstM=S
    print(f"{label:38s}: n={tot} maxD={maxD} (D<=3000 fails:{nofind}); minM={float(minM):.5f}")
    if worstD: print(f"      worst-D set: {worstD}")
print(f"\nDirichlet heuristic: D <~ 1/(M-1/14). loosest M=7/89 => 1/(7/89-1/14)={1/(7/89-1/14):.1f}")
