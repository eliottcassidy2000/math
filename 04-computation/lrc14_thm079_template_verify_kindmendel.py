#!/usr/bin/env python3
"""Flush-out verification: LRC(14) = THM-079 template (reduce-to-atom + Moon forcing).
Verifies the Move-B linchpins: tight locus {AP,GW} both NON-covering with denom-14 optima;
bounded covering core has M>1/14; apex-7 floor excludes covering from denom-14.
kind-mendel-2026-06-22-S8."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def M_argmax(S):
    cands=set()
    for si in S:
        for k in range(1,si): cands.add(F(2*k-1,2*si))
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for k in range(1,d): cands.add(F(k,d))
    best=F(0); args=[]
    for t in cands:
        if 0<t<=F(1,2):
            m=min(nrm(s*t) for s in S)
            if m>best: best=m; args=[t]
            elif m==best: args.append(t)
    return best,args

print("MOVE B linchpins (the LRC Moon/forcing step):")
print("1) tight locus {AP, GW}: M=1/14, optimum at denom-14, NON-covering")
for name,S in [('AP {1..13}',list(range(1,14))),('GW {1..11,13,24}',[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    M,args=M_argmax(S); dens=sorted({t.denominator for t in args})
    print(f"   {name}: M={M} optimum-denoms={dens} covering={is_covering(S)}")
print("2) apex-7 floor: covering forces a mult of 14 -> sits on observer at every t=a/14")
print("   (||14k * a/14|| = 0 < 1/14 for all a) => covering excluded from denom-14 optima")
print("3) bounded covering core has M>1/14 (the value):")
random.seed(0); worst=(F(1),None); cnt=0
for _ in range(30000):
    S=sorted(random.sample(range(1,23),13))
    if gall(S)!=1 or not is_covering(S): continue
    cnt+=1; M,_=M_argmax(S)
    if M<worst[0]: worst=(M,S)
print(f"   {cnt} bounded covering sets in [1,22]: min M = {worst[0]}={float(worst[0]):.5f} > 1/14={1/14:.5f}: {worst[0]>F(1,14)}")
print("\nCRUX (*): tight (M=1/14) => optimum at apex-blocked (denom-14) point => covering forces M>1/14.")
print("VERIFIED for {AP,GW}; (*) in general = the three-gap tight-locus characterization = the open Moon step.")
