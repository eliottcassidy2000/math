#!/usr/bin/env python3
# k=9: how tight is the cap margin? Scan all "near-consec" perturbations exactly.
# Compare L_y(competitor) to L_y(consec) AND to cap. Find the WORST (closest-to-cap, closest-to-consec) competitors.
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

def Nat(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi<=lo: continue
        t=Nat(E,(lo+hi)/2); p[t]+=hi-lo
    return p
def gpoly(k): return [F(-(t-2)*(t-3)*(t-6),36) for t in range(7)]  # k=9,10
def Ly(E,k):
    p=dist(E); g=gpoly(k); return sum(p[t]*g[t] for t in range(7)), p
def primitive(E): return reduce(gcd,E)==1

k=9; cap=F(49426,100000)  # approx; use float compare
C=list(range(k)); LC,pC=Ly(C,k)
print(f"consec L_y = {LC} = {float(LC):.6f}; cap ~ 0.49426; margin={0.49426-float(LC):.6f}")
print(f"consec p = {[float(x) for x in pC]}")

# Exhaustive over all 9-subsets of {0..N} containing 0, primitive, max<=14 (matches verified region)
N=14
beats=[]; ties=[]; closest=[]
cnt=0
for sub in itertools.combinations(range(1,N+1),k-1):
    E=[0]+list(sub)
    if not primitive(E): continue
    cnt+=1
    L,_=Ly(E,k)
    d=LC-L
    if d<0: beats.append((E,L))
    elif d==0: ties.append(E)
    closest.append((d,E,L))
closest.sort()
print(f"\nscanned {cnt} primitive 9-sets, max<=14")
print(f"beats consec: {len(beats)}  exact ties: {len(ties)}")
print("ties (should be only dilations/none):", ties[:10])
print("\nclosest 8 competitors to consec (smallest L_y gap), with cap margin:")
for d,E,L in closest[1:9]:
    print(f"  gap-to-consec={float(d):.6f}  L_y={float(L):.6f}  cap-margin={0.49426-float(L):.6f}  E={E}")
