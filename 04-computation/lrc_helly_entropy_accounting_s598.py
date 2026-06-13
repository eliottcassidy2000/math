#!/usr/bin/env python3
"""Helly entropy accounting: the cascade |SAFE|=prod c_i (S545), its log = CLEARANCE
ENTROPY (additive); the H-entropy S_H=log H (S543, max at tight); the Helly number of
the looseness certificate (S599); iterated-log scale currency (S600). Extend: the
clearance-entropy / H-entropy DUALITY (both extremized at the tight witness) + the
bounded-Helly trapping. opus-2026-06-03-S598."""
from fractions import Fraction as F
from math import log2, log
from itertools import permutations, combinations
def dist(x): x%=1; return min(x,1-x)
def safe_measure_prefix(V,n,j):  # measure of safe set of first j runners at level 1/n
    THR=F(1,n); Vp=V[:j]; 
    if not Vp: return F(1)
    eps=set([F(0)])
    for v in Vp:
        for k in range(v+1):
            for s in(1,-1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); meas=F(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>=THR for v in Vp): meas+=ln
    return meas
def cascade(V,n):  # c_i = mu_i/mu_{i-1}; |SAFE|=prod; returns clearances and entropy
    mus=[safe_measure_prefix(V,n,j) for j in range(len(V)+1)]
    cs=[]
    for i in range(1,len(mus)):
        cs.append(mus[i]/mus[i-1] if mus[i-1]!=0 else F(0))
    return cs, mus[-1]
def helly_trap(V,n):  # min subfamily size s.t. their joint safe set is measure-0 (the looseness/tightness certificate)
    THR=F(1,n)
    for r in range(1,len(V)+1):
        for sub in combinations(range(len(V)),r):
            Vs=[V[i] for i in sub]
            m=safe_measure_prefix(tuple(Vs)+tuple(),n,len(Vs))
            if m==0: return r,Vs
    return len(V),list(V)
def main():
    print("CLEARANCE-ENTROPY ledger log|SAFE|=sum log c_i; and the Helly trap number")
    for name,V,n in [("AP n=6 tight",(1,2,3,4,5),6),("generic n=6 loose",(1,2,4,7,9),6),
                     ("AP n=7 tight",(1,2,3,4,5,6),7),("generic n=7",(1,3,5,8,11,12),7)]:
        cs,SAFE=cascade(V,n)
        ent=[(float(log2(c)) if c>0 else float('-inf')) for c in cs]
        r,sub=helly_trap(V,n) if SAFE==0 else (None,None)
        print(f"  {name}: |SAFE|={float(SAFE):.4f}; clearances={[str(c) for c in cs]}")
        print(f"     clearance-entropy log2 c_i = {[('%.2f'%e if e!=float('-inf') else '-inf') for e in ent]}; sum=log2|SAFE|={('%.3f'%log2(float(SAFE)) if SAFE>0 else '-inf')}")
        if SAFE==0: print(f"     HELLY trap number (min subfamily with measure-0 safe set) = {r}, runners {sub}")
if __name__=='__main__': main()
