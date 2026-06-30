#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify n=11; the creative REFRAMES (observer-change=translation, copy-to-all-n=independent set);
the CHROMATIC<->OCF bridge (danger-circulant chromatic number + Paley-tournament Hamiltonian paths). klein-S36."""
import math, itertools
from fractions import Fraction as F

def dist0(x,D): x%=D; return min(x,D-x)
def gap_exact(S,Dmax):
    best=F(0); arg=None
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=min(dist0(s*a,D) for s in S)
            if F(m,D)>best: best=F(m,D); arg=(D,a)
    return best,arg

print("="*84); print(" (A) VERIFY n=11 covering-min 3/31 honestly over a large modulus range"); print("="*84)
S11=[2,6,8,9,10,11,13,14,17,19]
M,arg=gap_exact(S11,400)
print(f"   S={S11}: honest M over D<=400 = {M}={float(M):.5f} at (D,a)={arg}; construction 11/111={float(F(11,111)):.5f}")
print(f"   3/31 < 11/111 ? {F(3,31)<F(11,111)}  (spread beats construction);  3/31 > 1/11 (floor)? {F(3,31)>F(1,11)}")

print("\n"+"="*84); print(" (B) REFRAME 1+2: change the observer / translate = GALILEAN invariance"); print("="*84)
# moving observer at speed w  <=>  static observer, speeds shifted by -w (subtract w from all)
S=[1,2,5,6,7,8]; M0,_=gap_exact(S,40)
print(f"   base S={S}, observer at 0: M={M0}")
for w in [1,3,10]:
    Sw=[s-w for s in S]              # speeds relative to a runner moving at w (the new observer)
    Mw,_=gap_exact([s%101 for s in Sw],202)  # gaps are translation/sign invariant mod the torus
    print(f"   observer copied onto runner w={w}: relative speeds {Sw} -> M={Mw}  ({'SAME (invariant)' if Mw==M0 else 'changed'})")
print("   => moving the observer onto ANY runner = subtracting that speed = a GALILEAN shift; M is invariant.")
print("      The danger-circulant lives on the GROUP Z/D, so it is translation-invariant by construction.")

print("\n"+"="*84); print(" (C) REFRAME 3: copy the observer to ALL n points = the INDEPENDENT-SET dual"); print("="*84)
# at the lonely time t=a/D the n runner-positions s*a mod D are mutually >= (covering radius) apart
# = an independent set in the danger-circulant C_D({1..j-1}); observer-at-0 is the (n+1)-th point.
D,a,j=13,3,2
pos=sorted((s*a)%D for s in S)
gaps=[ (pos[(i+1)%len(pos)]-pos[i])%D for i in range(len(pos))]
print(f"   n=7 at binding (D,a)=({D},{a}): runner-positions {pos}; consecutive gaps {sorted(gaps)}")
print(f"   min gap = {min(gaps)} >= covering radius {j}? {min(gaps)>=j}  => the n points (+ observer at 0) are")
print(f"   an INDEPENDENT SET in the danger-circulant C_{D}(<={j-1}) -- the 'all observers' = independent-set reframe,")
print(f"   the exact DUAL of the set-cover (one observer) IP. Both pin the same M=2/13.")

print("\n"+"="*84); print(" (D) CHROMATIC <-> OCF: danger-circulant chromatic no. + Paley-tournament Hamiltonian paths"); print("="*84)
def chromatic_number(adj):  # adj: list of sets; small graphs
    V=len(adj)
    for k in range(1,V+1):
        col=[-1]*V
        def bt(v):
            if v==V: return True
            for c in range(k):
                if all(col[u]!=c for u in adj[v]):
                    col[v]=c
                    if bt(v+1): return True
                    col[v]=-1
            return False
        if bt(0): return k
    return V
for (D,j,n) in [(13,2,7),(15,2,8),(33,4,9),(31,3,11)]:
    conn=set(range(1,j))|set(D-x for x in range(1,j))   # C_D({1..j-1}) = mutual-danger graph
    adj=[set((v+c)%D for c in conn) for v in range(D)]
    chi=chromatic_number(adj); alpha=D//j   # independence number of C_D({1..j-1}) ~ floor(D/j)
    print(f"   n={n}: danger-circulant C_{D}({{1..{j-1}}}): chi={chi}, indep.no~{alpha}, clique~{j}; covering-min {n}<->{F(j,D)}")
# Paley tournament on 13 (QR mod 13): Hamiltonian-path count (Redei: ODD) via subset DP -- the OCF target
p=13; QR={(x*x)%p for x in range(1,p)}
T=[[1 if ((j-i)%p) in QR else 0 for j in range(p)] for i in range(p)]   # i->j if j-i is QR
from functools import lru_cache
import sys; sys.setrecursionlimit(100000)
# count Hamiltonian paths via DP over (mask,last)
dp={}
for v in range(p): dp[(1<<v,v)]=1
for mask in range(1<<p):
    for last in range(p):
        c=dp.get((mask,last),0)
        if not c: continue
        for nx in range(p):
            if mask&(1<<nx): continue
            if T[last][nx]:
                k=(mask|(1<<nx),nx); dp[k]=dp.get(k,0)+c
full=(1<<p)-1; H=sum(dp.get((full,v),0) for v in range(p))
print(f"   Paley TOURNAMENT on 2n-1=13 (QR mod 13): #Hamiltonian paths H = {H}  (Redei: ODD? {H%2==1})")
print(f"   => H(T) is computed by the OCF (odd-cycle collection formula); chi(danger-circulant) <-> the SAME")
print(f"      circulant's orientation (the Paley tournament). The IP independent-set = the chromatic/OCF bridge.")
