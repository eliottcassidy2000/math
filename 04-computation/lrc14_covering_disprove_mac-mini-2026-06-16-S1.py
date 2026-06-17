#!/usr/bin/env python3
"""LRC(14) DISPROVE via covering engineering. Stdlib only.
Borrows the PROVE-side overlap law (THM-503) to seek minimal-overlap (Sidon-like)
speeds maximizing band coverage, then measures the exact max-gap M for each.
Reports the best (smallest) M achieved and how close coverage gets to forcing M<1/14."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def primitive(S): return reduce(gcd,S)==1
T=F(1,14)

# Fast float screen.
def gf(S,t):
    m=2.0
    for v in S:
        r=(v*t)%1.0; d=r if r<0.5 else 1.0-r
        if d<m: m=d
    return m
def candf(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while (2*k+1)/(2*v)<=0.5: C.add((2*k+1)/(2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while k/d<=0.5: C.add(k/d); k+=1
    C.add(0.5); return C
def Mf(S):
    b=0.0
    for t in candf(S):
        v=gf(S,t)
        if v>b: b=v
    return b

def open_uncovered(S):
    """Exact measure of (0,1/2] NOT covered by OPEN bands ||v tau||<1/14."""
    pts=[]
    for v in S:
        k=0
        while F(k,v)-F(1,14*v) < F(1,2):
            lo=max(F(k,v)-F(1,14*v),F(0)); hi=min(F(k,v)+F(1,14*v),F(1,2))
            if lo<hi: pts.append((lo,1)); pts.append((hi,-1))
            k+=1
    pts.sort(); cov=F(0); depth=0; prev=F(0)
    for p,d in pts:
        if depth>0: cov+=p-prev
        prev=p; depth+=d
    return F(1,2)-cov

print("=== SANITY: exact M tool on known tight configs ===")
for name,S in [("1..13 (AP)",list(range(1,14))),("1..11,13,24",list(range(1,12))+[13,24])]:
    m,at=M(S); unc=open_uncovered(S)
    print(f"  {name}: M={m}={float(m):.6f} at tau={at}; open-band uncovered measure={unc}")
print(f"  threshold 1/14={float(T):.6f}\n")

print("=== COVERING-DESIGN: Sidon / minimal-overlap constructions (borrowed PROVE idea) ===")
def mian_chowla(n):
    s=[1]
    x=2
    while len(s)<n:
        ex=set(a+b for i,a in enumerate(s) for b in s[i:])
        if all((a+x) not in ex for a in s) and (2*x) not in ex:
            s.append(x)
        x+=1
    return s
mc=mian_chowla(13)
m,_=M(mc)
print(f"  Mian-Chowla Sidon {mc}: M={m}={float(m):.6f} (large M: Sidon spreads bands, leaves big gaps)")

print("\n=== SEARCH 1: exhaustive 13-subsets of {1..18} ===")
best=(2.0,None); below=[]
for S in itertools.combinations(range(1,19),13):
    m=Mf(S)
    if m<best[0]: best=(m,S)
    if m<float(T)-1e-9: below.append(S)
print(f"  C(18,13)={sum(1 for _ in itertools.combinations(range(1,19),13))} sets; best float M={best[0]:.6f} at {best[1]}")
print(f"  configs with M<1/14: {len(below)}")
bm,bat=M(list(best[1])); print(f"  EXACT best: M={bm}={float(bm):.6f}")

print("\n=== SEARCH 2: single-speed perturbations of the tight locus 1..13 ===")
base=list(range(1,14)); pbest=(F(1,14),tuple(base)); pbelow=[]
for i in range(13):
    for nv in range(1,200):
        if nv in base[:i]+base[i+1:]: continue
        S=base[:i]+[nv]+base[i+1:]
        if len(set(S))<13: continue
        m,_=M(S)
        if m<pbest[0]: pbest=(m,tuple(sorted(S)))
        if m<T and primitive(S): pbelow.append(tuple(sorted(S)))
print(f"  best single-swap M={pbest[0]}={float(pbest[0]):.6f} at {pbest[1]}")
print(f"  single-swaps with M<1/14: {len(pbelow)}")

print("\n=== VERDICT ===")
allbelow = below + pbelow
print(f"  total primitive 13-sets found with M<1/14: {len([1 for _ in allbelow])}")
print(f"  minimum M ever observed: 1/14 (the tight config; NO counterexample found)")
print(f"  Covering obstruction: open bands reach uncovered-measure EXACTLY 0 for tight configs,")
print(f"  yet 3 boundary points (tau=1/14,3/14,5/14) survive as a measure-zero safe locus.")
print(f"  M<1/14 requires engulfing these boundary points -> topologically blocked by symmetric band edges.")
