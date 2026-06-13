#!/usr/bin/env python3
"""The real dichotomy is the AUGMENTATION index j=eps(c)=sum of coeffs, NOT support.
unbalanced (j!=0, observer-coupled): a+b=c (eps=1), c=2a (eps=1). balanced (j=0,
observer-blind, translation-invariant): 2a=b+c (eps=0, support-3 AP-triple!), a+b=c+d
(eps=0). Verify hardness tracks UNBALANCED count, is INDEPENDENT of balanced count.
Theorem: balanced <=> translation-invariant.  opus-2026-06-03-S581."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random, statistics
def dist(x): x%=1; return min(x,1-x)
def Mexact(V):
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (V[i]+V[j],abs(V[i]-V[j])):
                if D:
                    for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*V[i]))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def unbalanced(V):  # eps!=0 small relations: a+b=c and c=2a
    s=set(V); u=0
    for a,b in combinations(sorted(V),2):
        if a+b in s: u+=1          # a+b=c, eps=1
    for a in V:
        if 2*a in s: u+=1          # c=2a, eps=1
    return u
def balanced(V):    # eps=0 relations: 2a=b+c (AP-triple) and a+b=c+d
    s=set(V); bcount=0
    V=sorted(V)
    for a in V:                    # 2a=b+c, b<a<c
        for b in V:
            c=2*a-b
            if b<a<c and c in s: bcount+=1
    sums={}
    for a,b in combinations(V,2): sums.setdefault(a+b,0); sums[a+b]+=1
    bcount+= sum(x*(x-1)//2 for x in sums.values())   # a+b=c+d
    return bcount
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(11)
    print("THEOREM check: relation survives translation S->S+t  <=>  eps=0 (balanced)")
    for rel,coeffs in [("a+b=c",(1,1,-1)),("2a=b+c",(2,-1,-1)),("a+b=c+d",(1,1,-1,-1)),("c=2a",(-1,2))]:
        print(f"   {rel:10s}: eps={sum(coeffs):+d}  translation-invariant={sum(coeffs)==0}")
    print()
    print("STRATIFICATION: M-delta binned by UNBALANCED count and (separately) by BALANCED count")
    for k in [6,8]:
        delta=F(1,k+1); B=2*k+8
        byU={}; byB_givenU0={}
        for _ in range(9000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            u=unbalanced(V); bb=balanced(V); M=float(Mexact(V)-delta)
            byU.setdefault(min(u,3),[]).append(M)
            if u==0:                              # ISOLATE: no unbalanced => vary balanced
                byB_givenU0.setdefault(min(bb,3),[]).append(M)
        print(f"  k={k} by UNBALANCED count: "+"; ".join(f"u={u}{'+'if u==3 else''}:minMargin={min(v):+.4f}(n={len(v)})" for u,v in sorted(byU.items())))
        if byB_givenU0:
            print(f"  k={k} | u=0, by BALANCED count: "+"; ".join(f"b={b}{'+'if b==3 else''}:minMargin={min(v):+.4f}(n={len(v)})" for b,v in sorted(byB_givenU0.items())))
if __name__=='__main__': main()
