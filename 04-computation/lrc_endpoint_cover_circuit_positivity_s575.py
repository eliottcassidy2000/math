#!/usr/bin/env python3
"""Endpoint-cover CIRCUIT POSITIVITY. Component C_i=(a_i,b_i) of G(S') is coverable by
a v=nw arc iff ||nw*m_i|| <= 1/n - (nw/2)*l_i  (m_i midpoint, l_i length). Tight <=>
ALL components coverable by one integer nw. Circuit-positivity margin
   P(S) = max_i ( ||nw*m_i|| - (1/n - (nw/2) l_i) );  P>0  <=>  loose.
Verify the criterion (100% vs direct), and that P>0 for every mult-of-n config
(the circuit can never close), with the per-component winding tie. opus-2026-06-03-S575."""
from fractions import Fraction as F
from math import gcd
import random
def dist(x): x%=1; return min(x,1-x)
def G_components(Sp,n):
    THR=F(1,n); eps=set([F(0)])
    for u in Sp:
        for k in range(u+1):
            for s in(-1,1): eps.add(F(k*n+s,n*u)%1)
    pts=sorted(eps); comps=[]; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(u*mid)>THR for u in Sp): comps.append((a,b,ln))
    return comps
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def tight_direct(V,n):
    return not G_components(V,n)
def coverable(a,b,ln,n,v):
    # component (a,b) of G(S') fits one v-arc  <=>  ||v*m|| <= 1/n - (v/2)*ln
    # where v = nw is the SPEED (the multiple of n), m the midpoint, ln the length.
    m=(a+ln/2)%1  # true midpoint, handles the wrap-around component
    lhs=dist(v*m); rhs=F(1,n)-F(v,2)*ln
    return lhs<=rhs, lhs-rhs
def main():
    rng=random.Random(34)
    print("Circuit-positivity: criterion ||nw m_i||<=1/n-(nw/2)l_i  vs direct tight/loose")
    for n in [6,8,10,12,14]:
        m=n-1; tot=0; agree=0; Ppos=0; tight=0
        for _ in range(2500):
            others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            w=rng.randint(1,3); v=n*w
            if v in others: continue
            V=prim(tuple(sorted(others|{v})))
            if len(V)!=m or not any(x%n==0 for x in V): continue
            mults=[x for x in V if x%n==0]; vv=mults[0]; ww=vv//n
            Sp=tuple(x for x in V if x!=vv)
            comps=G_components(Sp,n)
            if not comps: continue
            tot+=1
            # circuit positivity margin over components (loose iff some component uncoverable)
            margins=[]
            allcov=True
            for (a,b,ln) in comps:
                cov,marg=coverable(a,b,ln,n,vv)
                margins.append(marg)
                if not cov: allcov=False
            P=max(margins)
            pred_tight=allcov
            actual=tight_direct(V,n)
            if pred_tight==actual: agree+=1
            if P>0: Ppos+=1
            if actual: tight+=1
        print(f"  n={n:2d}: {tot} mult-of-n; criterion agrees={agree}/{tot}; "
              f"circuit-positive (P>0 => loose)={Ppos}/{tot}; actually tight={tight}")
if __name__=='__main__': main()
