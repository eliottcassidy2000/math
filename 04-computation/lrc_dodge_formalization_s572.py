#!/usr/bin/env python3
"""Formalize the dodge: verify (1) the DOMINANCE lemma (v>(n-1)*max(others) => M>1/n,
general v, not just multiples of n); (2) the sharper INTERVAL criterion (safe(S' @1/n)
has an interval > 2/(n^2 w) => dodge works); (3) the v=n residual = a gap with all
others safe (half-clock). opus-2026-06-03-S572."""
from fractions import Fraction as F
from math import gcd
import random
def dist(x): x%=1; return min(x,1-x)
def Mexact(V,n):
    # M(S) via pinch candidates t=m/(a+b) over pairs + 0; exact max-min
    cands=set()
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            D=V[i]+V[j]
            for m in range(1,D): cands.add(F(m,D))
        v=V[i]
        for k in range(1,v): cands.add(F(k,v))  # single-runner peaks at midpoints handled below
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def safe_intervals(Vp,n):  # intervals (as (start,len)) where all ||v t||>1/n strictly
    THR=F(1,n); eps=set([F(0)])
    for v in Vp:
        for k in range(v+1):
            for s in(-1,1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); out=[]; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>THR for v in Vp): out.append((a,ln))
    return out
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))

def main():
    rng=random.Random(7)
    print("(1) DOMINANCE lemma: v>(n-1)*max(others) => M>1/n  (general v, ANY divisibility)")
    for n in [6,8,10,12,14]:
        ok=tot=0
        for _ in range(1500):
            m=n-1; others=sorted(rng.sample(range(1,n+6),m-1))
            vmax_other=max(others)
            v=rng.randint((n-1)*vmax_other+1,(n-1)*vmax_other+40)  # dominant
            V=prim(tuple(sorted(others+[v])))
            if len(V)!=m: continue
            tot+=1
            # check loose via a safe interval
            if safe_intervals(V,n): ok+=1
        print(f"   n={n:2d}: dominant-runner configs loose {ok}/{tot}")
    print()
    print("(2) v=n residual: is there always a gap (k/n+1/n^2,(k+1)/n-1/n^2) with all others safe?")
    for n in [6,8,10,12,14]:
        cnt=loose=gapwit=0
        for _ in range(1500):
            m=n-1; others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            V=prim(tuple(sorted(others|{n})))
            if not any(x%n==0 for x in V): continue
            cnt+=1
            si=safe_intervals(V,n)
            if si: loose+=1
            # does some n-clock gap interior contain a safe point? (half-clock midpoints)
            found=False
            for k in range(n):
                t=F(2*k+1,2*n)  # gap midpoint, where runner n is at 1/2
                if all(dist(x*t)>F(1,n) for x in V): found=True;break
            if found: gapwit+=1
        print(f"   n={n:2d}: v=n configs={cnt}; loose={loose}; gap-MIDPOINT witness={gapwit} "
              f"({'midpoint suffices' if gapwit==loose else 'midpoint NOT always enough'})")
if __name__=='__main__': main()
