#!/usr/bin/env python3
"""The Vitali connection, quantified. The danger of v=nw is a PERIODIC arc family
(period 1/(nw), radius 1/(n^2 w)); since G(S')=safe(S\{v}) is a finite union of
intervals, E ⊆ D_v iff every E-interval has length <= 2/(n^2 w) AND sits in a v-arc.
=> Criterion B' (max E-interval > 2/(n^2 w)) is an IFF-direction: it PROVES looseness
by Lebesgue density / Vitali covering. Quantify: of mult-of-n configs, how many are
loose via the PROVED long-interval criterion vs the residual (all-short alignment).
opus-2026-06-03-S573."""
from fractions import Fraction as F
from math import gcd
import random
def dist(x): x%=1; return min(x,1-x)
def G_intervals(Vp,n):  # components (start,len) of {t: ||v t||>1/n for all v in Vp}, exact
    THR=F(1,n); eps=set([F(0)])
    for v in Vp:
        for k in range(v+1):
            for s in(-1,1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); raw=[]; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>THR for v in Vp): raw.append((a,ln))
    return raw
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def loose(Vp,n):  # full config safe set nonempty-open?
    return bool(G_intervals(Vp,n))
def main():
    rng=random.Random(11)
    print("Mult-of-n configs: PROVED-by-Criterion-B' (max G(S\\v) interval > 2/(n^2 w)) vs residual")
    for n in [6,8,10,12,14]:
        m=n-1; tot=0; proved=0; resid_loose=0; resid_tight=0
        for _ in range(1200):
            others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            w=rng.randint(1,3); v=n*w
            if v in others: continue
            V=prim(tuple(sorted(others|{v})))
            if len(V)!=m or not any(x%n==0 for x in V): continue
            # recover the actual multiple in primitive coords
            mults=[x for x in V if x%n==0]
            if not mults: continue
            vv=mults[0]; ww=vv//n
            Sp=tuple(x for x in V if x!=vv)
            tot+=1
            gi=G_intervals(Sp,n)
            Lmax=max((ln for _,ln in gi), default=F(0))
            arc=F(2,n*n*ww)
            if Lmax>arc:
                proved+=1            # Criterion B' fires: PROVED loose (Vitali/density)
            else:
                if loose(V,n): resid_loose+=1   # all-short but still loose (alignment fails) = residual
                else: resid_tight+=1            # all-short AND covered => would be tight (measure 0)
        print(f"  n={n:2d}: {tot} mult-of-n configs; PROVED by B' (long interval)={proved} "
              f"({100*proved/max(tot,1):.1f}%); residual all-short loose={resid_loose}; "
              f"all-short tight(measure-0 w/ multiple)={resid_tight}")
if __name__=='__main__': main()
