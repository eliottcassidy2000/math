#!/usr/bin/env python3
"""New PROVED criterion from the translator: a G(S') component with BOTH endpoints
owned by speeds < n is uncoverable by any v=nw arc (small owner => endpoint = arc
center => both centers equal => a=b, impossible). So such a component => S loose for
EVERY w. Measure coverage (alone, and combined with Criterion B' long-interval)."""
from fractions import Fraction as F
from math import gcd
import random
def dist(x): x%=1; return min(x,1-x)
def G_components(Sp,n):
    THR=F(1,n); pts={}
    for u in Sp:
        for k in range(u+1):
            for eps in (1,-1):
                pts.setdefault(F(k*n+eps,n*u)%1,[]).append((u,k,eps))
    order=sorted(pts); comps=[]; L=len(order)
    for i in range(L):
        a=order[i]; b=order[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(u*mid)>THR for u in Sp): comps.append((a,b,ln,pts[a],pts[b]))
    return comps
def small_both(comps,n):
    for (a,b,ln,oa,ob) in comps:
        la=[o for o in oa if o[2]==1]; rb=[o for o in ob if o[2]==-1]
        if la and rb and min(o[0] for o in la)<n and min(o[0] for o in rb)<n:
            return True
    return False
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(57)
    print("Coverage of C' by PROVED criteria: small-owner-component (w-free) and B' (long interval)")
    for n in [6,8,10,12,14]:
        m=n-1; tot=0; sm=0; comb=0; resid=[]
        for _ in range(2000):
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
            sb=small_both(comps,n)
            Lmax=max((ln for _,_,ln,_,_ in comps),default=F(0)); longi=Lmax>F(2,n*n*ww)
            if sb: sm+=1
            if sb or longi: comb+=1
            else: resid.append((V,float(Lmax)))
        print(f"  n={n:2d}: {tot} mult-of-n; small-owner-component={sm} ({100*sm/max(tot,1):.1f}%); "
              f"small-owner OR B' = {comb} ({100*comb/max(tot,1):.1f}%); residual={len(resid)}"
              +(f"  e.g.{resid[0][0]}" if resid else ""))
if __name__=='__main__': main()
