#!/usr/bin/env python3
"""The circuit-to-gap functional Phi(C): the EXACT refinement of ECCP positivity P(S).
P(S) read only the sign (max term). Phi sums per-component uncovered phase:
each component C_i=(a_i,b_i) of G(S') maps in v-phase to (v a_i, v b_i); v-danger is the
band U_k (k-1/n, k+1/n); uncovered_i = 2h_i - sum_k |(v a_i,v b_i) cap band_k|;
   Phi(C) = (1/v) sum_i uncovered_i.
CLAIM (verify exactly):  G(v) := mu(safe set of S=S'u{v}) = Phi(C).  Kernel{Phi=0}=tight.
opus-2026-06-03-S576."""
from fractions import Fraction as F
from math import gcd, floor, ceil
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
        if all(dist(u*mid)>THR for u in Sp): comps.append((a,ln))  # (lo, length) on the line
    return comps
def safe_measure(V,n):
    THR=F(1,n); eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(-1,1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); meas=F(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>=THR for v in V): meas+=ln
    return meas
def interval_band_overlap(lo,hi,n):
    """sum over integers k of |(lo,hi) cap (k-1/n, k+1/n)|, lo<hi reals (Fractions)."""
    inv=F(1,n); tot=F(0)
    k0=floor(lo-inv); k1=ceil(hi+inv)
    for k in range(k0,k1+1):
        l=max(lo,k-inv); r=min(hi,k+inv)
        if r>l: tot+=r-l
    return tot
def Phi(Sp,v,n):
    """circuit-to-gap functional over the components of G(Sp), arc speed v."""
    tot=F(0)
    for (lo,ln) in G_components(Sp,n):
        a=lo; b=lo+ln
        vph_lo=v*a; vph_hi=v*b           # phase interval (v a, v b), length v*ln
        ov=interval_band_overlap(vph_lo,vph_hi,n)
        uncovered=v*ln-ov                # = 2h - overlap
        tot+=uncovered
    return tot/v
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(76)
    print("Verify  G(v)=mu(safe set) == Phi(C)  exactly, over multiple-of-n configs")
    for n in [6,8,10,12,14]:
        m=n-1; tot=0; ok=0; maxerr=F(0)
        for _ in range(900):
            others=set(rng.sample([x for x in range(1,n+8) if x%n!=0],m-1))
            w=rng.randint(1,3); v=n*w
            if v in others: continue
            V=prim(tuple(sorted(others|{v})))
            if len(V)!=m or not any(x%n==0 for x in V): continue
            mults=[x for x in V if x%n==0]; vv=mults[0]
            Sp=tuple(x for x in V if x!=vv)
            if not G_components(Sp,n): continue
            tot+=1
            direct=safe_measure(V,n); phi=Phi(Sp,vv,n)
            if direct==phi: ok+=1
            else: maxerr=max(maxerr,abs(direct-phi))
        print(f"  n={n:2d}: {tot} configs; Phi==mu(safe) exactly: {ok}/{tot}; max|err|={float(maxerr):.2e}")
if __name__=='__main__': main()
