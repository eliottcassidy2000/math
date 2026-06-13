#!/usr/bin/env python3
"""Sleeve saturation. Each runner v is a SLEEVE Sigma_v={t:||vt||<delta}, measure 2*delta.
Total sleeve U=union; safe=complement; deficit Delta(S)=mu(safe). measure-saturated:
mu(safe)=0 (worry-set). point-saturated: safe=empty (M<delta, counterexample; LRC: never).
Explore the SATURATION CURVE mu_j = mu(safe of first j sleeves at fixed level delta=1/(k+1)),
the recursion, and the spread-vs-overlap tension. opus-2026-06-03-S578."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def safe_measure(V,delta):
    eps=set([F(0)])
    for v in V:
        # boundary pts ||v t||=delta: t=(k+-delta)/v
        kmax=v
        for k in range(kmax+1):
            for s in(1,-1):
                eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); meas=F(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>=delta for v in V): meas+=ln
    return meas
def three_terms(V):
    s=set(V); return sum(1 for a,b in combinations(sorted(V),2) if a+b in s)
def sat_curve(V,delta):  # mu_j for j=1..k, runners added in given order
    return [float(safe_measure(V[:j],delta)) for j in range(1,len(V)+1)]
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    print("Saturation curves mu_j (safe measure of first j sleeves, level delta=1/(k+1))")
    for k in [6,8,10,12]:
        delta=F(1,k+1)
        AP=tuple(range(1,k+1))
        c=sat_curve(AP,delta)
        sat_idx=next((j+1 for j,m in enumerate(c) if m==0), None)
        print(f"  k={k:2d} AP: curve={['%.3f'%x for x in c]}; measure-saturates at sleeve #{sat_idx}")
    print()
    print("Spread-vs-overlap tension: deficit mu(safe) vs #3-term (overlap), random configs")
    rng=random.Random(7)
    for k in [8,10]:
        delta=F(1,k+1); B=2*k+6
        rows=[]
        for _ in range(4000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            rows.append((three_terms(V), float(safe_measure(V,delta))))
        # bin by 3-term
        import statistics
        bins={}
        for t3,mu in rows: bins.setdefault(min(t3,4),[]).append(mu)
        print(f"  k={k}: "+"; ".join(f"#3t={b}{'+' if b==4 else ''}: meanDeficit={statistics.mean(v):.4f},min={min(v):.4f}" for b,v in sorted(bins.items())))
if __name__=='__main__': main()
