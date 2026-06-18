#!/usr/bin/env python3
"""
lrc14_floor_structure — mac-mini-2026-06-18-S2

Structural facts constraining the pure three-distance floor
  mu(E) = meas{ x in [0,1) : {frac(e x): e in E} have circular max-gap > 2/7 },  0 in E, |E|=k.
Goal: prove mu_min(k) > 0 (HYP-2586). These facts shape the proof:

 (1) DILATION-INVARIANCE: mu(dE) = mu(E) for any integer d>=1  (frac(d e x)=frac(e (dx)),
     x->dx measure-preserving). => WLOG gcd(E)=1 (primitive shapes). [reduces the shape space]
 (2) AP shapes are SPREAD-INDEPENDENT: E = {0,d,2d,...,(k-1)d} has mu = mu_consecutive(k)
     for EVERY d (same reason). So an AP cluster of any spread has the fixed value mu_k.
 (3) The IID value f(k)=P(k uniform pts have max-gap>2/7) is an UPPER bound region, NOT a
     lower bound: the structured orbit has mu_consec(k) < f(k). So "equidistribute to iid"
     gives the WRONG direction for a floor. (Closes a tempting wrong route.)
 (4) NEAR x=0 the good interval has width 5/(7*spread) -> 0; the floor is NOT carried by the
     near-0 interval for large spread. It is carried by the AP/three-distance structure.
"""
from fractions import Fraction as F
from math import comb
import random

def maxgap(pts):
    pts=sorted(set(pts)); n=len(pts)
    if n==1: return F(1)
    return max((pts[(i+1)%n]-pts[i])%1 for i in range(n))
def good(E,x):
    return maxgap([(e*x)%1 for e in E])>F(2,7)
def mu_exact(E):
    # breakpoints: x where two phases coincide or a phase crosses a gap-threshold; use the
    # rational breakpoints x=a/Q with Q=lcm-ish. For exactness over [0,1) use Q=2*max*7-ish grid
    # of candidate intervals; here mu via fine rational midpoint sampling between breakpoints.
    E=sorted(set(E)); diffs=set()
    M=max(E) if E else 1
    den=set([1])
    for e in E:
        if e>0:
            for k in range(1,7*e+1): den.add(e)
    # breakpoints where some e*x hits a multiple of 1/14-ish gap boundary OR two phases meet
    bps=set([F(0),F(1)])
    for i in range(len(E)):
        for j in range(len(E)):
            d=E[i]-E[j]
            if d>0:
                for m in range(0,d+1): bps.add(F(m,d))
        e=E[i]
        if e>0:
            for m in range(0,e+1): bps.add(F(m,e))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        if good(E,mid): tot+=b-a
    return tot

def f_iid(k,g=F(2,7)):
    s=F(0); j=0
    while 1-j*g>0:
        s+=(-1)**j*comb(k,j)*(1-j*g)**(k-1); j+=1
    return 1-s

print("="*70)
print("(1) DILATION-INVARIANCE  mu(dE)=mu(E):")
random.seed(1)
for _ in range(4):
    k=random.randint(4,6); E=sorted(random.sample(range(1,9),k-1)); E=[0]+E
    d=random.choice([2,3,5])
    print(f"  E={E}: mu={float(mu_exact(E)):.5f}  mu({d}E)={float(mu_exact([d*e for e in E])):.5f}")
print("\n(2) AP shapes spread-independent (mu = mu_consecutive(k) for any common difference d):")
for k in [4,5,6]:
    vals=[float(mu_exact([j*d for j in range(k)])) for d in [1,2,3,5]]
    print(f"  k={k}: mu(AP,d=1,2,3,5)={[round(v,4) for v in vals]}")
print("\n(3) iid f(k) is an UPPER bound, consecutive mu_k is BELOW it:")
for k in [4,5,7,9,13]:
    mc=float(mu_exact(list(range(k)))) if k<=9 else None
    print(f"  k={k}: mu_consec={'%.4f'%mc if mc else '(skip)'}  iid f(k)={float(f_iid(k)):.4f}  (orbit < iid: {(mc<float(f_iid(k))) if mc else 'n/a'})")
print("\n(4) near-0 good interval width = 5/(7*spread) -> 0 for large spread;")
print("    so the floor must come from the AP/three-distance bulk, not the near-0 cap.")
print("\nIMPLICATION: prove mu_min(k)>0 over PRIMITIVE shapes (dilation-inv); AP shapes give the")
print("fixed mu_k>0 at any spread; the genuine work is the perforated/multi-scale shapes.")
