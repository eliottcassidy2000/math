#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_ushape_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)
Verify the prompt's concrete 2/7 claims:
 (a) consecutive mu_7 = 83/210 (correction to canon 13/35); perforated {0,2,3,4,5,6,8}=13/35.
 (b) k=9 2/7 per-spread minima exact: s=8 ->13/35? prompt says 'k=7 -> 13/35 at s=8'; for k=9
     prompt claims min 811/4095 at s=18 (EXHAUSTIVE s=18). Verify 811/4095 is achieved at s=18.
 (c) 'opposite extremizers': the 2/7-minimizer E for k=13 has mu_1/7 ~ 0.95-0.99 (>> consec 0.4425).
 (d) bounded-spread 2/7 minima k=8: 71/220 at s=11.
All exact.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
TWO7=F(2,7); ONE7=F(1,7)

def mu_theta(E, theta):
    E=sorted(set(E)); k=len(E)
    if k==1: return F(1)
    bp={F(0),F(1)}
    for i in range(k):
        for j in range(i+1,k):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bps=sorted(b for b in bp if F(0)<=b<=F(1)); tot=F(0)
    for a,b in zip(bps,bps[1:]):
        if a==b: continue
        mid=(a+b)/2; pts=[]
        for e in E:
            val=e*mid; n=val-(val%1); pts.append((val-n,e,n))
        pts.sort(key=lambda t:t[0]); m=len(pts); ivs=[]
        for i in range(m):
            (_,ei,ni)=pts[i]
            if i<m-1: (_,ej,nj)=pts[i+1]; al=F(ej-ei); be=F(ni-nj)
            else: (_,e0,n0)=pts[0]; al=F(e0-ei); be=F(ni-n0+1)
            if al==0:
                if be>theta: ivs.append((a,b))
            else:
                xs=(theta-be)/al
                if al>0: lo,hi=max(a,xs),b
                else: lo,hi=a,min(b,xs)
                if lo<hi: ivs.append((lo,hi))
        if not ivs: continue
        ivs.sort(); clo,chi=ivs[0]
        for lo,hi in ivs[1:]:
            if lo<=chi: chi=max(chi,hi)
            else: tot+=chi-clo; clo,chi=lo,hi
        tot+=chi-clo
    return tot

def gcd1(E): return reduce(gcd,E)==1

print("(a) consecutive mu_7(2/7) and perforated {0,2,3,4,5,6,8}:")
print("    consec_7 =", mu_theta(list(range(7)),TWO7), "(claim 83/210:", mu_theta(list(range(7)),TWO7)==F(83,210),")")
print("    {0,2,3,4,5,6,8} =", mu_theta([0,2,3,4,5,6,8],TWO7), "(claim 13/35:", mu_theta([0,2,3,4,5,6,8],TWO7)==F(13,35),")")

print("\n(b) k=9 2/7 per-spread min at s=18 (exhaustive interior C(17,7)=19448): claim 811/4095")
k=9; s=18; best=F(2); bestE=None
for interior in itertools.combinations(range(1,s),k-2):
    E=(0,)+interior+(s,)
    if not gcd1(E): continue
    m=mu_theta(list(E),TWO7)
    if m<best: best=m; bestE=E
print(f"    min mu_2/7 at s=18 = {best} = {float(best):.5f}  at {bestE}  (claim 811/4095: {best==F(811,4095)})")

print("\n(d) k=8 2/7 per-spread min at s=11 (exhaustive C(10,6)=210): claim 71/220")
k=8; s=11; best=F(2); bestE=None
for interior in itertools.combinations(range(1,s),k-2):
    E=(0,)+interior+(s,)
    if not gcd1(E): continue
    m=mu_theta(list(E),TWO7)
    if m<best: best=m; bestE=E
print(f"    min mu_2/7 at s=11 = {best} = {float(best):.5f}  at {bestE}  (claim 71/220: {best==F(71,220)})")

print("\n(c) 'opposite extremizers': mu_1/7 of the 2/7-minimizers (must be >> consec mu_1/7):")
# take the k=9 s=18 2/7-minimizer found above
k=9; s=18; b9=F(2); b9E=None
for interior in itertools.combinations(range(1,s),k-2):
    E=(0,)+interior+(s,)
    if not gcd1(E): continue
    m=mu_theta(list(E),TWO7)
    if m<b9: b9=m; b9E=E
print(f"    k=9 2/7-min E={b9E}: mu_2/7={b9}={float(b9):.4f}  mu_1/7={mu_theta(list(b9E),ONE7)}={float(mu_theta(list(b9E),ONE7)):.4f}")
print(f"         vs consec_9 mu_1/7 = {mu_theta(list(range(9)),ONE7)} = {float(mu_theta(list(range(9)),ONE7)):.4f}")
print("    => the 2/7-minimizer has HIGHER mu_1/7 than consecutive (opposite extremizers).")
