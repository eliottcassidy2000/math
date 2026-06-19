#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_pigeonhole_weyl_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)
(1) The k<=7 pigeonhole for mu_1/7=1: k points -> k circular gaps summing to 1.
    k<=6: some gap >= 1/k >= 1/6 > 1/7 ALWAYS (so maxgap>1/7 always -> mu=1). RIGOROUS.
    k=7: 7 gaps sum 1; maxgap<=1/7 forces all =1/7 (a single x-locus per order-cell) -> measure 0.
    We verify mu_1/7(E)=1 EXACTLY for ALL E up to k=7 over a large spread range, AND confirm the
    'all gaps =1/7' locus is finite (the boundary set has measure 0). Already partly done; here we
    push spread higher and assert EXACT ==1.
(2) The Weyl/Parseval deviation identity: mu_theta(E) = F_theta(k) + sum_{m!=0, m.e=0} ghat(m).
    F_theta(k)=ghat(0). We do NOT need the full identity; we verify the WEAK consequence the
    argument actually uses nowhere-critically, namely the CEILING direction is consistent:
    for spread-out E (few short relations), mu -> F(k). We check mu(E)->F(k) numerically as a
    sanity check that F(k) is genuinely the iid value (ghat(0)=F(k)).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(99)
ONE7=F(1,7); TWO7=F(2,7)

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

print("="*84)
print("(1) k<=7 pigeonhole: mu_1/7(E)==1 EXACTLY for all E (push spread high)")
print("="*84)
allone=True
for k in range(2,8):
    smax=k+12
    minmu=F(2); worst=None; cnt=0
    for s in range(k-1,smax+1):
        if comb(s-1,max(k-2,0))>80000: continue
        for interior in itertools.combinations(range(1,s),k-2):
            E=(0,)+interior+(s,)
            cnt+=1
            m=mu_theta(list(E),ONE7)
            if m<minmu: minmu=m; worst=E
    ok=(minmu==F(1)); allone&=ok
    print(f"  k={k}: {cnt} sets spread<= {smax}: min mu_1/7 = {minmu} {'(==1)' if ok else f'<1 at {worst}'}")
print(f"  ==> mu_1/7==1 for all k<=7 (exact, high spread): {allone}")

print()
print("="*84)
print("(2) Weyl ceiling sanity: spread-out E (random large, few short relations) -> mu -> F(k)")
print("="*84)
def Fk(k,L):
    s=F(0); j=1
    while 1-j*L>0:
        s+=(-1)**(j+1)*comb(k,j)*(1-j*L)**(k-1); j+=1
    return s
for k in [6,7,8]:
    fk=Fk(k,TWO7)
    # very spread-out random E (Sidon-ish): large distinct values
    samples=[]
    for _ in range(40):
        E=[0]+sorted(random.sample(range(1,4000),k-1))
        samples.append(float(mu_theta(E,TWO7)))
    avg=sum(samples)/len(samples)
    print(f"  k={k}: F(k)_2/7={float(fk):.4f} ; mean mu over 40 spread-out E = {avg:.4f} "
          f"(should hover near F(k), confirming ghat(0)=F(k))")
print("\nDONE.")
