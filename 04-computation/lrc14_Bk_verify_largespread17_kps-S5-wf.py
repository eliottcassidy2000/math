#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_largespread17_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)
Direct test of the prompt's structural claim for the 1/7 measure:
  'for theta=1/7, ANY larger spread RAISES mu_1/7 above consecutive' (opposite of 2/7).
We compute, per k=8..13, the per-spread MINIMUM of mu_1/7 (exhaustive where feasible, else
the densest-ruler family: perforated APs of that spread), and confirm it is >= consecutive,
and trends UP (not U-shaped down) as spread grows. A dip below consec at any spread refutes
the 1/7 spread bound.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(2024)
ONE7=F(1,7)

def mu17(E):
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
                if be>ONE7: ivs.append((a,b))
            else:
                xs=(ONE7-be)/al
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

print("Per-spread MIN mu_1/7 (exhaustive if feasible else perforated-AP family), trend check.")
print("If the 1/7 spread bound holds, min over spread s should be the consecutive value at s=k-1,")
print("and NON-DECREASING (no dip) as s grows.\n")
violated=False
for k in range(8,14):
    consec=mu17(list(range(k)))
    row=[]
    for s in [k-1, k, k+1, k+2, k+4, k+8, 2*k, 3*k, 5*k, 10*k]:
        if s<k-1: continue
        best=F(2)
        if comb(s-1,k-2)<=20000:
            for interior in itertools.combinations(range(1,s),k-2):
                E=(0,)+interior+(s,)
                if not gcd1(E): continue
                m=mu17(list(E))
                if m<best: best=m
        else:
            # densest-ruler / perforated-AP family of this spread + random sample
            cands=set()
            full=list(range(s+1))
            need=(s+1)-k
            if need>=0 and comb(s-1,need)<=20000:
                for drops in itertools.combinations(range(1,s),need):
                    E=tuple(x for x in full if x not in drops)
                    if len(E)==k and E[0]==0 and E[-1]==s: cands.add(E)
            for _ in range(3000):
                body=sorted(random.sample(range(1,s),k-2))
                E=tuple([0]+body+[s])
                if len(set(E))==k: cands.add(E)
            for E in cands:
                m=mu17(list(E))
                if m<best: best=m
        row.append((s,best))
        if best<consec: violated=True
    desc=" ".join(f"s{s}:{float(b):.3f}" for s,b in row)
    flag="" if all(b>=consec for _,b in row) else "  *** DIP BELOW CONSEC ***"
    print(f"  k={k:2d} consec={float(consec):.3f}  {desc}{flag}")
print(f"\n1/7 spread bound (min over spread == consecutive, no dip): {'HOLDS' if not violated else 'VIOLATED'}")
