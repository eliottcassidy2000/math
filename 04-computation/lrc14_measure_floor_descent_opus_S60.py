#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
AGGRESSIVE descent to MINIMIZE meas(lonely S) over primitive covering families (opus-S60).
Tests whether the uniform measure floor (inf meas > 0, ~0.005 from S60 scans) survives a hard search
allowing large magnitude. If the descent cannot beat ~0.005 even with speeds <= 80, the floor is
robustly positive. Descent in FLOAT; finalists re-verified EXACT.
"""
import sys, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
MAXV=80; BANDf=1.0/14.0
def saf(v): return [((k+BANDf)/v,(k+1-BANDf)/v) for k in range(v)]
def interf(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=A[i][0] if A[i][0]>B[j][0] else B[j][0]
        hi=A[i][1] if A[i][1]<B[j][1] else B[j][1]
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
_mc={}
def measf(S):
    key=tuple(sorted(set(S)))
    if key in _mc:return _mc[key]
    a=saf(key[0])
    for v in key[1:]:
        a=interf(a,saf(v))
        if not a:_mc[key]=0.0;return 0.0
    m=sum(h-l for l,h in a);_mc[key]=m;return m
def cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def prim(S): return reduce(gcd,S)==1
def sae(v):
    B=Fr(1,14);return [((Fr(k)+B)/v,(Fr(k+1)-B)/v) for k in range(v)]
def intere(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]);hi=min(A[i][1],B[j][1])
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
def mease(S):
    S=sorted(set(S));a=sae(S[0])
    for v in S[1:]:
        a=intere(a,sae(v))
        if not a:return Fr(0)
    return sum(h-l for l,h in a)
POOL=list(range(1,MAXV+1))
def val(S):
    S=sorted(set(S))
    if len(S)!=13 or not cov(S) or not prim(S): return None
    return measf(S)
def descend(S):
    S=sorted(set(S));cur=val(S)
    if cur is None:return None,None
    imp=True
    while imp:
        imp=False
        for i in range(13):
            for w in POOL:
                if w in S:continue
                T=S[:i]+S[i+1:]+[w]
                v=val(T)
                if v is not None and 0<v<cur-1e-9:
                    S=sorted(set(T));cur=v;imp=True;break
            if imp:break
    return cur,S
seeds=[sorted([2,4,5,6,8,10,14,16,18,20,22,24,26]),
       sorted([1,2,3,4,5,7,8,9,10,11,12,13,14])]
random.seed(3);tries=0
while len(seeds)<400 and tries<300000:
    tries+=1
    S=random.sample([x for x in POOL if x<=40],13)
    if any(s%14==0 for s in S) and cov(S) and prim(S): seeds.append(sorted(S))
best=(9.9,None);lm={}
for s in seeds:
    v,S=descend(s)
    if v is None:continue
    lm[tuple(S)]=v
    if v<best[0]:best=(v,S)
Sb=best[1]
print("="*96);print(" MINIMIZE meas(lonely S), primitive covering, speeds<=%d, %d seeds"%(MAXV,len(seeds)));print("="*96)
print(f"  MIN meas (float): {best[0]:.6f}   S={Sb}")
print(f"  EXACT re-verify:  {float(mease(Sb)):.6f}   v_max={max(Sb)}  primitive={prim(Sb)} covering={cov(Sb)}")
print("\n  lowest 12 distinct local minima (EXACT meas):")
for S,v in sorted(lm.items(),key=lambda kv:kv[1])[:12]:
    print(f"    meas={float(mease(list(S))):.6f}  vmax={max(S):>3}  S={list(S)}")
print(f"\n  => min meas found = {best[0]:.5f} (v_max={max(Sb)}). If this stays ~0.005 and does NOT")
print("     shrink with larger allowed magnitude, the uniform measure floor inf meas>0 is robust.")
print("DONE.")
