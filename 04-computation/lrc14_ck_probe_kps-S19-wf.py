#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Fast focused probe of the two load-bearing premises:
  (P1) K_1 saturation: single block {0..L-1} alone stays O(1) as L grows.
  (P2) per-cluster attributed contribution stays O(1) when embedded in a
       multi-cluster core (the S18 part (C) gave 4.787 -- much > K_1=1.108).
  (P3) hunt: does worst w|Delta_w| stay <= c*(k-2) with c~2.72?
Small, fast, exact.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

def G0(y):
    y = y - int(y)
    if y < 0: y += 1
    if y < F(1,7): return y*F(6,7)
    return F(6,49) - (y-F(1,7))*F(1,7)
def _f(y):
    y=y-int(y)
    if y<0: y+=1
    return y
def breakpoints(Ep):
    Ep=sorted(set(e for e in Ep if e!=0)); bp={F(0),F(1)}
    for e in Ep:
        for j in range(7):
            c=F(j,7); m=0
            while True:
                xv=(c+m)/e
                if xv>=1: break
                if xv>=0: bp.add(xv)
                m+=1
    return sorted(b for b in bp if 0<=b<1)
def cells(Ep):
    Ep=sorted(set(e for e in Ep if e!=0)); bp=breakpoints(Ep); out=[]
    for lo,hi in zip(bp,bp[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        hit=set(int(_f(e*mid)*7) for e in Ep)
        miss=set(range(1,7))-hit
        out.append((lo,hi,next(iter(miss)) if len(miss)==1 else None))
    return out
def wDsigned(Ep,w):
    D=F(0)
    for lo,hi,s in cells(Ep):
        if s is None: continue
        D+=G0(w*hi-F(s,7))-G0(w*lo-F(s,7))
    return D
def wD(Ep,w): return abs(float(wDsigned(Ep,w)))
def prim(E):
    nz=[e for e in E if e!=0]
    return len(nz)>0 and reduce(gcd,nz)==1
def lcm(a,b): return a*b//gcd(a,b)
def lcml(xs): return reduce(lcm,xs,1)
def clusters(Ep,ratio=3):
    es=sorted(e for e in Ep if e>0)
    if not es: return []
    cl=[[es[0]]]
    for a,b in zip(es,es[1:]):
        if b>a*ratio: cl.append([b])
        else: cl[-1].append(b)
    return cl
def nclu(Ep): return len(clusters(Ep))

print("="*70)
print("(P1) K_1 saturation: single block {0..L-1}, denser/larger w-scan")
print("="*70)
K1=0.0; wit=None
for L in (6,7,8,9,10,11,12):
    core=list(range(L))
    nz=[e for e in core if e>0]; Lall=lcml(nz)
    cand=set(range(2,1000))
    for m in range(1,40):
        if Lall*m<200000: cand.add(Lall*m)
    worst,ww=0.0,0
    for w in sorted(cand):
        v=wD(core,w)
        if v>worst: worst,ww=v,w
    if worst>K1: K1,wit=worst,(core,ww)
    print(f"  L={L:3d}  worst w|D|={worst:.4f}  w*={ww}  (lcm_all={Lall})")
print(f"  => K_1={K1:.4f} at {wit}")

print()
print("="*70)
print("(P3) HUNT: 3-cluster geometric resonant cores, w divisible by all scales")
print("="*70)
print(f"  {'core':<48}{'k':>3}{'r':>3}{'worst':>9}{'/(k-2)':>8}{'/(r-1)':>8}  w*")
worstg=0.0; gw=None
designs=[
    [0,1,2,30,31,32,900,901,902],
    [0,1,2,3,40,41,42,43,1600,1601],
    [0,1,2,20,21,22,400,401,402],
    [0,1,2,12,13,14,144,145,146],
    [0,1,2,10,11,12,100,101,102],
    [0,1,2,3,4,30,31,32,33,34],
    [0,1,2,15,16,17,225,226],
    [0,1,2,3,50,51,52,53,2500,2501],
]
for core in designs:
    core=sorted(set(core))
    if not prim(core):
        print(f"  {str(core):<48} not primitive"); continue
    nz=[e for e in core if e>0]
    cl=clusters(core); reps=[c[0] for c in cl if c[0]>0]
    Lrep=lcml(reps); Lall=lcml(nz)
    cand=set(range(2,300))
    for base in (Lrep,Lall):
        for m in range(1,25):
            if base*m<500000: cand.add(base*m)
    worst,ww=0.0,0
    for w in sorted(cand):
        v=wD(core,w)
        if v>worst: worst,ww=v,w
    k=len(core); r=nclu(core)
    if worst>worstg: worstg,gw=worst,(core,ww)
    print(f"  {str(core):<48}{k:>3}{r:>3}{worst:>9.4f}{worst/max(k-2,1):>8.3f}{worst/max(r-1,1):>8.3f}  {ww}")
print(f"  => worst={worstg:.4f} at {gw}")
print(f"     claim says C(k)<=2.72*(k-2); for the worst core k={len(gw[0])} bound={2.72*(len(gw[0])-2):.2f}")
