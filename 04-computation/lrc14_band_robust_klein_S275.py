#!/usr/bin/env python3
"""
lrc14_band_robust_klein_S275.py
===============================
klein-2026-07-13-S275. Firm up test (C): max Phi over the moderate-w band for MANY primitive
7-clusters (structured worst-cases + random), w in [26, 6*diam]. If max Phi stays well < cap9,
the k=8 row is robustly TRUE on the band (the residual is rigor, not truth).
Prime grid (moderate Ng; band Phi-values don't need resonance resolution).
"""
import math
from math import gcd
from functools import reduce
NG=200003
def gcdall(xs): return reduce(gcd,[x for x in xs if x],0)
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def moments3(E):
    s1=s2=s3=0
    for k in range(1,NG):
        N=7-bin(occ_of(E,k/NG)).count("1"); s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;return s1/n,s2/n,s3/n
def Phi(E):m1,m2,m3=moments3(E);return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
CAP9=1979/4004

# deterministic pseudo-random primitive 7-clusters (no Math.random; use a fixed LCG on index)
def lcg(seed):
    s=seed
    while True:
        s=(1103515245*s+12345)&0x7fffffff; yield s
def rand_prim7(idx, dmax):
    g=lcg(idx*7919+1); elts={0}
    while len(elts)<7:
        r=next(g); elts.add(1+r%dmax)
    C=sorted(elts)
    if gcdall(C)!=1: C[1]=1  # force primitivity cheaply
    return sorted(set(C)) if len(set(C))==7 else rand_prim7(idx+1,dmax)

structured=[
  [0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,1,2,4,6,8,10],[0,2,4,6,8,10,12],
  [0,1,2,3,4,8,16],[0,1,3,6,10,15,21],[0,1,2,3,5,8,13],[0,5,10,15,20,25,31],
]
supPhi=0; supArg=None; cnt=0
def scan(C):
    global supPhi,supArg,cnt
    if gcdall(C)!=1 or len(set(C))!=7: return
    d=max(C)
    for w in [26,30,37,d+1, (3*d)//2, 2*d, 3*d, 6*d]:
        if w<26: continue
        E=sorted(set(C+[w]))
        if len(E)!=8 or gcdall(E)!=1: continue
        p=Phi(E); cnt+=1
        if p>supPhi: supPhi=p; supArg=(C,w)
for C in structured: scan(C)
for i in range(40):
    for dmax in (12,25,45):
        scan(rand_prim7(i,dmax))
print(f"scanned {cnt} primitive 8-cores on the band 26<=w<=6*diam")
print(f"MAX Phi on band = {supPhi:.5f} at {supArg}")
print(f"cap9 = {CAP9:.5f}, margin +{CAP9-supPhi:.5f}")
print(f"=> k=8 row {'ROBUSTLY TRUE on the band (margin>0.1)' if CAP9-supPhi>0.1 else 'check'}")
print("\ndone.")
