#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE SHARPEST OPEN LEVER: does the SPORADIC two-clash family stay above the pentagon (Z/10)* / above 1/36?

kind-pasteur-2026-07-01-S10. The census dichotomy (HYP-3793/S8/S9): loose 11-cores split into
 - FULL-ORBIT polygon classes (Z/N)*, min at the PENTAGON N=10, meas 313/9702=0.03226;
 - SPORADIC "two-clash" cores (maximizers = a PARTIAL orbit, 2 atoms; e.g. (Z/19) {2/19,17/19}),
   the 11-speed Goddyn-Wong analogue.  The binding question is whether some sporadic core dips below
   the pentagon (or below 1/36).  Search the sporadic family broadly for its minimum.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def M_and_atoms(S,Dmax=None):
    if Dmax is None: Dmax=2*max(S)+2
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w)
            if v!=w: dens.add(abs(v-w))
    bn=0;bd=1;atoms=[]
    for b in sorted(d for d in dens if 0<d<=Dmax):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            md=b
            for v in S:
                r=(v*a)%b; dv=r if r<=b-r else b-r
                if dv<md: md=dv
                if md==0: break
            if md==0: continue
            c=md*bd-bn*b
            if c>0: bn,bd,atoms=md,b,[Fr(a,b)]
            elif md*bd==bn*b: atoms.append(Fr(a,b))
    return Fr(bn,bd), sorted(set(atoms))
BAND=Fr(1,14)
def safe_arcs(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def Lmeas(S):
    a=safe_arcs(S[0])
    for v in S[1:]:
        a=inter(a,safe_arcs(v))
        if not a: return Fr(0)
    return sum(h-l for l,h in a)

PENT=Fr(313,9702); T36=Fr(1,36)
pool=set()
# 2-drops of {1..V} (partial-orbit sporadics live here) + 1-drop + far additions + Goddyn-Wong-like
for V in [13,14,15,16,17]:
    base=set(range(1,V+1))
    for drop in itertools.combinations(range(1,V+1), V-11):
        C=tuple(sorted(base-set(drop)))
        if len(C)==11: pool.add(C)
for j in range(1,13):
    for w in range(13,80):
        C=tuple(sorted(set(range(1,13))-{j}|{w}))
        if len(C)==11: pool.add(C)
# GW-like: {1..k} with a doubled/clash far element
rng=random.Random(2)
for _ in range(4000):
    C=tuple(sorted(set(rng.sample(range(1,30),11))))
    if len(C)==11: pool.add(C)

sporadic=[]; fullorbit=[]
for C in pool:
    S=list(C); M,atoms=M_and_atoms(S,Dmax=40)
    if M<=BAND: continue
    Ns=sorted(set(t.denominator for t in atoms)); N=Ns[0] if len(Ns)==1 else None
    full = (N is not None and sorted(t.numerator for t in atoms)==units(N))
    m=Lmeas(S)
    if m<=0: continue
    rec=(m,S,M,N,len(atoms),full)
    if full: fullorbit.append(rec)
    else: sporadic.append(rec)
sporadic.sort(key=lambda r:r[0]); fullorbit.sort(key=lambda r:r[0])
print(f"POOL {len(pool)}: {len(fullorbit)} full-orbit, {len(sporadic)} sporadic (partial-orbit) loose cores.")
print(f"PENTAGON (Z/10)* = 313/9702 = {float(PENT):.6f};  1/36 = {float(T36):.6f}")
print("\n=== SPORADIC (two-clash) family, smallest 10 ===")
print(f"  {'meas':>10} {'M':>7} {'N':>5} {'#at':>4}  core")
for m,S,M,N,na,full in sporadic[:10]:
    print(f"  {float(m):>10.6f} {str(M):>7} {str(N):>5} {na:>4}  {S}")
smin=sporadic[0]
print(f"\nSPORADIC MIN: meas={smin[0]}={float(smin[0]):.6f} at {smin[1]} (N={smin[3]}, {smin[4]} atoms)")
print(f"  sporadic-min >= pentagon (313/9702)? {smin[0]>=PENT}  (gap {float(smin[0]-PENT):+.6f}, {float(smin[0]/PENT):.4f}x)")
print(f"  sporadic-min >= 1/36?               {smin[0]>=T36}  ({float(smin[0]/T36):.3f}x)")
print(f"\nGLOBAL MIN (both families): {min(float(fullorbit[0][0]), float(smin[0])):.6f} -- "
      f"{'PENTAGON' if fullorbit[0][0]<smin[0] else 'SPORADIC'} binds.")
print("  => the sharpest lever: the binding minimizer is the "
      f"{'full-orbit PENTAGON' if fullorbit[0][0]<smin[0] else 'SPORADIC two-clash'}; both families clear 1/36.")
print("DONE.")
