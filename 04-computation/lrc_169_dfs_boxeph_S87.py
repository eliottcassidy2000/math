#!/usr/bin/env python3
"""
DFS around 169 = 13^2: the M<1/13 => rho>=13 crux  (boxeph-2026-07-18-S87)
==========================================================================
Goal: M<1/13 => rho>=13 closes LRC(14) (with THM-1008: rho>=13 => M>=1/14).
Equivalent: rho<13 covering => M>=1/13.  The M<1/13 families are deep-well-like
and cover 13 via a FAR multiple.  169=13^2 is the boundary; explore why.

BRANCHES:
 A) the 13-multiple threshold: for which k does '13k in V' force M>=1/13?
 B) structure of ALL M<1/13 families (rho, killer, core, continued fraction)
 C) continued-fraction / Ostrowski: is a1>=13 for M<1/13?
 D) the witness for '13k, k<=13 => M>=1/13'
 E) does M<1/13 => rho>=13 hold, and what is the tightest rho?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def exact_M_full(V):
    """M and its maximizer (a,q) via pair-sum denominators (THM-999)."""
    best=F(0); bt=None; qs=set([13,14])
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for dd in range(1,s+1):
                if s%dd==0: qs.add(dd)
    for q in qs:
        for a in range(1,q):
            if gcd(a,q)==1:
                m=min(min((v*a)%q,q-(v*a)%q) for v in V); c=F(m,q)
                if c>best: best=c; bt=(a,q,m)
    return best,bt
def cov(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1
def rho(V):
    Vs=sorted(V); return F(Vs[-1],Vs[-2])
def cont_frac(fr, depth=6):
    out=[]; x=fr
    for _ in range(depth):
        a=x.numerator//x.denominator; out.append(a); x=x-a
        if x==0: break
        x=1/x
    return out

print("="*88)
print("BRANCH A: for which k does '13k in V' force M >= 1/13 ?  (169 = 13*13 boundary)")
print("="*88)
# families {1..12, 13k}: when is M < 1/13?  (the {1..12} core is the tightest)
for k in range(1, 20):
    V=sorted(list(range(1,13))+[13*k])
    if len(set(V))!=13:
        print(f"  k={k:2d}: 13k={13*k} already in 1..12"); continue
    M,bt=exact_M_full(V)
    print(f"  k={k:2d}: 13k={13*k:3d}  M={str(M):9s}({float(M):.5f})  <1/13:{M<F(1,13)}  "
          f"cover={cov(V)} rho={float(rho(V)):.1f}  maximizer a/q={bt[0]}/{bt[1]} val={bt[2]}")

print("\n" + "="*88)
print("BRANCH B+E: hunt ALL M<1/13 families; record rho, killer, core, continued fraction")
print("="*88)
random.seed(50)
below=[]
def consider(V):
    V=sorted(set(V))
    if len(V)!=13 or not prim(V) or not cov(V): return
    M,bt=exact_M_full(V)
    if M<F(1,13):
        below.append((M,V,float(rho(V)),cont_frac(M),bt))
# structured deep-well-like + random
for core_size in [11,12]:
    core=list(range(1,core_size+1))
    for w in range(150, 400):
        consider(core+[w] if core_size==12 else core+[13,w])
for k in range(1,7):
    consider(list(range(1,13))+[182*k])
for _ in range(4000):
    consider(random.sample(range(1,220),13))
# adversarial: two large elements
for w1 in [169,182,195,364]:
    for w2 in range(150,400,7):
        consider(list(range(1,12))+[w1,w2])
below.sort(key=lambda t:t[0])
print(f"found {len(below)} M<1/13 families")
minrho=min((r for _,_,r,_,_ in below), default=None)
print(f"MIN rho among M<1/13 families: {minrho}  (claim: all >= 13)")
print("  smallest-M families (M, rho, cont-frac(M), maximizer q):")
for M,V,r,cf,bt in below[:16]:
    print(f"    M={str(M):9s} rho={r:6.1f} cf={cf} q={bt[1]}  V={V}")

print("\n" + "="*88)
print("BRANCH C: continued fraction of M for M<1/13 families -- is a1>=13 always?")
print("="*88)
a1s = sorted(set(cf[1] if len(cf)>1 else None for _,_,_,cf,_ in below))
print(f"  distinct a1 (first partial quotient of M) over M<1/13 families: {a1s}")
print(f"  all a1 >= 13 ? {all((cf[1]>=13) for _,_,_,cf,_ in below if len(cf)>1)}")

print("\n" + "="*88)
print("BRANCH D: the WITNESS for '13k in V, k<=13 => M>=1/13' -- what t certifies?")
print("="*88)
for k in [1,2,3,7,13]:
    V=sorted(list(range(1,13))+[13*k])
    if len(set(V))!=13: continue
    M,bt=exact_M_full(V)
    print(f"  {{1..12,{13*k}}}: M={M} witness t={bt[0]}/{bt[1]}  (does dilated-sieve d exist?)")
