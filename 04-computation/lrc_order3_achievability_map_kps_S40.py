#!/usr/bin/env python3
"""
kps-S40: the ORDER-3 achievability map across N -- the deeper analog of opus's mod-30 mediant
gate (HYP-4516: mediant achievable <=> N==1 mod 6 AND 5-nmid-(3N+2), a binder congruence mod
primes {2,3,5}).  mac-mini THM-632 proved N=12 misses the MEDIANT by parity; the OPEN frontier
is the deeper orders (kps S39: order-3 members are dilated APs).

Order-3 in-gap values at N: 4/(4N+3) and 5/(5N+3).  For each N we search dilated-AP + boundary
defect families (the S39 structure) for these values and record achievability, then look for the
N-periodic (mod 2,3,5,6,10,15,30) pattern -- the order-3 analog of the mediant's N==1 mod 6.

Integrates: S39 dilated-AP structure; opus mod-30 binder gate; the residue/witness picture.
"""
from fractions import Fraction
from itertools import combinations
import numpy as np
from math import gcd, isqrt
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd)
def isprime(n): return n>=2 and all(n%p for p in range(2,isqrt(n)+1))
def fact(n):
    f=[];d=2
    while d*d<=n:
        while n%d==0: f.append(d); n//=d
        d+=1
    if n>1: f.append(n)
    return f

def achievable(N, target):
    """search dilated-AP + boundary-defect families of N speeds for M==target (order-3 struct)."""
    q=target.denominator
    Hmax=min(2*q+4, 90)
    for d in range(2,8):                       # AP spacing
        for L in range(max(3,N-3), N):         # AP length (few defects)
            need=N-L
            for a in range(1, d+2):
                ap=[a+i*d for i in range(L)]
                if ap[-1]>Hmax: continue
                pool=[x for x in range(1,Hmax+1) if x not in ap]
                # bias defects near AP endpoints (the S39 structure)
                near=[x for x in pool if any(abs(x-e)<=d for e in (ap[0],ap[-1]))]
                cand = near if len(near)>=need else pool
                combos=list(combinations(cand, need))
                if len(combos)>1500: combos=combos[::max(1,len(combos)//1500)]
                for defs in combos:
                    v=sorted(set(ap)|set(defs))
                    if len(v)!=N or reduce(gcd,v)!=1: continue
                    if Mw(v)==target:
                        return v
    return None

print("=== ORDER-3 achievability map (dilated-AP search), N=5..18 ===\n", flush=True)
print(f"  {'N':>3}{'N%6':>5}{'N%30':>6}{'4/(4N+3)':>11}{'q=4N+3 fact':>16}{'ach?':>6}   {'5/(5N+3)':>11}{'q=5N+3 fact':>16}{'ach?':>6}", flush=True)
res={}
for N in range(5,19):
    t4=Fraction(4,4*N+3); t5=Fraction(5,5*N+3)
    w4=achievable(N,t4); w5=achievable(N,t5)
    res[N]=(w4 is not None, w5 is not None)
    print(f"  {N:>3}{N%6:>5}{N%30:>6}{str(t4):>11}{str(fact(4*N+3)):>16}{('YES' if w4 else 'no'):>6}   {str(t5):>11}{str(fact(5*N+3)):>16}{('YES' if w5 else 'no'):>6}", flush=True)

print("\n=== pattern analysis: which N give an order-3 member? ===", flush=True)
ach=[N for N in res if res[N][0] or res[N][1]]
print(f"  N with an order-3 member (4/(4N+3) OR 5/(5N+3)): {ach}", flush=True)
for m in (2,3,5,6,10,15,30):
    classes=sorted(set(N%m for N in ach))
    print(f"    N mod {m}: {classes}", flush=True)
print("\n  contrast: mediant (order-2) is N==1 mod 6 AND 5-nmid-(3N+2) (opus mod-30 gate).", flush=True)
print("  Look for the analogous congruence for order-3.  N=12 should be 'no' (matches (G)).", flush=True)
