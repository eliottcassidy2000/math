#!/usr/bin/env python3
"""
kps-S39: the ORDER-3 non-monotonicity -- why is N=6 first gap NONEMPTY (order-3 value 5/33,
kps S34) but N=12's order-3 (4/51,5/63) EMPTY (kps S38)?  mac-mini closed the MEDIANT (order-2):
F(N) gap member iff N==1 mod 6 (HYP-4572).  This goes one order deeper.

Order-3 in-gap values at N: s/(Ns+3), 3<s<6 => 4/(4N+3), 5/(5N+3).
Part 1: find the N=6 witness for 5/33 by residue search (speeds mod 33 in {5..28} => M>=5/33 at
        t=1/33; then check exact M) -- minimal-height representatives.
Part 2: NON-MONOTONICITY MAP -- for N=5..16, structured family (AP{1..N-2}+defect+outlier) gap hit
        and its order (mediant k=2 vs deeper k=3).
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
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0)
        qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd),(bc,bd)
def isprime(n): return n>=2 and all(n%p for p in range(2,isqrt(n)+1))
def order_of(M,N):
    p,q=M.numerator,M.denominator; return p, q-N*p

print("=== PART 1: N=6 witness for 5/33 (residue search, speeds mod 33 in {5..28}) ===", flush=True)
N=6; target=Fraction(5,33); band=list(range(5,29))   # dist>=5 mod 33
found=[]
for combo in combinations(band,6):
    v=list(combo)
    if reduce(gcd,v)!=1: continue
    M,(c,q)=Mw(v)
    if M==target:
        res=[(x*c)%q for x in v]; dist=[min(r,q-r) for r in res]
        found.append((v,c,q,dist))
        if len(found)<=6:
            print(f"  {v}: M=5/33 at t={c}/{q}; residues {res} dists {dist}", flush=True)
print(f"  total 6-subsets of {{5..28}} with M=5/33: {len(found)}", flush=True)
if found:
    v=found[0][0]
    print(f"  EXAMPLE structure {v}: diffs {[v[i+1]-v[i] for i in range(5)]}; is it AP+defects?", flush=True)

print("\n=== PART 2: NON-MONOTONICITY MAP -- structured gap members by N and order ===", flush=True)
print(f"  {'N':>3}{'N%6':>5}{'gap':>15}  {'best structured gap member':>42}{'k':>4}", flush=True)
rows=[]
for N in range(5,17):
    lo,hi=Fraction(1,N+1),Fraction(2,2*N+1)
    best=None
    for swap in range(N-2,N+5):
        b=sorted(set(range(1,N-1))|{swap})
        if len(b)!=N-1: continue
        mu,(cb,D)=Mw(b)
        if mu<=hi: continue
        for j in range(1,8):
            for x in (j*D-1,j*D,j*D+1):
                if x<=max(b): continue
                w=sorted(b+[x])
                if len(set(w))!=N or reduce(gcd,w)!=1: continue
                M,_=Mw(w)
                if lo<M<hi and (best is None or M<best[1]):
                    best=(w,M)
    if best:
        s,k=order_of(best[1],N); rows.append((N,N%6,k))
        print(f"  {N:>3}{N%6:>5}{f'(1/{N+1},2/{2*N+1})':>15}  {str(best[0])+'='+str(best[1]):>42}{k:>4}", flush=True)
    else:
        rows.append((N,N%6,None))
        print(f"  {N:>3}{N%6:>5}{f'(1/{N+1},2/{2*N+1})':>15}  {'(no structured gap member)':>42}", flush=True)

print("\n  by N%6:", flush=True)
for r in range(6):
    ns=[str(x[0]) for x in rows if x[1]==r]
    ks=[str(x[2]) for x in rows if x[1]==r]
    print(f"    N=={r} mod 6: N={ns}  orders k={ks}", flush=True)
print("\nREADING: mediant (k=2) at N==1 mod 6 (mac-mini). Deeper orders fill other residues where", flush=True)
print("they can. N=12 (==0 mod 6) shows no structured gap member here => consistent with EMPTY;", flush=True)
print("N=6 (==0 mod 6) DID (5/33 order-3) -- so N%6 alone does NOT decide order-3 (the real crux).", flush=True)
