#!/usr/bin/env python3
"""
Verify the HYP-2293 counterexample V=(2,3,4,6,8) at n=6 (r=5 movers):
Gstar(V) = 3/19 < 1/6.  Three independent methods + the optimal cut.
monad-explorer-2026-06-06-S2d.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce

def norm(x):
    f=x-int(x);  f=f+1 if f<0 else f;  return min(f,1-f)
def cand(W):
    W=[w for w in W if w]; ms=set(abs(w) for w in W)
    for a,b in combinations(W,2): ms.add(abs(a+b)); ms.add(abs(a-b))
    ms.discard(0); T=set()
    for m in ms:
        for a in range(1,2*m): T.add(F(a,2*m))
    return T
def gap_exact(W):
    W=[w for w in W if w]; 
    if not W: return F(1,2)
    return max(min(norm(w*t) for w in W) for t in cand(W))
def gap_float(W,N):
    W=[abs(w) for w in W if w]; best=0.0
    for a in range(1,N):
        t=a/N; m=1.0
        for w in W:
            f=(w*t)%1.0; d=f if f<0.5 else 1-f
            if d<m: m=d
            if m<=best: break
        if m>best: best=m
    return best
def sdiff(V,eps): return [eps[i]*V[i]-eps[j]*V[j] for i,j in combinations(range(len(V)),2)]

V=(2,3,4,6,8); n=len(V)+1
print(f"V={V}  gcd={reduce(gcd,V)}  n={n}  floor 1/n={F(1,n)}={1/n:.5f}")
best=F(0); bestcut=None
for tail in product([1,-1],repeat=len(V)-1):
    eps=(1,)+tail
    g=gap_exact(sdiff(V,eps))
    if g>best: best=g; bestcut=eps
print(f"\nGstar (exact, max over all cuts) = {best} = {float(best):.5f}")
print(f"  optimal sign pattern eps = {bestcut}")
A=[V[i] for i in range(n-1) if bestcut[i]==1]; B=[V[i] for i in range(n-1) if bestcut[i]==-1]
print(f"  cut: A(+)={A}  B(-)={B}")
print(f"  relative-speed multiset D = {sorted(abs(d) for d in sdiff(V,bestcut))}")
print(f"  gcd of |D| = {reduce(gcd,[abs(d) for d in sdiff(V,bestcut) if d])}")
# independent float check of the OPTIMAL cut at high resolution
gf=gap_float(sdiff(V,bestcut), 5_000_000)
print(f"\nindependent float grid (N=5e6) on optimal cut: {gf:.6f}  (exact {float(best):.6f})")
print(f"3/19 = {float(F(3,19)):.6f}")
print(f"\nVERDICT: Gstar={best} {'<' if best<F(1,n) else '>='} 1/n={F(1,n)}  "
      f"=> HYP-2293 {'REFUTED' if best<F(1,n) else 'holds here'}")
# why: every cut's D shares a common factor? check gcd over ALL cuts
print("\nper-cut |D| gcd (shows imprimitivity of relative-speed set):")
gcds={}
for tail in product([1,-1],repeat=len(V)-1):
    eps=(1,)+tail; d=[abs(x) for x in sdiff(V,eps) if x]
    g=reduce(gcd,d); gcds[g]=gcds.get(g,0)+1
print(f"  distribution of gcd(|D|) over the {2**(len(V)-1)} cuts: {dict(sorted(gcds.items()))}")
