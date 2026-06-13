#!/usr/bin/env python3
"""
The consecutive-block extremal family for the observer-inclusive signed mutual
LRC.  monad-explorer-2026-06-06-S3.

Empirically inf_S Gstar_full(S) = 2/(4n-5), achieved by the mid-range block
   B_n = {n-1, n, n+1, ..., 2n-3}   (n-1 consecutive speeds, max = 2n-3).
This script:
  - confirms Gstar_full(B_n) = 2/(4n-5) for n up to N_MAX,
  - extracts the REALIZING cut (A,B) and the witness time t*,
  - prints the binding speeds (those w with ||w t*|| = M) to expose the mechanism.
"""
from fractions import Fraction as F
from itertools import combinations, product

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1-f)

def candidate_times(Ws):
    T=set()
    for w in Ws:
        for a in range(0,w): T.add(F(2*a+1,2*w))
    for x,y in combinations(Ws,2):
        for d in (x+y,x-y):
            if d==0: continue
            d=abs(d)
            for a in range(1,d): T.add(F(a,d))
    return T

def maximin_arg(Ws):
    best=F(0); bt=None
    for t in candidate_times(Ws):
        ok=True; m=F(1,2)
        for w in Ws:
            nv=norm(w*t)
            if nv<=best: ok=False; break
            if nv<m: m=nv
        if ok and m>best: best=m; bt=t
    return best,bt

def cutWs(V,side):
    r=len(V); W=list(V)
    for i,j in combinations(range(r),2):
        W.append(abs(V[i]-V[j]) if side[i]==side[j] else V[i]+V[j])
    return sorted(set(W))

def gstar_full_arg(V):
    r=len(V); best=F(-1); bs=None; bt=None
    for tail in product([0,1],repeat=r-1):
        side=(0,)+tail
        Ws=cutWs(V,side)
        g,t=maximin_arg(Ws)
        if g>best: best=g; bs=side; bt=t
    return best,bs,bt

def main():
    print("Consecutive-block family B_n={n-1,..,2n-3}: Gstar_full vs 2/(4n-5)")
    print("="*70)
    for n in range(3,11):
        V=tuple(range(n-1,2*n-2))   # n-1, ..., 2n-3
        pred=F(2,4*n-5)
        g,side,t=gstar_full_arg(V)
        A=[V[i] for i in range(len(V)) if side[i]==0]
        Bs=[V[i] for i in range(len(V)) if side[i]==1]
        Ws=cutWs(V,side)
        binding=sorted(w for w in Ws if norm(w*t)==g)
        ok = (g==pred)
        print(f" n={n:2d} B_n={V}")
        print(f"     Gstar_full={g}={float(g):.6f}  2/(4n-5)={pred}  {'MATCH' if ok else 'DIFF!'}")
        print(f"     cut A={A} | B={Bs}   t*={t}={float(t):.5f}  (4n-5={4*n-5})")
        print(f"     binding speeds (||w t*||=M): {binding}")

if __name__=="__main__":
    main()
