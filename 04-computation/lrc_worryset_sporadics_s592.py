#!/usr/bin/env python3
"""INVESTIGATE the worry-set sporadics: exhaustively find tight configs (M=1/n) in a box,
classify AP vs sporadic, and characterize: transversal mod 2n-1? round? chi (THM-402=>2?).
opus-2026-06-03-S592 round3."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def safe_measure_pos(V,n):
    THR=F(1,n); eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(1,-1): eps.add(F(k*n+s,n*v)%1)
    pts=sorted(eps); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>THR for v in V): return True
    return False
def floor_tight(V,n):  # measure-0 but lonely at some j/n
    return (not safe_measure_pos(V,n)) and any(all(dist(v*F(j,n))>=F(1,n) for v in V) for j in range(1,n))
def is_transversal(V,n):
    C=2*n-1; shells=set()
    for v in V:
        r=v%C
        if r==0: return False
        shells.add(min(r,C-r))
    return len(shells)==len(V) and len(V)==n-1
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    print("Worry-set (floor-tight M=1/n) configs in box; AP vs sporadic; transversal mod 2n-1?")
    for n in [6,8,10]:
        m=n-1; B={6:14,8:16,10:18}[n]
        tights=[]
        for combo in combinations(range(1,B+1),m):
            V=prim(combo)
            if V!=tuple(sorted(combo)): continue
            if floor_tight(V,n): tights.append(V)
        AP=tuple(range(1,n))
        spor=[V for V in tights if V!=AP]
        tvs=sum(1 for V in tights if is_transversal(V,n))
        print(f"  n={n}: floor-tight={len(tights)} (AP present={AP in tights}); sporadics={len(spor)}: {spor[:5]}")
        print(f"       transversals mod {2*n-1}: {tvs}/{len(tights)}; non-transversal sporadics: {[V for V in spor if not is_transversal(V,n)][:4]}")
if __name__=='__main__': main()
