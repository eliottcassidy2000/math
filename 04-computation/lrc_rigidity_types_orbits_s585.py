#!/usr/bin/env python3
"""Rigidity as ORBIT DIMENSION + types by which group acts. Witness set: rigid=0-dim
(finite orbit, worry-set) vs flexible=1-dim (interval, loose). Decompose the rigid
witness orbit under (Z/n)* (cyclic) and <x2> (doubling). opus-2026-06-03-S585."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def units(n): return [u for u in range(1,n) if gcd(u,n)==1]
def witness_clock(V,n):  # lonely clock points j in 1..n-1
    return [j for j in range(1,n) if all(dist(v*F(j,n))>=F(1,n) for v in V)]
def doubling_orbits(S,n):  # partition S (subset of Z/n) into <x2> orbits
    S=set(S); orbits=[]; seen=set()
    for x in sorted(S):
        if x in seen: continue
        o=[]; y=x
        while y not in seen:
            seen.add(y); o.append(y); y=(2*y)%n
            if y not in S: break
        orbits.append(o)
    return orbits
def main():
    print("RIGIDITY = witness orbit dimension; sub-orbit structure (cyclic units vs doubling)")
    for n in [7,8,12,14,16]:
        AP=tuple(range(1,n)); W=witness_clock(AP,n); U=units(n)
        print(f"  n={n}: AP witness clock-points = {W}  (= units (Z/n)*? {W==U}); |orbit|=phi(n)={len(U)}")
        print(f"       under DOUBLING x2 mod n, the unit witnesses split into orbits: {doubling_orbits(W,n)}")
    print()
    print("Flexible (loose) example -> positive-measure witness (dim 1):")
    for V in [(1,4,5),(2,3,7,8)]:
        n=len(V)+1
        # safe interval check
        THR=F(1,n); eps=set([F(0)])
        for v in V:
            for k in range(v+1):
                for s in (1,-1): eps.add(F(k*n+s,n*v)%1)
        pts=sorted(eps); meas=F(0); L=len(pts)
        for i in range(L):
            a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
            if all(dist(v*mid)>=THR for v in V): meas+=ln
        print(f"  V={V} (n={n}): safe-measure={float(meas):.4f} -> {'FLEXIBLE(dim1)' if meas>0 else 'RIGID(dim0)'}")
if __name__=='__main__': main()
