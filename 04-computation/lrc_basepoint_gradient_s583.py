#!/usr/bin/env python3
"""Observer is illusory: compute the LRC value M_p from EVERY basepoint p in the n-point
geometry P. M_p = M({p_q - p_p : q!=p}) = loneliness margin seen from point p. The
gradient M_p vs position; the observer M_0; is the extreme the hardest? round 2."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def Mval(V):  # V = nonzero speeds (a centered config), M = max_t min ||v t||
    V=[v for v in V if v!=0]
    if not V: return F(1,2)
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (abs(V[i]+V[j]),abs(V[i]-V[j])):
                if D:
                    for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*abs(V[i])))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def main():
    print("ROUND 2: M_p from every basepoint p (delta=1/n). gradient + which observer is hardest")
    examples={
      "AP n=7 {0..6}":[0,1,2,3,4,5,6],
      "AP n=9 {0..8}":[0,1,2,3,4,5,6,7,8],
      "n=7 {0,1,2,3,5,7,11}":[0,1,2,3,5,7,11],
      "n=8 {0,1,2,3,4,5,6,7} (AP)":[0,1,2,3,4,5,6,7],
      "n=6 {0,2,3,7,8,10}":[0,2,3,7,8,10],
    }
    for name,P in examples.items():
        P=sorted(set(P)); n=len(P); delta=F(1,n)
        Mps=[]
        for p in P:
            V=[q-p for q in P if q!=p]
            Mps.append(Mval(V))
        margins=[float(m-delta) for m in Mps]
        hardest=P[margins.index(min(margins))]
        print(f"  {name}: delta={float(delta):.3f}")
        print(f"     M_p - delta by basepoint: {[ (P[i],round(margins[i],4)) for i in range(n)]}")
        print(f"     observer(0) margin={margins[0]:+.4f}; hardest basepoint={hardest} (margin {min(margins):+.4f})")
if __name__=='__main__': main()
