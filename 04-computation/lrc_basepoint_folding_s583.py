#!/usr/bin/env python3
"""The gradient IS folding: basepoint p sees magnitudes {|q-p|}; equal magnitudes fold
(||-vt||=||vt||). Extreme observer sees n-1 distinct (one-sided, tight); center sees
~(n-1)/2 distinct (two-sided -> folded, loose). M_p tracks the distinct-magnitude count.
Ties to augmentation (one-sided=unbalanced) & 2-adic doubling. round 3."""
from fractions import Fraction as F
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def Mval(V):
    V=[abs(v) for v in V if v!=0]
    if not V: return F(1,2)
    cands=set()
    for v in V: cands.add(F(1,2*v))
    for a,b in combinations(V,2):
        for D in (a+b,abs(a-b)):
            if D:
                for m in range(1,D): cands.add(F(m,D))
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def main():
    print("ROUND 3: M_p vs #distinct magnitudes seen from basepoint p (AP {0..n-1})")
    for n in [7,8,9,14]:
        P=list(range(n)); delta=F(1,n); rows=[]
        for p in P:
            mags=sorted(set(abs(q-p) for q in P if q!=p))
            rows.append((p,len(mags),float(Mval([q-p for q in P if q!=p])-delta)))
        print(f"  n={n} (delta={float(delta):.3f}): (basepoint, #distinct-magnitudes, M_p-delta)")
        print(f"     {[(p,d,round(m,3)) for p,d,m in rows]}")
        # the center: half the magnitudes => M ~ 1/(ceil((n-1)/2)+1)
        c=rows[n//2]
        print(f"     center p={c[0]}: {c[1]} distinct magnitudes; M_center-delta={c[2]:+.3f}; "
              f"predicted loose level 1/({c[1]}+1)={1/(c[1]+1):.3f} vs delta={float(delta):.3f}")
if __name__=='__main__': main()
