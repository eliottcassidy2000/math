#!/usr/bin/env python3
"""
COVERING-MIN = n/Phi6(n)?  (mac-mini-2026-07-03-S30)
Conjecture: for LRC(n), the covering-min (min view M over covering (n-1)-families) = n/(n^2-n+1) = n/Phi6(n),
achieved by the tight family T_n = {1,...,n-2, (n-1)*n}, with the binding pair (1, (n-1)n) where
(n-1)*n ≡ -1 (mod Phi6(n))  [since (n-1)n = n^2-n = Phi6(n)-1 ≡ -1].  The apex speed (n-1)n = lcm(n-1,n)
covers the two hard moduli n-1, n with ONE runner; {1..n-2} covers 2..n-2. So T_n is covering, and the
binding pair 1 <-> -1 (mod Phi6(n)) mirrors, giving M = n/Phi6(n).
VERIFY: (a) M(T_n) = n/Phi6(n) exactly; (b) T_n is the covering-min (search small-speed covering families);
(c) the -1 mod Phi6(n) binding; (d) g (# gaps at optimum) <= 3.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd,xs)
def nd(x):
    x=x%1; return min(x,1-x)
def is_covering(sp,n): return all(any(v%q==0 for v in sp) for q in range(2,n+1))
def M_view(sp, D):
    best=F(0); bt=F(0)
    for a in range(1,D):
        t=F(a,D); m=min(nd(v*t) for v in sp)
        if m>best: best,bt=m,t
    return best,bt
def ngaps(sp,t):
    ph=sorted(set(float((v*t)%1) for v in sp)|{0.0})
    gs=set(round(ph[(i+1)%len(ph)]-ph[i] if i<len(ph)-1 else 1-ph[-1],6) for i in range(len(ph)))
    return len(gs)

if __name__=="__main__":
    print("COVERING-MIN(n) vs n/Phi6(n), Phi6(n)=n^2-n+1.  T_n = {1..n-2, (n-1)n}.")
    print("="*90)
    print(f"{'n':>3} {'Phi6(n)':>8} {'n/Phi6(n)':>12} {'M(T_n) exact':>14} {'match?':>7} {'apex mod Phi6':>13} {'g@opt':>6} {'search-min match?':>18}")
    for n in range(4,11):
        phi6=n*n-n+1
        target=F(n,phi6)
        apex=(n-1)*n
        Tn=list(range(1,n-1))+[apex]
        D=phi6*6
        M,t=M_view(Tn,D)
        apex_mod=apex%phi6
        g=ngaps(Tn,t)
        # search small-speed covering (n-1)-families for the min M (verify T_n is the covering-min)
        rng=random.Random(1000+n)
        minM=(F(1),None)
        cands=[Tn]
        hi=max(2*apex, 40)
        for _ in range(4000):
            k=n-1
            sp=sorted(set(rng.sample(range(1,min(hi,3*apex+1)),k)))
            if len(sp)==k: cands.append(sp)
        for sp in cands:
            if len(sp)!=n-1 or gcd_all(sp)!=1 or not is_covering(sp,n): continue
            m,_=M_view(sp, phi6*3)
            if m<minM[0]: minM=(m,sp)
        search_match = (minM[0]==target)
        print(f"{n:>3} {phi6:>8} {float(target):>12.6f} {str(M):>14} {str(M==target):>7} "
              f"{'-1' if apex_mod==phi6-1 else str(apex_mod):>13} {g:>6} {str(search_match)+' '+str(minM[1] if not search_match else ''):>18}")
    print("\n=> if M(T_n)==n/Phi6(n) for all n AND apex ≡ -1 (mod Phi6): the covering-min is the CYCLOTOMIC Phi6")
    print("   value, tight family {1..n-2,(n-1)n}, binding pair (1,-1) mod Phi6(n). 183=Phi6(14) explained.")
    print("   n=14: 14/183 > 1/14  <=>  183 < 14*14=196  <=>  Phi6(14)=183 < n^2  (covering margin exists).")
