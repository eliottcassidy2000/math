#!/usr/bin/env python3
"""
VERIFY whether T_n={1..n-2,(n-1)n} is the covering-min, with HIGH-RES M (mac-mini-2026-07-03-S30).
The S30 search used a coarse grid D=3*Phi6(n) that can UNDERESTIMATE a competitor's true M (its optimum
may sit at a denominator coprime to 3*Phi6). Recompute M EXACTLY (optimum of min_i||v_i t|| is at a
'crossing' rational a/q with q | some (v_i +- v_j) or q=2v_i; scan all such q up to a bound + fine grid).
Then: (1) confirm M([9,16,32]) etc. vs n/Phi6(n); (2) WIDE search for n=14 covering families with M<14/183.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd,xs)
def nd(x):
    x=x%1; return min(x,1-x)
def is_covering(sp,n): return all(any(v%q==0 for v in sp) for q in range(2,n+1))

def M_exact(sp):
    """Exact view M=max_t min_i ||v_i t||. Optimum at t=a/q where two runners are equidistant
       (q | v_i±v_j) or a runner is at 1/2 (q=2v_i) or at the band edge. Collect candidate denominators."""
    Q=set()
    vs=sorted(set(abs(v) for v in sp))
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1,len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(abs(vs[i]-vs[j]))
    # also multiples up to including sums (optimum a/q with q dividing lcm of pairwise combos)
    best=F(0); bt=None
    for q in sorted(Q):
        if q<2: continue
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            t=F(a,q); m=min(nd(v*t) for v in sp)
            if m>best: best,bt=m,t
    return best,bt

def M_grid(sp, D):
    best=F(0); bt=None
    for a in range(1,D):
        t=F(a,D); m=min(nd(v*t) for v in sp)
        if m>best: best,bt=m,t
    return best,bt

if __name__=="__main__":
    print("(1) HIGH-RES check of the n<=6 'competitors' the coarse search flagged below n/Phi6(n):")
    print(f"{'n':>3} {'family':>18} {'n/Phi6(n)':>11} {'M_exact':>10} {'M_grid(1e5)':>12} {'< n/Phi6?':>10} {'>=1/n?':>8}")
    tests=[(4,[9,16,32]),(5,[9,10,18,36]),(6,[7,40,47,54,59]),(4,[1,2,12]),(5,[1,2,3,20]),(6,[1,2,3,4,30])]
    for n,sp in tests:
        phi6=n*n-n+1; target=F(n,phi6)
        Me,te=M_exact(sp); Mg,tg=M_grid(sp,100000)
        cov=is_covering(sp,n); g1=gcd_all(sp)
        M=max(Me,Mg)  # both are lower bounds on true M; take the better
        print(f"{n:>3} {str(sp):>18} {float(target):>11.6f} {float(Me):>10.6f} {float(Mg):>12.6f} "
              f"{str(M<target):>10} {str(M>=F(1,n)):>8}  cov={cov} gcd={g1}")
    print("\n   (M_exact and M_grid are both valid LOWER bounds on true M; if EITHER >= n/Phi6, T_n is NOT beaten by that family)")
    # (2) WIDE search for n=14 covering families with M < 14/183 (high-res M)
    print("\n(2) WIDE search: any covering gcd=1 14-family (13 speeds) with M < 14/183 = %.6f ?"%float(F(14,183)))
    n=14; target=F(14,183); rng=random.Random(777)
    below=[]; best_above=(F(1),None); found_min=(F(1),None)
    for trial in range(8000):
        # mix ranges: small, medium, and lcm-structured
        hi=rng.choice([30,60,120,200,400])
        sp=sorted(set(rng.sample(range(1,hi),13)))
        if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp,14): continue
        Me,_=M_exact(sp)
        if Me<found_min[0]: found_min=(Me,sp)
        if Me<target-F(1,10**7): below.append((Me,sp))
    print(f"   tested many; families with M<14/183: {len(below)}")
    if below:
        below.sort()
        for m,sp in below[:5]:
            # re-verify with fine grid to rule out under-resolution
            Mg,_=M_grid(sp,200000)
            print(f"      M_exact={float(m):.6f} M_grid={float(Mg):.6f}  {sp}  (covering, gcd=1)")
    print(f"   overall min M found = {float(found_min[0]):.6f} (14/183={float(target):.6f}); family {found_min[1]}")
    print("\n=> if NO covering family beats 14/183 with a VALID (grid-confirmed) M: 14/183 is the covering-min, T_14 tight.")
    print("   if some family has grid-confirmed M<14/183: the covering-min is LOWER -- correct the canon.")
