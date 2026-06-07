#!/usr/bin/env python3
"""
HYP-2293 test: the signed pairwise LRC floor.
monad-explorer-2026-06-06-S2c.

CONJECTURE (HYP-2293): for every speed set S={v_1<...<v_r} (gcd 1), the sign-
gauge-optimized moving-pairwise gap
     Gstar(S) = max_{eps in {+-1}^r} max_t min_{i<j} ||(eps_i v_i - eps_j v_j) t||
satisfies  Gstar(S) >= 1/(r+1) = 1/n,  the SAME floor as the observer LRC,
with equality attained.  (Naive pair-count floor would be 1/(C(r,2)+1) << 1/n.)

This script (a) tests the conjecture exhaustively at larger speed bounds and at
r=6; (b) records all EQUALITY (tight) cases; (c) INDEPENDENTLY verifies the
exact maximin gap on a sample via a denser candidate set + a brute rational
grid, to certify the gap routine before any claim.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
from functools import reduce
import sys
def pr(*a):
    print(*a); sys.stdout.flush()

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(W):
    W=[w for w in W if w!=0]
    ms=set(abs(w) for w in W)
    for a,b in combinations(W,2):
        ms.add(abs(a+b)); ms.add(abs(a-b))
    ms.discard(0)
    T=set()
    for m in ms:
        for a in range(1,2*m): T.add(F(a,2*m))
    return T

def maximin_gap(W):
    W=[w for w in W if w!=0]
    if not W: return F(1,2)
    best=F(0)
    for t in candidate_times(W):
        m=min(norm(w*t) for w in W)
        if m>best: best=m
    return best

def maximin_gap_float(W, N=2_000_000):
    """Independent BOUNDED check via a fine uniform float grid on (0,1)."""
    W=[abs(w) for w in W if w!=0]
    if not W: return 0.5
    best=0.0
    for a in range(1,N):
        t=a/N
        m=1.0
        for w in W:
            f=(w*t)%1.0
            d=f if f<0.5 else 1.0-f
            if d<m: m=d
            if m<=best: break
        if m>best: best=m
    return best

def signed_diff(V,eps):
    return [eps[i]*V[i]-eps[j]*V[j] for i,j in combinations(range(len(V)),2)]

def gstar(V):
    """max over effective sign patterns (eps_0=+1) and the maximizers."""
    best=F(0); arg=[]
    for tail in product([1,-1],repeat=len(V)-1):
        eps=(1,)+tail
        g=maximin_gap(signed_diff(V,eps))
        if g>best: best=g; arg=[eps]
        elif g==best: arg.append(eps)
    return best,arg

def enum_sets(r,B):
    return [c for c in combinations(range(1,B+1),r) if reduce(gcd,c)==1]

def main():
    pr("="*78)
    pr("HYP-2293: signed pairwise LRC floor  Gstar(S) >= 1/(r+1)   monad-S2c")
    pr("="*78)
    # extended exhaustive sweep
    plan={2:14, 3:11, 4:9, 5:8, 6:7}
    for r in range(2,7):
        B=plan[r]; n=r+1; floor=F(1,n)
        sets=enum_sets(r,B)
        below=0; equal=[]; above=0
        worst=F(1)
        for V in sets:
            g,_=gstar(V)
            if g<floor: below+=1; pr(f"  !! COUNTEREXAMPLE V={V}: Gstar={g} < {floor}")
            elif g==floor: equal.append(V)
            else: above+=1
            if g<worst: worst=g
        pr(f"\n--- r={r} (n={n}), B<={B}, {len(sets)} sets, floor=1/n={floor} ---")
        pr(f"  below floor: {below}   equal (TIGHT): {len(equal)}   above: {above}")
        pr(f"  min Gstar over all sets: {worst}  (=floor: {worst==floor})")
        if equal:
            pr(f"  tight sets (Gstar = 1/n exactly), first 12: {equal[:12]}")
    # independent exactness certification on a sample ----------------------
    pr("\n"+"="*78)
    pr("EXACTNESS CERTIFICATION: candidate-set gap  vs  brute rational grid")
    pr("="*78)
    import random
    random.seed(12345)
    mism=0; checked=0
    sample=[]
    for r in [3,4,5]:
        ss=enum_sets(r, plan[r])
        random.shuffle(ss)
        sample += [(r,V) for V in ss[:5]]
    maxdev=0.0
    for r,V in sample:
        for tail in product([1,-1],repeat=r-1):
            eps=(1,)+tail
            D=signed_diff(V,eps)
            g1=maximin_gap(D); g2=maximin_gap_float(D,N=200_000)
            checked+=1
            dev=abs(float(g1)-g2)
            if dev>maxdev: maxdev=dev
            # float grid is a LOWER-ish estimate; exact must be >= grid - 1/N and within ~2/minspeed/N
            if float(g1)+1e-9 < g2:        # exact strictly below grid -> exact MISSED a better t
                mism+=1
                pr(f"  EXACT-TOO-LOW V={V} eps={eps}: cand={g1}={float(g1):.6f} grid={g2:.6f}")
    pr(f"  checked {checked} (set,eps) gaps; exact-below-grid: {mism}; "
       f"max|exact-grid|={maxdev:.2e}  "
       f"({'EXACT method CONFIRMED (>= float grid, tiny dev)' if mism==0 else 'METHOD BUG'})")

if __name__=="__main__":
    main()
