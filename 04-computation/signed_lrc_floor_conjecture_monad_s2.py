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

def maximin_gap_grid(W, mult=4):
    """Independent check: brute rational grid t=a/Q, Q=mult*lcm(2*|w|)."""
    W=[w for w in W if w!=0]
    if not W: return F(1,2)
    Q=1
    for w in W: Q=lcm(Q,2*abs(w))
    Q*=mult
    best=F(0)
    for a in range(1,Q):
        t=F(a,Q)
        m=min(norm(w*t) for w in W)
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
    print("="*78)
    print("HYP-2293: signed pairwise LRC floor  Gstar(S) >= 1/(r+1)   monad-S2c")
    print("="*78)
    # extended exhaustive sweep
    plan={2:14, 3:11, 4:9, 5:8, 6:7}
    for r in range(2,7):
        B=plan[r]; n=r+1; floor=F(1,n)
        sets=enum_sets(r,B)
        below=0; equal=[]; above=0
        worst=F(1)
        for V in sets:
            g,_=gstar(V)
            if g<floor: below+=1; print(f"  !! COUNTEREXAMPLE V={V}: Gstar={g} < {floor}")
            elif g==floor: equal.append(V)
            else: above+=1
            if g<worst: worst=g
        print(f"\n--- r={r} (n={n}), B<={B}, {len(sets)} sets, floor=1/n={floor} ---")
        print(f"  below floor: {below}   equal (TIGHT): {len(equal)}   above: {above}")
        print(f"  min Gstar over all sets: {worst}  (=floor: {worst==floor})")
        if equal:
            print(f"  tight sets (Gstar = 1/n exactly), first 12: {equal[:12]}")
    # independent exactness certification on a sample ----------------------
    print("\n"+"="*78)
    print("EXACTNESS CERTIFICATION: candidate-set gap  vs  brute rational grid")
    print("="*78)
    import random
    random.seed(12345)
    mism=0; checked=0
    sample=[]
    for r in [3,4,5]:
        ss=enum_sets(r, plan[r])
        random.shuffle(ss)
        sample += [(r,V) for V in ss[:15]]
    for r,V in sample:
        for tail in product([1,-1],repeat=r-1):
            eps=(1,)+tail
            D=signed_diff(V,eps)
            g1=maximin_gap(D); g2=maximin_gap_grid(D,mult=3)
            checked+=1
            if g1!=g2:
                mism+=1
                print(f"  MISMATCH V={V} eps={eps}: cand={g1} grid={g2}")
    print(f"  checked {checked} (set,eps) gaps; mismatches: {mism}  "
          f"({'EXACT method CONFIRMED' if mism==0 else 'METHOD BUG'})")

if __name__=="__main__":
    main()
