#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: named families + the observer-inclusive variant.
monad-explorer-2026-06-06-S3.

(1) Gstar(AP) for AP = {1,2,...,n-1}: a concrete, named UPPER bound on the inf
    (and a check that the pure AP is NOT the minimizer / is above the inf).
(2) Gstar_full: the OBSERVER-INCLUSIVE mutual gap.  Pairs with the observer
    (speed 0) give the bare speed |v_i| (sign-invariant), so
        W_full(A,B) = {v_i : all i} u W(A,B).
    Gstar_full(S) = max over cuts M(W_full) <= min(M_obs(S), Gstar(S)).
    Reports inf_S Gstar_full(S) per n and compares to 1/n.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(Ws):
    T = set()
    for w in Ws:
        for a in range(0, w):
            T.add(F(2*a+1, 2*w))
    for x, y in combinations(Ws, 2):
        for d in (x+y, x-y):
            if d == 0: continue
            d = abs(d)
            for a in range(1, d):
                T.add(F(a, d))
    return T

def maximin(W):
    Ws = sorted(set(abs(w) for w in W if w != 0))
    if not Ws: return F(1,2)
    best = F(0)
    for t in candidate_times(Ws):
        ok = True; m = F(1,2)
        for w in Ws:
            nv = norm(w*t)
            if nv <= best: ok=False; break
            if nv < m: m = nv
        if ok and m > best: best = m
    return best

def cutW(V, side, obs=False):
    r=len(V); W=list(V) if obs else []
    for i,j in combinations(range(r),2):
        W.append(abs(V[i]-V[j]) if side[i]==side[j] else V[i]+V[j])
    return W

def gstar(V, obs=False):
    r=len(V); best=F(-1); bs=None
    for tail in product([0,1],repeat=r-1):
        side=(0,)+tail
        g=maximin(cutW(V,side,obs))
        if g>best: best=g; bs=side
    return best,bs

def enum_sets(r,B):
    return [c for c in combinations(range(1,B+1),r) if reduce(gcd,c)==1]

def main():
    print("="*72)
    print("(1) Gstar(AP={1..n-1}) — named upper bound on inf_S Gstar")
    print("="*72)
    for n in range(3,11):
        V=tuple(range(1,n))
        g,side=gstar(V)
        A=[V[i] for i in range(len(V)) if side[i]==0]
        B=[V[i] for i in range(len(V)) if side[i]==1]
        print(f"  n={n:2d} AP={V}: Gstar={g}={float(g):.5f}  (1/n={float(F(1,n)):.5f}, "
              f"1/(n-1)={float(F(1,n-1)):.5f})  best cut A={A}|B={B}")
    print("\n"+"="*72)
    print("(2) inf_S Gstar_full (observer-inclusive mutual gap) vs 1/n")
    print("="*72)
    cfg={2:14,3:12,4:10,5:9,6:8}
    for r in range(2,7):
        B=cfg[r]; n=r+1; floor=F(1,n)
        sets=enum_sets(r,B)
        inf_full=F(2); arg=None; below=0; tight=0
        for V in sets:
            g,_=gstar(V,obs=True)
            if g<inf_full: inf_full=g; arg=V
            if g<floor: below+=1
            if g==floor: tight+=1
        print(f"  n={n} (B<={B}, {len(sets)} sets): inf Gstar_full={inf_full}={float(inf_full):.5f} "
              f"(n*inf={float(n*inf_full):.4f}) at {arg}; below 1/n: {below}; tight: {tight}")

if __name__=="__main__":
    main()
