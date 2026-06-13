#!/usr/bin/env python3
"""
arc_dependence_tiling_lrc_s520o.py  (oracle-2026-06-01-S520o)

Test the claim: arc-flips are NOT independent; the essence is a recursive tree of
sub-rankings (tiling model), and the LRC-realizable tournaments are a thin
DEPENDENT slice, not an arc-subcube.

(1) DEPENDENCE LAYERS: 2^C(n,2) labeled  ->  A000568(n) iso  ->  2*Fib(n-2) circular.
(2) ARC-FLIP IS NOT LOCAL: for circular (round, runner-realizable) tournaments,
    what fraction of single-arc-flips STAY circular?  (low => the realizable set
    is not an arc-subcube => arcs are dependent.)
(3) H = #Hamiltonian paths = #rankings the arcs support: quantifies the hidden
    dependence (arcs shared among H rankings). Show H over the circular menu.
"""
import random
from itertools import permutations, combinations
from functools import lru_cache

def half_turn(pts):
    n=len(pts); A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and 0<(pts[i]-pts[j])%1.0<0.5: A[i][j]=1
    return tuple(map(tuple,A))
def canon(A):
    n=len(A); best=None
    for p in permutations(range(n)):
        f=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best
def H(A):
    n=len(A); full=(1<<n)-1
    @lru_cache(None)
    def dp(m,l):
        if m==full: return 1
        return sum(dp(m|(1<<x),x) for x in range(n) if not (m>>x)&1 and A[l][x])
    return sum(dp(1<<s,s) for s in range(n))
def is_round(A):
    """round/local: exists cyclic order where each out-set is a consecutive arc."""
    n=len(A)
    for p in permutations(range(n)):
        ok=True
        for a in range(n):
            outs=[ (b-a)%n for b in range(n) if a!=b and A[p[a]][p[b]] ]
            outs=sorted(outs)
            # consecutive run starting at 1?
            if outs and outs!=list(range(1,len(outs)+1)): ok=False; break
        if ok: return True
    return False
def menu(n, s=40000, seed=0):
    rng=random.Random(seed); raw=set()
    for _ in range(s): raw.add(half_turn(sorted(rng.random() for _ in range(n))))
    reps={}
    for A in raw:
        c=canon(A)
        if c not in reps: reps[c]=A
    return reps
def flips(A):
    n=len(A)
    for i in range(n):
        for j in range(i+1,n):
            B=[list(r) for r in A]; B[i][j],B[j][i]=B[j][i],B[i][j]
            yield tuple(map(tuple,B))

A000568={3:2,4:4,5:12,6:56,7:456}
fib=[1,1,2,3,5,8,13]
def main():
    print("Arc dependence / tiling / LRC (oracle-2026-06-01-S520o)\n")
    print("(1) DEPENDENCE LAYERS:")
    print(" n | 2^C(n,2) labeled | A000568 iso | 2*Fib(n-2) circular")
    for n in (3,4,5,6,7):
        print(f" {n} | {2**(n*(n-1)//2):>10} | {A000568[n]:>6} | {2*fib[n-2]:>4}")
    print("   => the 'independent arc' count is astronomically larger than the")
    print("      essence (iso classes), and the runner-realizable slice is Fibonacci-thin.\n")

    print("(2) ARC-FLIP IS NOT LOCAL (fraction of single-arc-flips staying circular):")
    for n in (5,6,7):
        reps=menu(n)
        circ=set(reps)
        tot=0; stay=0
        for c,A in reps.items():
            for B in flips(A):
                tot+=1
                if canon(B) in circ: stay+=1
        print(f" n={n}: {stay}/{tot} = {stay/tot:.3f} of single-arc-flips of a circular")
        print(f"        tournament stay circular  (C(n,2)={n*(n-1)//2} arcs/class)")
    print("   => the realizable (round) set is NOT closed under independent arc-flips:")
    print("      it is a thin dependent slice, not a coordinate subcube. A flip mostly")
    print("      throws you OUT of the realizable family -> arcs are correlated.\n")

    print("(3) H = #rankings (Hamiltonian paths) the arcs support = the dependence:")
    for n in (5,6,7):
        reps=menu(n)
        Hs=sorted(set(H(A) for A in reps.values()))
        print(f" n={n}: circular menu H-values = {Hs}  (H counts the rankings sharing the arcs)")
    print("   => the same arc-set supports H different Ham-path rankings; they share arcs,")
    print("      so a flip perturbs ALL of them at once. H is exactly this hidden coupling")
    print("      (and was the S26 loneliness meter).")

if __name__=="__main__": main()
