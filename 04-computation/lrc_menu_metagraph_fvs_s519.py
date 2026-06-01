#!/usr/bin/env python3
"""
lrc_menu_metagraph_fvs_s519.py  (oracle-2026-06-01-S519)

Try to prove LRC via the walk-on-the-metagraph concept, with a REFRAME:
 - The runner clock is a CLOSED walk in G_n through the circular-tournament menu
   (the 2*Fib(n-2) realizable classes, S518).
 - Adjacent clock cells differ by ONE wall-crossing = ONE arc flip = a G_n edge.
 - So a runner system = a CLOSED WALK in the "menu metagraph" M_n (menu classes,
   edges = single-arc-flip transitions that stay in the menu).
 - LONELY classes (here: transitive, the bunched/empty-semicircle = 1/2-gap end).
 REFRAME / KEY TEST: are the lonely classes a FEEDBACK VERTEX SET of M_n
 (does M_n minus the lonely classes become ACYCLIC)?  If yes, EVERY closed walk
 hits a lonely class -> LRC (at this resolution).  We test it.
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
def is_trans(A): return sorted(sum(r) for r in A)==list(range(len(A)))

def flips(A):
    """all single-arc-flip neighbors (reverse one edge)."""
    n=len(A); out=[]
    for i in range(n):
        for j in range(i+1,n):
            B=[list(r) for r in A]
            B[i][j],B[j][i]=B[j][i],B[i][j]
            out.append(tuple(map(tuple,B)))
    return out

def menu(n, samples=40000, seed=0):
    rng=random.Random(seed); raw=set()
    for _ in range(samples):
        raw.add(half_turn(sorted(rng.random() for _ in range(n))))
    cl={}
    for A in raw:
        c=canon(A)
        if c not in cl: cl[c]=A
    return cl  # canon -> representative adj

def main():
    for n in (5,6,7):
        cl=menu(n)
        cset=set(cl)
        # build metagraph: edge c1-c2 if some flip of rep(c1) canonicalizes to c2 in menu
        import collections
        adj=collections.defaultdict(set)
        for c,A in cl.items():
            for B in flips(A):
                cb=canon(B)
                if cb in cset and cb!=c:
                    adj[c].add(cb); adj[cb].add(c)
        V=list(cl); idx={c:i for i,c in enumerate(V)}
        E=sum(len(adj[c]) for c in V)//2
        lonely=[c for c in V if is_trans(cl[c])]
        # remove lonely, test acyclic (forest) on remaining
        rem=[c for c in V if c not in set(lonely)]
        remset=set(rem)
        # count components & edges among rem
        seen=set(); comps=0; cyc=False
        # union-find for forest test
        par={c:c for c in rem}
        def find(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        edges_rem=0
        for c in rem:
            for d in adj[c]:
                if d in remset and idx[c]<idx[d]:
                    edges_rem+=1
                    ra,rb=find(c),find(d)
                    if ra==rb: cyc=True
                    else: par[ra]=rb
        ncomp=len({find(c) for c in rem})
        forest = (edges_rem == len(rem)-ncomp)  # forest iff E = V - components
        Hs={H(cl[c]) for c in V}
        print(f"n={n}: menu |V|={len(V)} edges={E}  lonely(transitive)={len(lonely)}")
        print(f"   M_n minus lonely: |V|={len(rem)} edges={edges_rem} comps={ncomp}  "
              f"ACYCLIC(forest)={forest}  => lonely is FVS: {forest}")
        print(f"   H-values in menu: {sorted(Hs)}")
        # also: is the WHOLE menu metagraph connected? does it even have cycles?
        par2={c:c for c in V}
        def f2(x):
            while par2[x]!=x: par2[x]=par2[par2[x]]; x=par2[x]
            return x
        for c in V:
            for d in adj[c]:
                if idx[c]<idx[d]: par2[f2(c)]=f2(d)
        comp_all=len({f2(c) for c in V})
        print(f"   whole M_n: components={comp_all}, has_cycle={E>len(V)-comp_all}")
        print()

if __name__=="__main__": main()
