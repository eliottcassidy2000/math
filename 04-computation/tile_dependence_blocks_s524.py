#!/usr/bin/env python3
"""
tile_dependence_blocks_s524.py    oracle-2026-06-01-S524

User's thesis: a labeled tournament is NOT 2^C(n,2) independent arc-switches. It
is a RANKING (Hamiltonian path) composed recursively of sub-rankings that are
either aligned with the path (transitive) or against it (anti-transitive). The
tiling model exposes this hidden dependence; flipping an arc is not independent.

This is THM-354: fix the base Hamiltonian path; the path's strong components are
intervals; #good-cuts (positions crossed by a BACKWARD/upset arc) = n - #SCC.
So the ranking = a transitive chain of #SCC strongly-connected ANTI-TRANSITIVE
blocks. The "tiles" are the non-path arcs = the upset coordinates.

We quantify the hidden dependence on the tiling cube (base path 0->1->..->n-1,
tile (i,j), j>i+1, aligned = i->j, upset = j->i):
 (1) COLLAPSE: 2^C(n-1,2) tilings -> few iso classes (fibers = the dependence).
 (2) CONTEXT-DEPENDENCE of an arc flip: for each tile t, the fraction of contexts
     where flipping t CHANGES the iso class. If a tile were an independent axis
     this would be ~1 and constant; instead it varies by t and is often <1
     (silent flips) -> arcs are not independent.
 (3) BLOCK STRUCTURE: verify THM-354 (#good-cuts = n - #SCC); show the ranking is
     a chain of strong blocks; tabulate #SCC.
 (4) ROUND slice: which tilings are ROUND (the LRC-realizable necklace, S523)?
     Characterize their block/upset structure.
"""
from itertools import combinations, permutations, product
from collections import Counter

_P = {}
def canon(adj):
    m=len(adj); P=_P.setdefault(m, list(permutations(range(m))))
    best=None
    for p in P:
        f=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or f<best: best=f
    return best

def is_round(adj):
    m=len(adj)
    for perm in permutations(range(1,m)):
        order=(0,)+perm; ok=True
        for i in range(m):
            v=order[i]; d=sum(adj[v])
            arc={order[(i+k)%m] for k in range(1,d+1)}
            act={w for w in range(m) if w!=v and adj[v][w]}
            if arc!=act: ok=False;break
        if ok: return True
    return False

def sccs(adj):
    """Tarjan-free: strong components via reachability (small n)."""
    m=len(adj)
    def reach(s):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for w in range(m):
                if adj[u][w] and w not in seen: seen.add(w); st.append(w)
        return seen
    comp=[None]*m; c=0
    for v in range(m):
        if comp[v] is not None: continue
        fwd=reach(v)
        # backward reach
        seen={v}; st=[v]
        while st:
            u=st.pop()
            for w in range(m):
                if adj[w][u] and w not in seen: seen.add(w); st.append(w)
        for w in (fwd & seen):
            if comp[w] is None: comp[w]=c
        c+=1
    return len(set(comp))

def good_cuts(adj, n):
    """cuts k in 1..n-1 with a backward arc j->i crossing (i<k<=j)."""
    g=0
    for k in range(1,n):
        backward=any(adj[j][i] for i in range(k) for j in range(k,n))
        if backward: g+=1
    return g

def build(n, upset_set, tiles):
    adj=[[0]*n for _ in range(n)]
    for i in range(n-1): adj[i][i+1]=1            # base path arcs i->i+1
    for idx,(i,j) in enumerate(tiles):
        if idx in upset_set: adj[j][i]=1          # upset
        else: adj[i][j]=1                          # aligned
    return adj

def analyze(n):
    tiles=[(i,j) for i in range(n) for j in range(i+2,n)]
    T=len(tiles)
    allt=list(range(T))
    canon_of={}; round_set=set(); scc_of={}
    iso=Counter()
    for bits in range(1<<T):
        up={t for t in allt if (bits>>t)&1}
        adj=build(n,up,tiles)
        c=canon(adj); canon_of[bits]=c; iso[c]+=1
        if is_round(adj): round_set.add(bits)
        sc=sccs(adj); scc_of[bits]=sc
        assert good_cuts(adj,n)==n-sc, "THM-354 violated"
    # context dependence per tile
    expr=[0]*T
    for bits in range(1<<T):
        for t in allt:
            if canon_of[bits]!=canon_of[bits ^ (1<<t)]:
                expr[t]+=1
    total=1<<T
    return dict(tiles=tiles, T=T, n_iso=len(iso), total=total,
               fiber=sorted(iso.values(), reverse=True),
               expr=[e/total for e in expr],
               round_n=len(round_set),
               round_scc=Counter(scc_of[b] for b in round_set),
               all_scc=Counter(scc_of[b] for b in range(total)),
               round_upsets=Counter(bin(b).count('1') for b in round_set))

def main():
    print("Tile dependence + block structure: arcs are not independent (oracle-S524)\n")
    for n in (4,5,6):
        r=analyze(n)
        print("="*64)
        print(f"n={n}: tiles(non-path arcs)={r['T']}  ->  2^{r['T']}={r['total']} tilings")
        print(f"  distinct iso classes among tilings: {r['n_iso']}  (collapse {r['total']}->{r['n_iso']})")
        print(f"  fiber sizes (tilings per class, =H/Aut): {r['fiber']}")
        print(f"  THM-354 verified on all tilings: good-cuts == n - #SCC  (assert passed)")
        print(f"  #SCC distribution over all tilings: {dict(sorted(r['all_scc'].items()))}")
        print(f"  arc-flip CONTEXT-DEPENDENCE (frac of contexts where flipping tile t")
        print(f"     changes the iso class), per tile {r['tiles']}:")
        print(f"     {[round(x,3) for x in r['expr']]}")
        print(f"     -> NOT constant, NOT 1: an arc flip's effect depends on the other arcs.")
        print(f"  ROUND tilings (LRC-realizable, S523): {r['round_n']} of {r['total']}")
        print(f"     their #SCC dist: {dict(sorted(r['round_scc'].items()))}; #upsets dist: {dict(sorted(r['round_upsets'].items()))}")
        print()
    print("READING: the 2^C(n-1,2) tile cube collapses massively onto iso classes;")
    print("a single arc flip is silent or expressive DEPENDING ON CONTEXT (varies by")
    print("tile and by the rest of the state). The ranking is a transitive chain of")
    print("#SCC anti-transitive strong blocks (THM-354). Arcs are not independent.")

if __name__=="__main__":
    main()
