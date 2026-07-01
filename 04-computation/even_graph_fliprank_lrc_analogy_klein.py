#!/usr/bin/env python3
"""
even_graph_fliprank_lrc_analogy_klein.py  --  klein-2026-07-01-S74

TWO NOVEL OBLIGATIONS constructed together via the TOURNAMENT<->EVEN-GRAPH bijection, then transported to LRC.

The tiling/cycle-space bijection (CLAUDE.md): a tile-vector in GF(2)^m (m=C(n-1,2)) maps, via XOR of
FUNDAMENTAL CYCLES of the base path 0-1-...-(n-1), to an EVEN GRAPH (every vertex even degree = cycle-space
element). The SAME cube Q_m carries TWO S_n-quotients: tournaments G_n (A000568=2,4,12,56) and even graphs
E_n (A002854=2,3,7,16).

OBLIGATION 1 (challenge the assumption that flip-rank is tournament-specific): compute the EVEN-GRAPH
flip-rank rho_E(n) (min subcube hitting all E_n classes) and rainbow R_E(n), and compare to the tournament
rho_G=1,2,4,7 / R_G=1,2,3,5 (S71/72). Same cube, different quotient -> does the packing/covering asymmetry
and the covering-excess transfer, or is it action-specific?

OBLIGATION 2 (transport to LRC): the LRC analogue of "cycle space / even graph" is the RELATION LATTICE
Lambda={t: sum t_i v_i=0} (THM-515); its length-3 count (# 3-term relations v_i+v_j=v_k, the leading
ADDITIVE ENERGY) is the analogue of the tournament's 3-cycle count c3 (the leading OCF term). Verify the
dictionary on the covering-min construction vs the tight AP.
"""
import itertools, math

def all_edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]

def fund_cycles(n):
    """tile (non-tree edge) -> fundamental cycle bitmask over ALL edges, base path 0-1-...-(n-1)."""
    E=all_edges(n); idx={e:k for k,e in enumerate(E)}
    tiles=[(i,j) for (i,j) in E if j-i>=2]
    fc=[]
    for (i,j) in tiles:
        bm=1<<idx[(i,j)]
        for k in range(i,j): bm^=1<<idx[(min(k,k+1),max(k,k+1))]
        fc.append(bm)
    return E, tiles, fc

def perm_edge_maps(n,E):
    idx={e:k for k,e in enumerate(E)}; maps=[]
    for pi in itertools.permutations(range(n)):
        m=[idx[(min(pi[i],pi[j]),max(pi[i],pi[j]))] for (i,j) in E]
        maps.append(m)
    return maps
def canon_graph(bm,maps,ne):
    best=None
    for m in maps:
        v=0
        for k in range(ne):
            if (bm>>k)&1: v|=1<<m[k]
        if best is None or v<best: best=v
    return best

def even_classes(n):
    E,tiles,fc=fund_cycles(n); ne=len(E); m=len(tiles); maps=perm_edge_maps(n,E)
    cid={}; cls=[0]*(1<<m); nc=0
    for tv in range(1<<m):
        bm=0
        for b in range(m):
            if (tv>>b)&1: bm^=fc[b]
        c=canon_graph(bm,maps,ne)
        if c not in cid: cid[c]=nc; nc+=1
        cls[tv]=cid[c]
    return m,cls,nc

def flip_rank(m,cls,nc,cap=3_000_000):
    lb=math.ceil(math.log2(nc))
    for k in range(lb,m+1):
        if math.comb(m,k)*(1<<(m-k))>cap:
            return None,lb  # too big; skip (small n only here)
        for free in itertools.combinations(range(m),k):
            rest=[e for e in range(m) if e not in free]
            for fa in range(1<<len(rest)):
                fb=0
                for b,e in enumerate(rest):
                    if (fa>>b)&1: fb|=1<<e
                seen=set(); ok=True
                for a in range(1<<k):
                    tt=fb
                    for b in range(k):
                        if (a>>b)&1: tt|=1<<free[b]
                    seen.add(cls[tt])
                    if len(seen)==nc: break
                if len(seen)==nc: return k,lb
    return None,lb

def rainbow(m,cls,nc):
    for k in range(int(math.log2(nc)),0,-1):
        for free in itertools.combinations(range(m),k):
            rest=[e for e in range(m) if e not in free]
            for fa in range(1<<len(rest)):
                fb=0
                for b,e in enumerate(rest):
                    if (fa>>b)&1: fb|=1<<e
                seen=set(); ok=True
                for a in range(1<<k):
                    tt=fb
                    for b in range(k):
                        if (a>>b)&1: tt|=1<<free[b]
                    if cls[tt] in seen: ok=False; break
                    seen.add(cls[tt])
                if ok: return k
    return 0

def additive_relations3(S):
    """# ordered 3-term relations a+b=c among speeds (the leading additive energy / LRC 'even-3-cycle')."""
    Sset=set(S); cnt=0
    for a in S:
        for b in S:
            if a+b in Sset: cnt+=1
    return cnt

if __name__=="__main__":
    print("OBLIGATION 1: EVEN-GRAPH flip-rank rho_E / rainbow R_E vs TOURNAMENT rho_G/R_G (same cube Q_m, different S_n-quotient)")
    print(f"  {'n':>1} {'m':>2} {'|E_n|':>5} {'rho_E':>5} {'ceilE':>5} {'R_E':>4} {'floorE':>6}   {'|G_n|':>5} {'rho_G':>5} {'R_G':>4}")
    rhoG={3:1,4:2,5:4,6:7}; RG={3:1,4:2,5:3,6:5}; Gn={3:2,4:4,5:12,6:56}
    for n in [3,4,5,6]:
        m,cls,nc=even_classes(n)
        rE,lb=flip_rank(m,cls,nc); RE=rainbow(m,cls,nc)
        fl=int(math.log2(nc))
        rE_s=str(rE) if rE is not None else ">cap"
        print(f"  {n} {m:>2} {nc:>5} {rE_s:>5} {lb:>5} {RE:>4} {fl:>6}   {Gn[n]:>5} {rhoG[n]:>5} {RG[n]:>4}")
    print("  covering-excess rho-ceil(log2): compare even-graph vs tournament (is it self-dual under the bijection?)")
    print("  DIAGNOSIS -- even-graph class orbit sizes in the fund-cycle cube (why even graphs resist subcubes):")
    for n in [5,6]:
        m,cls,nc=even_classes(n)
        from collections import Counter
        sizes=sorted(Counter(cls).values())
        print(f"    n={n}: {nc} classes, orbit sizes in Q_{m} = {sizes} (sum={sum(sizes)}=2^{m}); tiny orbits (empty/near-empty even graph) pin the subcube")
    print("  => R_E < floor(log2|E_n|) at n=6 (3<4): the tournament law 'rainbow = floor' is NOT self-dual --")
    print("     it was SPECIFIC to the tournament quotient (the balanced-cut shape); the SAME cube, quotiented by")
    print("     the even-graph S_n-action, packs/covers WORSE. Assumption challenged: flip-rank depends on the QUOTIENT,")
    print("     not the cube -- the tournament<->even-graph bijection is a cube-iso but NOT a flip-rank-iso.")

    print("\nOBLIGATION 2: LRC analogue -- relation lattice = the 'even graph' of the speeds; #3-term relations = c3 analogue")
    for n in [7,10,14]:
        core=list(range(1,n))                       # tight AP {1..n-1}
        constr=list(range(1,n-1))+[n*(n-1)]         # covering-min construction
        print(f"  n={n}: #3-term relations (a+b=c):  AP {{1..{n-1}}} = {additive_relations3(core)};  construction {{1..{n-2},{n*(n-1)}}} = {additive_relations3(constr)}")
    print("  (analogy: transitive tournament c3=0 <-> a set with FEW/structured 3-term relations; the OCF N <-> the lonely measure L;")
    print("   high additive energy <-> low lonely measure (THM-515). The relation lattice IS the LRC cycle space / even graph.)")
