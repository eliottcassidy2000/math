#!/usr/bin/env python3
"""
S632 — when ties (resonances) are deleted, a tournament becomes a partial oriented graph.
WHICH oriented graphs appear in our proofs?  Claim: the TIE graph (binding/resonant pairs) is a
PROXIMITY graph of the residues on the underlying abelian group:
  - LRC: runners at residues r_i = v_i*a (mod m) on the discrete circle Z/m; tie iff circular
    distance <= band  ==>  a CIRCULAR-ARC / CIRCULANT graph (the discrete unit-distance graph on Z/m).
  - unit distance: points in R^2, tie iff |p_i-p_j|=1  ==>  the planar unit-distance (Cayley) graph.
So both tie graphs are UNIT-DISTANCE / proximity graphs of the ambient group; the relevant invariants
are chromatic number (Hadwiger-Nelson / sieve coloring) and independence (the lonely/corrector set).
"""
from fractions import Fraction as Fr
from math import gcd
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_brute

def norm(x):
    f = x - (x.numerator // x.denominator); return f if f <= Fr(1, 2) else 1-f

def tie_graph_at_witness(speeds, n):
    """at the optimal witness t*=a/m, edge (i,j) if pair is BINDING: ||(v_i-v_j)t*|| <= gap.
       include observer 0. returns (m,a, adjacency, residues)."""
    V = [0]+list(speeds); gap = gap_brute(speeds)
    # find a witness t* achieving the gap (search shells)
    best = None
    for m in range(2, 4*n):
        for a in range(1, m):
            if gcd(a, m) != 1: continue
            md = min(norm(Fr(a*v, m)) for v in V if v != 0)  # observer gap = min over runners
            if md == gap:
                best = (m, a); break
        if best: break
    if not best: return None
    m, a = best
    res = [(a*v) % m for v in V]
    band = m*gap                      # tie threshold in residue units (gap = band/m)
    N = len(V); adj = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(i+1, N):
            d = (res[i]-res[j]) % m; d = min(d, m-d)
            if d <= band + 1e-9:      # binding/tie pair
                adj[i][j] = adj[j][i] = 1
    return m, a, adj, res, gap

def describe(adj):
    N = len(adj); deg = [sum(r) for r in adj]
    edges = sum(deg)//2
    # circulant test on the residue order is implicit; report degree sequence + connectivity + chromatic (greedy ub)
    # greedy chromatic
    order = sorted(range(N), key=lambda v: -deg[v]); color = {}
    for v in order:
        used = {color[u] for u in range(N) if adj[v][u] and u in color}
        c = 0
        while c in used: c += 1
        color[v] = c
    chi_ub = max(color.values())+1
    from collections import Counter
    return edges, dict(sorted(Counter(deg).items())), chi_ub

if __name__ == "__main__":
    print("LRC TIE GRAPHS (binding pairs at the optimal witness, incl. observer)")
    print("="*70)
    cases = [("AP n=5 {1,2,3,4}", [1,2,3,4], 5),
             ("AP n=6 {1..5}", [1,2,3,4,5], 6),
             ("AP n=7 {1..6}", [1,2,3,4,5,6], 7),
             ("sporadic tight (1,3,4,7) n=5", [1,3,4,7], 5),
             ("loose (1,2,4,8) n=5", [1,2,4,8], 5)]
    for name, S, n in cases:
        r = tie_graph_at_witness(S, n)
        if not r: print(f"  {name}: no shell witness found"); continue
        m, a, adj, res, gap = r
        E, degs, chi = describe(adj)
        # is it the cycle C_{n} (each vertex deg 2)?
        iscycle = all(d == 2 for d in [sum(row) for row in adj])
        print(f"  {name}: gap={gap} witness a/m={a}/{m} residues(incl obs)={sorted(res)}")
        print(f"     TIE graph: {E} edges, degree-dist {degs}, chromatic(ub)={chi}, is-cycle={iscycle}")
    print("\nINTERPRETATION:")
    print("  AP tie graph = the CYCLE C_n (consecutive runners on the clock) = circulant <±1>.")
    print("  General config tie graph = circular-arc proximity graph of residues on Z/m (a CIRCULANT")
    print("  when residues are an AP). This is the discrete unit-distance graph on the circle Z/m.")
    print("  Relevant invariants: chromatic number (= Hadwiger-Nelson on Z/m / the sieve coloring)")
    print("  and independence number (= the simultaneously-lonely set = the corrector).")
