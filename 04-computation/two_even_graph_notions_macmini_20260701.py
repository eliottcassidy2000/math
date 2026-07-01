#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S82 -- TWO DISTINCT "EVEN GRAPH" NOTIONS (resolving the repo's equinumerosity confusion).

The repo's tournaments-and-even-graphs.md flagged "A000568 = #even graphs NOT confirmed (4!=3)". Resolution:
there are TWO different "even graph" definitions, and the repo used the wrong one for the equinumerosity.

 (I)  EVEN-DEGREE graph (repo's E_n): every vertex has even degree = the cycle space of K_n.
      Count up to iso = A002854 = 1,2,3,7,16,54 (= two-graphs = switching classes, Mallows-Sloane 1975).

 (II) AUTOMORPHISM-PARITY even graph (Royle-Praeger-Glasby-Freedman-Devillers 2022, arXiv:2204.01947):
      fix a reference orientation (i->j for i<j). An automorphism g of graph X REVERSES edge {u,v} (u<v)
      iff u^g > v^g. X is ODD if some g in Aut(X) reverses an ODD number of edges, EVEN otherwise.
      THEOREM 1.1: #even graphs on n = #tournaments on n = A000568 (2,4,12,56,456).
      Mechanism: #graphs = #tournaments + #odd graphs (Cauchy-Frobenius).

This script verifies (II) gives A000568 and #graphs = #tournaments + #odd, for n=3,4,5(,6).
"""
import itertools
from math import comb

def graph_canon(edges_mask, n, idx, perms):
    best = None
    for p in perms:
        m = 0
        for (i, j), k in idx.items():
            if (edges_mask >> k) & 1:
                a, b = (p[i], p[j]) if p[i] < p[j] else (p[j], p[i])
                m |= (1 << idx[(a, b)])
        if best is None or m < best: best = m
    return best

def reverses_odd(edges_mask, g, idx):
    """# edges {u,v} (u<v) in the graph with u^g > v^g, taken mod 2 -- is it odd for automorphism g?"""
    cnt = 0
    for (i, j), k in idx.items():
        if (edges_mask >> k) & 1:
            if g[i] > g[j]: cnt ^= 1        # parity of reversed edges
    return cnt  # 1 if odd

def is_automorphism(edges_mask, g, idx, n):
    for (i, j), k in idx.items():
        bit = (edges_mask >> k) & 1
        a, b = (g[i], g[j]) if g[i] < g[j] else (g[j], g[i])
        if ((edges_mask >> idx[(a, b)]) & 1) != bit: return False
    return True

A000568 = {3:2, 4:4, 5:12, 6:56}
A000088 = {3:4, 4:11, 5:34, 6:156}   # graphs up to iso

print(f"{'n':>2} {'#graphs':>8} {'#tournaments':>12} {'#even(II)':>10} {'#odd':>6} {'even=tourn?':>12} {'#graphs=t+odd?':>14}")
for n in range(3, 6):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    P = len(pairs); perms = list(itertools.permutations(range(n)))
    # enumerate graph iso classes
    reps = {}
    for mask in range(1 << P):
        c = graph_canon(mask, n, idx, perms)
        if c not in reps: reps[c] = mask
    ngraphs = len(reps)
    neven = nodd = 0
    for c, mask in reps.items():
        auts = [g for g in perms if is_automorphism(mask, g, idx, n)]
        odd = any(reverses_odd(mask, g, idx) for g in auts)
        if odd: nodd += 1
        else: neven += 1
    print(f"{n:>2} {ngraphs:>8} {A000568[n]:>12} {neven:>10} {nodd:>6} "
          f"{str(neven==A000568[n]):>12} {str(ngraphs==A000568[n]+nodd):>14}")

print("\n=> (II) EVEN graphs (no automorphism reverses an odd #edges) = tournaments = A000568.")
print("   The repo's E_n (A002854 even-DEGREE / switching classes) is a DIFFERENT object -- do not conflate.")
print("   Equinumerosity mechanism: #graphs = #tournaments + #odd graphs (Cauchy-Frobenius; RPGFD 2022).")
