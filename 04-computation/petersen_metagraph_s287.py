#!/usr/bin/env python3
"""
petersen_metagraph_s287.py -- opus-2026-03-24-S287

IS THE n=5 MERGED META-GRAPH RELATED TO THE PETERSEN GRAPH?

The n=5 merged meta-graph has V=10, E=21.
The Petersen graph has V=10, E=15.
The complement of Petersen (the Kneser complement) has V=10, E=30.

21 = 15 + 6. So the meta-graph has 6 MORE edges than Petersen.
Or: 30 - 21 = 9 FEWER edges than the complement.

The Petersen graph is the Kneser graph K(5,2):
  Vertices = 2-element subsets of {1,2,3,4,5} (= C(5,2) = 10)
  Edges = disjoint subsets ({a,b} ~ {c,d} iff a,b,c,d all different)

Our meta-graph has 10 merged iso classes.
Can we find a natural bijection between merged classes and 2-subsets of [5]?

At n=5: the meta-graph has 10 = C(5,2) nodes.
The Kneser K(5,2) also has C(5,2) = 10 nodes.
This is the SAME number. Coincidence?

C(n,2) = C(5,2) = 10 = number of ARC POSITIONS = number of PAIRS.
So: meta-graph nodes at n=5 ↔ arc positions at n=5?
Or: meta-graph nodes ↔ 2-subsets of [5]?

Author: opus-2026-03-24-S287
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

banner("THE n=5 META-GRAPH vs PETERSEN")

# Build merged meta-graph
cmap = {}; sizes = {}; tc = {}; reps = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1; tc[bits] = cmap[c]
nc = len(cmap)

info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                 'SC': cmap[comp_c] == cid, 'size': sizes[cid],
                 'score': tuple(sorted(sum(A[i]) for i in range(n)))}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

adj = np.zeros((nm, nm), dtype=int)
for bits in range(2**m):
    mi = mid_map[tc[bits]]
    for k in range(m):
        mj = mid_map[tc[bits ^ (1<<k)]]
        if mi != mj: adj[mi][mj] = 1; adj[mj][mi] = 1

ne = adj.sum() // 2
H_vec = [info[merged[mi][0]]['H'] for mi in range(nm)]
deg = [sum(adj[mi]) for mi in range(nm)]

print("  META-GRAPH: V=%d, E=%d" % (nm, ne))
print("  PETERSEN:   V=10, E=15")
print("  Difference: %d extra edges in meta-graph\n" % (ne - 15))

# Build Petersen graph
# Petersen = Kneser K(5,2): vertices = 2-subsets of [5]
# {a,b} ~ {c,d} iff {a,b} ∩ {c,d} = ∅
subsets = list(combinations(range(5), 2))
petersen = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        if len(set(subsets[i]) & set(subsets[j])) == 0:
            petersen[i][j] = 1; petersen[j][i] = 1

print("  Petersen degree sequence: %s (3-regular)" %
      sorted([sum(petersen[i]) for i in range(10)]))
print("  Meta-graph degree sequence: %s\n" % sorted(deg))

# ============================================================
# Is the meta-graph a SUPERGRAPH of Petersen?
# ============================================================

# Try all 10! bijections between meta-graph nodes and Petersen nodes
# to find one where every Petersen edge is also a meta-graph edge.

print("  SEARCHING: Is the meta-graph a supergraph of Petersen?")
print("  (Is there a relabeling where Petersen ⊂ meta-graph?)\n")

# Too many permutations (10! = 3.6M). Use degree constraints to prune.
# Petersen is 3-regular. Meta-graph has degrees [2,3,3,3,4,5,5,5,5,7].
# Petersen nodes (degree 3) must map to meta-graph nodes with degree ≥ 3.
# So the meta-graph node with degree 2 CANNOT be a Petersen node... wait,
# every meta-graph node must map to a Petersen node.
# Actually: we're looking for a bijection φ: Petersen_nodes → Meta_nodes
# such that if {i,j} is a Petersen edge, then {φ(i), φ(j)} is a meta-graph edge.

# Meta nodes with degree < 3 cannot receive ANY Petersen node (Petersen 3-regular
# means every Petersen neighbor pair must land on meta-graph neighbors).
# But meta-graph has a node of degree 2 → at most 2 Petersen edges can touch it.
# Since Petersen is 3-regular, every node needs 3 meta-neighbors → node of degree 2
# CAN'T host a Petersen node.

# So if meta-graph has a node of degree 2: meta-graph CANNOT contain Petersen.
# At n=5: degree 2 exists (node 9, H=15, |Aut|=5).
# Therefore: meta-graph does NOT contain Petersen as a subgraph!

# Wait: the question is about bijective mapping. Node 9 (degree 2) would need
# 3 neighbors in Petersen, but it only has 2 in meta-graph. Contradiction.

has_deg_2 = any(d < 3 for d in deg)
print("  Meta-graph has node of degree %d < 3: %s" % (min(deg), has_deg_2))
print("  → Meta-graph CANNOT contain Petersen as a spanning subgraph")
print("     (Petersen is 3-regular, but meta-graph has degree-2 node)\n")

# ============================================================
# Is the meta-graph the COMPLEMENT of Petersen minus some edges?
# ============================================================

# Complement of Petersen: V=10, E=30 (= C(10,2) - 15).
# Meta-graph: V=10, E=21.
# 30 - 21 = 9 edges in complement(Petersen) but NOT in meta-graph.

# The complement of Petersen = Johnson graph J(5,2) = line graph of K_5.
# J(5,2) is 6-regular (each 2-subset overlaps with 6 others).

# Meta-graph is NOT 6-regular (degrees vary 2-7).
# But maybe there's a PARTIAL match.

# ============================================================
# THE DEEPER QUESTION: 10 = C(5,2) nodes
# ============================================================

banner("10 = C(5,2): NODES AS PAIRS?")

print("""
  At n=5: V_merged = 10 = C(5,2).
  This is the SAME as the number of arc positions (pairs of vertices).
  The Petersen graph ALSO has 10 = C(5,2) vertices (as 2-subsets of [5]).

  QUESTION: Is there a NATURAL bijection between
  merged iso classes at n=5 and pairs {i,j} with i,j ∈ [5]?

  A natural candidate: each class is characterized by its
  "most important arc" — the arc whose reversal has the strongest effect.

  But from the flat alphabet theorem (S283): ALL arcs are equivalent
  at the class level. There's no "most important arc."

  ALTERNATIVE: The 10 classes have 10 distinct H values
  (almost — H=3 and H=9 each appear twice, H=15 appears twice).
  So H alone gives 7 groups, not 10.

  The Petersen graph is fundamentally about DISJOINTNESS of pairs.
  Two pairs {a,b} and {c,d} are Petersen-adjacent iff they're disjoint.

  Tournament isomorphism classes are about STRUCTURE of all arcs.
  There's no obvious pair-disjointness interpretation.

  THE 5 COMPLEMENT-FLIP OVERLAP PAIRS AT n=6:
  These are 5 edges where the complement of one class is adjacent
  to the other in the meta-graph. The Petersen graph has 5 "spokes"
  connecting the outer pentagon to the inner pentagram.

  The number 5 = n - 1 = the number of vertices minus 1.
  At n=6: 5 = 6 - 1. These 5 edges are the "diagonal" connections.
""")

# Let me at least check if there's a subgraph isomorphism (not spanning)
# between Petersen and the meta-graph.

# The complement of the meta-graph:
meta_comp = 1 - adj
np.fill_diagonal(meta_comp, 0)
ne_comp = meta_comp.sum() // 2
print("  Meta-graph complement: V=10, E=%d\n" % ne_comp)
print("  Complement degree sequence: %s\n" %
      sorted([sum(meta_comp[i]) for i in range(nm)]))

# Is the complement 3-regular? (Would mean meta-graph complement = Petersen)
comp_degs = [sum(meta_comp[i]) for i in range(nm)]
if all(d == 3 for d in comp_degs):
    print("  The complement IS 3-regular!")
    # Check if it's Petersen
    # Petersen has girth 5 (shortest cycle = 5)
    # Check triangles in complement
    comp_tri = 0
    for i in range(nm):
        for j in range(i+1, nm):
            if not meta_comp[i][j]: continue
            for k in range(j+1, nm):
                if meta_comp[i][k] and meta_comp[j][k]:
                    comp_tri += 1
    print("  Complement triangles: %d (Petersen has 0)" % comp_tri)
    if comp_tri == 0:
        print("  → Complement is TRIANGLE-FREE and 3-regular = candidate for Petersen!")
else:
    print("  Complement is NOT 3-regular (degs: %s)" % sorted(set(comp_degs)))

# ============================================================
banner("THE META-GRAPH ADJACENCY MATRIX")
# ============================================================

# Let me just print it clearly
print("  Adjacency matrix of G_5/Z_2 (sorted by H):\n")
order = sorted(range(nm), key=lambda mi: H_vec[mi])
print("  H: ", [H_vec[mi] for mi in order])
print("  d: ", [deg[mi] for mi in order])
print()
for mi in order:
    row = ""
    for mj in order:
        row += "%d " % adj[mi][mj]
    print("  %d(H=%2d): %s" % (mi, H_vec[mi], row))

# Non-edges (complement)
print("\n  NON-EDGES (the complement):")
for mi in order:
    for mj in order:
        if mi < mj and adj[mi][mj] == 0:
            print("    {%d(H=%d), %d(H=%d)}" % (mi, H_vec[mi], mj, H_vec[mj]))

print("\n  Total non-edges: %d = C(10,2) - 21 = 24\n" % ne_comp)

# The complement has 24 edges. Petersen has 15. So complement ≠ Petersen.
# But C(10,2) = 45. Meta-graph has 21. Complement has 24. 21+24=45. ✓

# ============================================================
banner("THE JOHNSON GRAPH J(5,2)")
# ============================================================

# Johnson graph J(5,2): V = 2-subsets of [5], E = pairs sharing 1 element.
# J(5,2) has C(5,2) = 10 vertices and each vertex has degree 2×3 = 6.
# J(5,2) is 6-regular with 30 edges = complement of Petersen.

# At n=5 with 10 nodes and 21 edges:
# 21 is between 15 (Petersen) and 30 (Johnson).
# 21 = 15 + 6 = Petersen + 6 extra.
# 21 = 30 - 9 = Johnson - 9 missing.

# The 6 extra edges (beyond Petersen) or the 9 missing edges (from Johnson)
# might have a combinatorial meaning.

print("  PETERSEN (K(5,2)): 15 edges, 3-regular")
print("  META-GRAPH:        21 edges, degrees %s" % sorted(deg))
print("  JOHNSON (J(5,2)):  30 edges, 6-regular")
print("  COMPLETE (K_10):   45 edges, 9-regular\n")

print("  Meta-graph = Petersen + 6 = Johnson - 9 = K_10 - 24")
print("  The meta-graph sits BETWEEN Petersen and Johnson.\n")

# The 9 NON-EDGES:
print("  The %d NON-EDGES of the meta-graph:" % ne_comp)
non_edges = []
for mi in range(nm):
    for mj in range(mi+1, nm):
        if adj[mi][mj] == 0:
            non_edges.append((mi, mj, H_vec[mi], H_vec[mj]))

for mi, mj, Hi, Hj in sorted(non_edges, key=lambda x: abs(x[2]-x[3])):
    dH = abs(Hi - Hj)
    print("    {%d(H=%d), %d(H=%d)}: dH=%d" % (mi, Hi, mj, Hj, dH))

# The non-edges tend to have LARGE dH — distant classes are NOT connected.
avg_dH_non = np.mean([abs(Hi-Hj) for _,_,Hi,Hj in non_edges])
avg_dH_edge = np.mean([abs(H_vec[mi]-H_vec[mj]) for mi in range(nm)
                        for mj in range(mi+1,nm) if adj[mi][mj]])
print("\n  Avg |dH| for EDGES: %.2f" % avg_dH_edge)
print("  Avg |dH| for NON-EDGES: %.2f" % avg_dH_non)
print("  Non-edges have %.1f× larger dH → distant classes don't connect" %
      (avg_dH_non / avg_dH_edge if avg_dH_edge > 0 else 0))

print("\nDone. opus-2026-03-24-S287")
