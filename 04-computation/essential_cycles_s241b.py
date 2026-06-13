#!/usr/bin/env python3
"""
essential_cycles_s241b.py -- opus-2026-03-23-S241b

WHAT ARE THE ESSENTIAL 1-CYCLES?

At n=5, beta_1 = 2 means two independent cycles that are not boundaries.
Let's find them explicitly.

A 1-cycle is a closed walk in the graph that is not the boundary of
any set of triangles. Equivalently, a cycle in the graph that doesn't
bound a "disk" of triangles.

We can find generators of H_1 by computing the kernel of ∂_1 modulo
the image of ∂_2. The kernel of ∂_1 is all cycle vectors; the image
of ∂_2 is all boundary vectors (sums of triangle boundaries).

The quotient H_1 = ker(∂_1) / im(∂_2) has dimension beta_1.

Author: opus-2026-03-23-S241b
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

def is_sc(A, n):
    return canon(A,n) == canon(complement(A,n),n)

# Build for n=5
n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
cmap = {}; sizes = {}; reps = {}; tc = {}
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
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c], 'SC': is_sc(A, n)}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

F_mat = np.zeros((nc, nc), dtype=int)
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m): F_mat[ci][tc[bits^(1<<k)]] += 1

adj = [[False]*nm for _ in range(nm)]
for ci in range(nc):
    mi = mid_map[ci]
    for cj in range(nc):
        if F_mat[ci][cj] > 0 and mid_map[cj] != mi:
            adj[mi][mid_map[cj]] = True
            adj[mid_map[cj]][mi] = True

# H and SC labels
H_label = {mi: info[merged[mi][0]]['H'] for mi in range(nm)}
SC_label = {mi: info[merged[mi][0]]['SC'] for mi in range(nm)}

# ============================================================
banner("ESSENTIAL 1-CYCLES AT n=5")
# ============================================================

# Build edges and triangles
edges = [(i,j) for i in range(nm) for j in range(i+1,nm) if adj[i][j]]
triangles = []
for i in range(nm):
    for j in range(i+1,nm):
        if not adj[i][j]: continue
        for k in range(j+1,nm):
            if adj[i][k] and adj[j][k]:
                triangles.append((i,j,k))

ne = len(edges)
nt = len(triangles)
edge_idx = {e: i for i, e in enumerate(edges)}

print("  V=%d, E=%d, T=%d" % (nm, ne, nt))

# Build ∂_1: C_1 → C_0 (edges → vertices)
d1 = np.zeros((nm, ne), dtype=float)
for j, (a, b) in enumerate(edges):
    d1[a, j] = 1
    d1[b, j] = -1

# Build ∂_2: C_2 → C_1 (triangles → edges)
d2 = np.zeros((ne, nt), dtype=float)
for j, (a, b, c) in enumerate(triangles):
    # Boundary of (a,b,c) = (b,c) - (a,c) + (a,b)
    e_ab = edge_idx[(a, b)]
    e_ac = edge_idx[(a, c)]
    e_bc = edge_idx[(b, c)]
    d2[e_ab, j] = 1
    d2[e_ac, j] = -1
    d2[e_bc, j] = 1

# ker(∂_1) = null space of d1
# im(∂_2) = column space of d2
# H_1 = ker(∂_1) / im(∂_2)

# Find ker(∂_1): using SVD
U1, s1, V1h = np.linalg.svd(d1)
rank1 = np.sum(s1 > 1e-10)
nullity1 = ne - rank1
print("  rank(∂_1) = %d, nullity(∂_1) = %d" % (rank1, nullity1))

# Null space of d1 = last (ne - rank1) rows of V1h
null_d1 = V1h[rank1:, :]  # shape (nullity1, ne)

# Find im(∂_2): column space of d2
rank2 = np.linalg.matrix_rank(d2)
print("  rank(∂_2) = %d" % rank2)
print("  β_1 = nullity(∂_1) - rank(∂_2) = %d - %d = %d" % (nullity1, rank2, nullity1 - rank2))

# Find representatives of H_1 generators
# Project the null space basis of d1 modulo the column space of d2
# Combine: ker(∂_1) has basis null_d1[i], and we need to find cycles
# not in im(∂_2).

# Alternative: find a cycle basis of the graph
# A cycle basis has exactly β_1 = E - V + 1 = ne - nm + 1 = 21 - 10 + 1 = 12 fundamental cycles.
# Wait: for a connected graph, the cycle space has dimension E - V + 1 = nullity of ∂_1.
# That's 12 here. But β_1 = 2 because im(∂_2) has rank 10, so 12 - 10 = 2.

# So 10 of the 12 fundamental cycles are boundaries of triangles,
# and 2 are the essential (non-bounding) cycles.

# Let me find an explicit spanning tree and fundamental cycles.
# BFS from vertex 0
from collections import deque
tree_edges = set()
visited = {0}
queue = deque([0])
parent = {0: None}
while queue:
    v = queue.popleft()
    for j in range(nm):
        if adj[v][j] and j not in visited:
            visited.add(j)
            parent[j] = v
            tree_edges.add((min(v,j), max(v,j)))
            queue.append(j)

non_tree = [e for e in edges if e not in tree_edges]
print("\n  Spanning tree has %d edges" % len(tree_edges))
print("  Non-tree edges: %d (= fundamental cycles)" % len(non_tree))

# Each non-tree edge (a,b) creates a fundamental cycle:
# the unique tree path from a to b, plus the edge (a,b).
def tree_path(a, b):
    """Find path from a to b in the spanning tree."""
    # BFS in tree
    prev = {a: None}
    q = deque([a])
    while q:
        v = q.popleft()
        if v == b:
            path = []
            while v is not None:
                path.append(v)
                v = prev[v]
            return path[::-1]
        for w in range(nm):
            if w not in prev and (min(v,w), max(v,w)) in tree_edges:
                prev[w] = v
                q.append(w)
    return None

# Build fundamental cycles
fund_cycles = []
for e in non_tree:
    a, b = e
    path = tree_path(a, b)
    if path:
        cycle = list(path) + [path[0]]  # close the cycle
        fund_cycles.append(cycle)

print("\n  FUNDAMENTAL CYCLES:")
for i, cycle in enumerate(fund_cycles):
    H_vals = [H_label[v] for v in cycle[:-1]]
    SC_vals = ['SC' if SC_label[v] else 'NS' for v in cycle[:-1]]
    length = len(cycle) - 1
    print("    C%d (length %d): %s" % (i, length, ' → '.join('%d(H=%d,%s)' % (v, H_label[v], 'SC' if SC_label[v] else 'NS') for v in cycle[:-1])))

# Now: which of these are boundaries? A cycle is a boundary if it's
# a sum of triangle boundaries. Encode each cycle as a vector in Z^E.
cycle_vecs = []
for cycle in fund_cycles:
    vec = np.zeros(ne)
    for k in range(len(cycle)-1):
        a, b = cycle[k], cycle[k+1]
        e = (min(a,b), max(a,b))
        idx = edge_idx[e]
        # Orientation: if a < b, coefficient +1; else -1
        if a < b:
            vec[idx] += 1
        else:
            vec[idx] -= 1
    cycle_vecs.append(vec)

cycle_matrix = np.array(cycle_vecs)  # shape (12, 21)

# Check which cycles are in im(∂_2)
# A cycle c is a boundary iff c = ∂_2 * x for some x
# Solve d2 * x = c for each cycle
# If solvable, it's a boundary; otherwise, it's essential.

print("\n  TESTING WHICH CYCLES ARE BOUNDARIES:")
essential_idx = []
for i, vec in enumerate(cycle_vecs):
    # Solve d2 * x = vec using least squares
    x, res, rank_d2, sv = np.linalg.lstsq(d2, vec, rcond=None)
    residual = np.linalg.norm(d2 @ x - vec)
    is_boundary = residual < 1e-8
    if not is_boundary:
        essential_idx.append(i)
    status = "BOUNDARY" if is_boundary else "ESSENTIAL"
    print("    C%d: %s (residual=%.2e)" % (i, status, residual))

print("\n  ESSENTIAL CYCLES (not boundaries): %d" % len(essential_idx))
for i in essential_idx:
    cycle = fund_cycles[i]
    print("    C%d: %s" % (i, ' → '.join('%d(H=%d)' % (v, H_label[v]) for v in cycle[:-1])))

    # What non-tree edge creates this cycle?
    nte = non_tree[i]
    print("      Non-tree edge: %s (H=%d,%d)" % (str(nte), H_label[nte[0]], H_label[nte[1]]))

# ============================================================
banner("CYCLE ANALYSIS")
# ============================================================

print("  The %d essential cycles generate H_1(Δ(G_5/Z_2)).\n" % len(essential_idx))

# What H-values do essential cycles traverse?
for i in essential_idx:
    cycle = fund_cycles[i]
    H_vals = [H_label[v] for v in cycle[:-1]]
    print("  Cycle C%d:" % i)
    print("    Vertices: %s" % [v for v in cycle[:-1]])
    print("    H values: %s" % H_vals)
    print("    H range: [%d, %d]" % (min(H_vals), max(H_vals)))
    print("    SC nodes: %d, NS nodes: %d" %
          (sum(1 for v in cycle[:-1] if SC_label[v]),
           sum(1 for v in cycle[:-1] if not SC_label[v])))
    dH_vals = [abs(H_label[cycle[k]] - H_label[cycle[k+1]]) for k in range(len(cycle)-1)]
    print("    dH along cycle: %s" % dH_vals)
    print()

# ============================================================
banner("THE TOPOLOGICAL OBSTRUCTION")
# ============================================================

print("""
  The essential 1-cycles represent TOPOLOGICAL OBSTRUCTIONS in
  tournament space. They are loops in the meta-graph that cannot
  be "filled in" by any collection of triangles.

  This means: there exist paths through tournament iso classes
  that form IRREDUCIBLE LOOPS — sequences of single-arc-reversals
  that return to the starting class but whose "interior" cannot
  be tiled by triangles (3-cliques in the meta-graph).

  At n=5, there are EXACTLY 2 independent such obstructions.
  At n=6, there are 15 — a dramatic increase.

  These essential cycles might correspond to:
  - Score sequence barriers (tournaments with same H but unreachable
    via triangle-filling paths)
  - SC/NS boundary phenomena (cycles crossing the bipartition)
  - H-gradient reversals (cycles that go up-down-up-down)

  The 2-cycles (β_2 = 7 at n=6) are even more mysterious:
  they represent "hollow spheres" in the clique complex — 2D surfaces
  that cannot be filled by any 3D body (collection of K_4 subgraphs).
""")

print("Done. opus-2026-03-23-S241b")
