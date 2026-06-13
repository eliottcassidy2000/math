#!/usr/bin/env python3
"""
IS the blue line skeleton = the quotient of the GS cube by antipodal + isomorphism?

The GS cube is {0,1}^k where k = gs_dof.
Flip = antipodal map (complement all bits).
Isomorphism identifies GS tilings in the same class.

The skeleton edges connect classes that are Hamming-distance-1 neighbors
in the GS cube. But wait: the skeleton edges connect classes whose GS tilings
are antipodal pairs. That's Hamming distance k (maximal!), not 1.

Actually, let me reconsider. The skeleton edges come from GS flip pairs:
T and flip(T) = antipodal(T). These are at Hamming distance k.
But the SKELETON EDGES are between CLASSES, not between cube vertices.
Two classes are adjacent in the skeleton iff some GS tiling in one class
is the flip of a GS tiling in the other class.

The skeleton is NOT the quotient cube graph. It's the flip-adjacency graph
on the quotient.

But there's another graph: the CUBE ADJACENCY on the quotient.
Two classes are cube-adjacent if some GS tiling in one class is
Hamming-distance-1 from a GS tiling in the other class.

Compare: skeleton (flip adjacency) vs cube adjacency.

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs, fixed, seen = [], [], set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen: continue
        ti, tj = n-1-j, n-1-i
        if ti > tj: ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx); seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx)); seen.add(idx); seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_dof = len(pairs) + len(fixed)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Build class database
canon_db = {}
class_list = []
bits_to_class = {}
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({'sc': is_SC(A, n), 't3': count_3cycles(A, n),
                          'scores': score_seq(A, n), 'gs_tilings': set()})
    bits_to_class[bits] = canon_db[c]

for bits in gs_tilings:
    class_list[bits_to_class[bits]]['gs_tilings'].add(bits)

sc_indices = [i for i, c in enumerate(class_list) if c['sc']]

# Map GS tilings to free-bit index
gs_to_free = {bits: idx for idx, bits in enumerate(gs_tilings)}

# Build CUBE adjacency graph on classes (Hamming-1 neighbors in free bits)
cube_adj = defaultdict(set)
for idx1 in range(2**gs_dof):
    c1 = bits_to_class[gs_tilings[idx1]]
    for bit in range(gs_dof):
        idx2 = idx1 ^ (1 << bit)
        c2 = bits_to_class[gs_tilings[idx2]]
        if c1 != c2:
            cube_adj[c1].add(c2)
            cube_adj[c2].add(c1)

# Build FLIP adjacency graph (the skeleton)
flip_adj = defaultdict(set)
for bits in gs_tilings:
    c_from = bits_to_class[bits]
    flipped = bits ^ ((1 << m) - 1)
    c_to = bits_to_class[flipped]
    if c_from != c_to:
        flip_adj[c_from].add(c_to)
        flip_adj[c_to].add(c_from)

print(f"GS CUBE ADJACENCY vs FLIP ADJACENCY at n={n}")
print(f"  k={gs_dof}, m={m}")
print()

print(f"  CUBE ADJACENCY (Hamming-1 neighbors on GS cube):")
for i in sorted(sc_indices):
    nbrs = sorted(cube_adj[i] & set(sc_indices))
    print(f"    Class {i} (t3={class_list[i]['t3']}): cube neighbors = {nbrs}")

print(f"\n  FLIP ADJACENCY (skeleton):")
for i in sorted(sc_indices):
    nbrs = sorted(flip_adj[i] & set(sc_indices))
    print(f"    Class {i} (t3={class_list[i]['t3']}): flip neighbors = {nbrs}")

# Compare
print(f"\n  COMPARISON:")
all_same = True
for i in sc_indices:
    cube_n = sorted(cube_adj[i] & set(sc_indices))
    flip_n = sorted(flip_adj[i] & set(sc_indices))
    match = cube_n == flip_n
    if not match:
        all_same = False
        print(f"    Class {i}: cube={cube_n}, flip={flip_n} -> DIFFERENT")
    else:
        print(f"    Class {i}: both={cube_n}")

print(f"\n  Cube adjacency = Flip adjacency: {all_same}")

# The cube adjacency has MORE edges than the flip adjacency
cube_edges = set()
flip_edges = set()
for i in sc_indices:
    for j in cube_adj[i]:
        if j in set(sc_indices) and i < j:
            cube_edges.add((i, j))
    for j in flip_adj[i]:
        if j in set(sc_indices) and i < j:
            flip_edges.add((i, j))

print(f"\n  Cube edges (SC-SC): {len(cube_edges)}")
print(f"  Flip edges (SC-SC): {len(flip_edges)}")
print(f"  Cube \\ Flip: {cube_edges - flip_edges}")
print(f"  Flip \\ Cube: {flip_edges - cube_edges}")

# The cube adjacency graph on classes: single-bit-flip on free coordinates
# This is the LOCAL structure (changing one edge of the tournament)
# The flip adjacency graph: antipodal map (changing ALL non-backbone edges)
# These are fundamentally different: one is local, the other is global

# The cube adjacency graph changes t3 by at most a few
# Let's check: for cube-adjacent classes, what's the t3 difference?
print(f"\n  t3 changes across cube edges:")
for (i, j) in sorted(cube_edges):
    dt3 = abs(class_list[i]['t3'] - class_list[j]['t3'])
    print(f"    {i}(t3={class_list[i]['t3']}) -- {j}(t3={class_list[j]['t3']}): |dt3|={dt3}")
