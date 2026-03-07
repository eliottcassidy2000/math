#!/usr/bin/env python3
"""
GS Cube Geometry: the GS subspace as a hypercube quotient.

GS tilings form a k-dimensional cube {0,1}^k (free bits).
The flip map is the antipodal map on this cube.
Isomorphism classes partition the cube into orbits.

This gives a natural metric on the skeleton: the Hamming distance
between GS tilings of different classes = geodesic distance on the cube.

Questions:
1. What is the diameter of the GS cube quotient?
2. Is the quotient metric tree-like?
3. How does the fiber size (#GS tilings per class) vary across the cube?
4. The GS cube as a COVER of the skeleton

Also: explore the "folding" structure.
- At odd n, flip maps each vertex to its antipodal point
- No fixed points (since flip changes t3 parity)
- The quotient {0,1}^k / flip is a real projective space RP^{k-1}
- The skeleton embeds in this quotient projective space

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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
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

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 5
m = num_tiling_bits(n)
print(f"GS CUBE GEOMETRY at n={n}")
print(f"{'='*60}")

pairs, fixed = tiling_transpose_pairs(n)
gs_dof = len(pairs) + len(fixed)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

print(f"  GS DOF (k): {gs_dof}")
print(f"  GS cube: {2**gs_dof} vertices = {len(gs_tilings)} GS tilings")

# Build class database
canon_db = {}
class_list = []
bits_to_class = {}
gs_to_freebits = {}

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({
            'rep': A, 'sc': is_SC(A, n),
            'scores': score_seq(A, n),
            't3': count_3cycles(A, n),
            'gs_tilings': set()
        })
    bits_to_class[bits] = canon_db[c]

# Map GS tilings to free-bit coordinates
for idx, bits in enumerate(gs_tilings):
    class_idx = bits_to_class[bits]
    class_list[class_idx]['gs_tilings'].add(bits)
    gs_to_freebits[bits] = idx  # The free_val index IS the coordinate

# The GS cube has vertices 0,...,2^k-1
# Two vertices are adjacent in the cube iff they differ in exactly 1 free bit
# The iso-class assignment colors the cube vertices

print(f"\n  GS cube vertex coloring (by iso class):")
for idx in range(2**gs_dof):
    bits = gs_tilings[idx]
    cls = bits_to_class[bits]
    t3 = class_list[cls]['t3']
    free_str = format(idx, f'0{gs_dof}b')
    if class_list[cls]['sc']:
        print(f"    free={free_str} -> class {cls} (SC, t3={t3})")

# Hamming distance between free-bit coordinates of different classes
print(f"\n  Inter-class Hamming distances (on GS cube):")
sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
for i, ci in enumerate(sc_indices):
    for j, cj in enumerate(sc_indices):
        if j <= i:
            continue
        min_dist = gs_dof + 1
        for bi in class_list[ci]['gs_tilings']:
            for bj in class_list[cj]['gs_tilings']:
                fi = gs_to_freebits[bi]
                fj = gs_to_freebits[bj]
                d = bin(fi ^ fj).count('1')
                if d < min_dist:
                    min_dist = d
        if min_dist <= 2:
            print(f"    Class {ci}(t3={class_list[ci]['t3']}) <-> "
                  f"Class {cj}(t3={class_list[cj]['t3']}): min Hamming dist = {min_dist}")

# GS cube with flip: check that flip = antipodal map on free bits
print(f"\n  Flip = antipodal on GS cube?")
for idx in range(min(8, 2**gs_dof)):
    bits = gs_tilings[idx]
    flipped = flip_tiling(bits, m)
    flip_idx = gs_tilings.index(flipped) if flipped in gs_tilings else -1
    expected_antipodal = (2**gs_dof - 1) ^ idx  # flip all free bits
    print(f"    free={idx:0{gs_dof}b} -> flip free={flip_idx:0{gs_dof}b}, "
          f"expected antipodal={expected_antipodal:0{gs_dof}b}, "
          f"match={'YES' if flip_idx == expected_antipodal else 'NO'}")

# The quotient: GS cube / flip
# At odd n, flip has no fixed points (bipartite skeleton)
# Quotient has 2^(k-1) elements
# Each element is a pair {x, antipodal(x)}
# This is the k-dimensional CROSS-POLYTOPE quotient = RP^{k-1}
print(f"\n  GS cube quotient:")
print(f"    2^k = {2**gs_dof}, quotient size = {2**gs_dof // 2}")
print(f"    Number of classes in quotient = number of skeleton edges + self-loops")

# Visualize the GS cube as a 4D hypercube projected to 2D
# The free bits are: pairs[0] (bit 0), pairs[1] (bit 1), fixed[0] (bit 2), fixed[1] (bit 3)
print(f"\n  Free bit assignment:")
for k, (idx1, idx2) in enumerate(pairs):
    edges_list = []
    for i2 in range(n):
        for j2 in range(i2+2, n):
            edges_list.append((i2, j2))
    print(f"    Bit {k}: pair ({edges_list[idx1]}, {edges_list[idx2]})")
for k, fidx in enumerate(fixed):
    edges_list = []
    for i2 in range(n):
        for j2 in range(i2+2, n):
            edges_list.append((i2, j2))
    print(f"    Bit {len(pairs)+k}: fixed edge {edges_list[fidx]}")

# How many orbits (iso classes) per Hamming weight level?
print(f"\n  Orbits per Hamming weight of free bits:")
for hw in range(gs_dof + 1):
    classes_at_hw = set()
    for idx in range(2**gs_dof):
        if bin(idx).count('1') == hw:
            bits = gs_tilings[idx]
            classes_at_hw.add(bits_to_class[bits])
    sc_at_hw = [c for c in classes_at_hw if class_list[c]['sc']]
    nsc_at_hw = [c for c in classes_at_hw if not class_list[c]['sc']]
    print(f"    HW={hw}: {len(sc_at_hw)} SC classes, {len(nsc_at_hw)} NSC classes "
          f"(total {len(classes_at_hw)}) out of C({gs_dof},{hw})={len([idx for idx in range(2**gs_dof) if bin(idx).count('1')==hw])} vertices")
