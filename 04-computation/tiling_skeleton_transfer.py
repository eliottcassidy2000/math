#!/usr/bin/env python3
"""
The Tiling Skeleton: connecting tilings via single-tile flips,
colored by isomorphism class, annotated with transfer matrix data.

The "skeleton of blue pairs" refers to the graph structure where:
- Vertices = tilings in {0,1}^m (tournaments containing the base path)
- Edges = pairs differing by a single tile flip (Hamming distance 1)
- Colors = isomorphism classes of the resulting tournament

THM-030 gives us: M[a,b] = M[b,a] for every tournament.
What does this tell us about how M changes as we flip tiles?

Key questions:
1. How does M[a,b] change when we flip a single tile?
2. Is there a pattern in the "delta M" for blue pairs?
3. How does the skeleton connect isomorphism classes?
4. What is the role of grid-symmetry in the skeleton?
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tiling_to_tournament(bits, n):
    """Convert tiling bits to tournament adjacency matrix.
    Base path: n-1 -> n-2 -> ... -> 0 (vertices labeled 0..n-1)
    Tiles: pairs (a,b) with a >= b+2, listed in lexicographic order.
    Bit 0 = a->b (forward), bit 1 = b->a (backward).
    """
    A = [[0]*n for _ in range(n)]
    # Base path: i -> i-1 for i = n-1, n-2, ..., 1
    # So vertex i beats vertex i-1 (i.e., A[i][i-1] = 1)
    for i in range(1, n):
        A[i][i-1] = 1

    # Non-path tiles
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()

    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1  # backward: b -> a
        else:
            A[a][b] = 1  # forward: a -> b

    return A

def tournament_canonical(A):
    """Simple canonical form: sorted score sequence + hash of sorted adjacency."""
    n = len(A)
    scores = tuple(sorted([sum(row) for row in A]))
    # For more precise isomorphism, check all permutations
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def ham_path_count(A):
    n = len(A)
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

# =====================================================================
# ANALYSIS AT n=4
# =====================================================================
n = 4
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)
num_tilings = 2**m

print(f"n={n}: {m} tiles, {num_tilings} tilings")
print(f"Tiles: {tiles}")
print()

# Compute all tilings
tiling_data = {}
iso_classes = defaultdict(list)

for bits in range(num_tilings):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    canon = tournament_canonical(A)
    M = transfer_matrix(A)
    tiling_data[bits] = {
        'A': A, 'H': H, 'canon': canon, 'M': M,
        'bits': format(bits, f'0{m}b')
    }
    iso_classes[canon].append(bits)

n_classes = len(iso_classes)
print(f"Isomorphism classes: {n_classes}")

# Label classes
class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx
    members = iso_classes[canon]
    H = tiling_data[members[0]]['H']
    M = tiling_data[members[0]]['M']
    print(f"  Class {idx}: H={H}, |class|={len(members)}, tr(M)={np.trace(M)}, Sigma={M.sum()-np.trace(M)}")

# =====================================================================
# SKELETON: single-tile-flip connections
# =====================================================================
print()
print("=" * 70)
print("TILING SKELETON (single-tile flip connections)")
print("=" * 70)

# Count inter-class vs intra-class edges
intra_edges = 0
inter_edges = 0
inter_class_connections = defaultdict(int)

for bits in range(num_tilings):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:  # avoid double counting
            c1 = class_labels[tiling_data[bits]['canon']]
            c2 = class_labels[tiling_data[neighbor]['canon']]
            if c1 == c2:
                intra_edges += 1
            else:
                inter_edges += 1
                edge = (min(c1,c2), max(c1,c2))
                inter_class_connections[edge] += 1

print(f"\nTotal skeleton edges: {intra_edges + inter_edges}")
print(f"  Intra-class (same iso class): {intra_edges}")
print(f"  Inter-class (different iso class): {inter_edges}")

print(f"\nInter-class connections:")
for (c1, c2), count in sorted(inter_class_connections.items()):
    H1 = tiling_data[iso_classes[sorted(iso_classes.keys())[c1]][0]]['H']
    H2 = tiling_data[iso_classes[sorted(iso_classes.keys())[c2]][0]]['H']
    print(f"  Class {c1} (H={H1}) <-> Class {c2} (H={H2}): {count} edges")

# =====================================================================
# TRANSFER MATRIX CHANGES UNDER SINGLE TILE FLIPS
# =====================================================================
print()
print("=" * 70)
print("HOW DOES M CHANGE UNDER SINGLE TILE FLIPS?")
print("=" * 70)

# For each tile, compute delta_M = M(flipped) - M(original)
print("\nSample: tiling 0 (all forward arcs)")
bits = 0
A0 = tiling_data[0]['A']
M0 = tiling_data[0]['M']
H0 = tiling_data[0]['H']
print(f"  Tiling 0: bits={tiling_data[0]['bits']}, H={H0}")
print(f"  M = {M0.tolist()}")
print(f"  Class = {class_labels[tiling_data[0]['canon']]}")

for tile_idx in range(m):
    neighbor = 1 << tile_idx
    An = tiling_data[neighbor]['A']
    Mn = tiling_data[neighbor]['M']
    Hn = tiling_data[neighbor]['H']
    delta_M = Mn - M0
    a, b = tiles[tile_idx]
    cn = class_labels[tiling_data[neighbor]['canon']]
    print(f"\n  Flip tile {tile_idx} ({a},{b}): H {H0}->{Hn}, class->{cn}")
    print(f"    delta_M = {delta_M.tolist()}")
    print(f"    delta_M symmetric? {np.allclose(delta_M, delta_M.T)}")
    print(f"    delta_tr = {np.trace(delta_M)}, delta_Sigma = {delta_M.sum()-np.trace(delta_M)}")

# =====================================================================
# THE "BLUE PAIR" STRUCTURE
# =====================================================================
print()
print("=" * 70)
print("BLUE PAIR STRUCTURE: self-flip connections")
print("=" * 70)
print()
print("A 'blue pair' connects a tiling T to flip(T) (complement all bits).")
print("Self-flip: T and flip(T) are in the same isomorphism class.")
print()

flip_connections = {}
for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits  # complement all bits
    if flip_bits >= bits:  # avoid double counting
        c1 = class_labels[tiling_data[bits]['canon']]
        c2 = class_labels[tiling_data[flip_bits]['canon']]
        is_self_flip = (c1 == c2)
        flip_connections[(bits, flip_bits)] = {
            'self_flip': is_self_flip,
            'c1': c1, 'c2': c2,
            'H1': tiling_data[bits]['H'],
            'H2': tiling_data[flip_bits]['H']
        }

self_flip_pairs = [(b, f) for (b, f), d in flip_connections.items() if d['self_flip']]
non_self_flip = [(b, f) for (b, f), d in flip_connections.items() if not d['self_flip']]

print(f"Total flip pairs: {len(flip_connections)}")
print(f"Self-flip pairs (same class): {len(self_flip_pairs)}")
print(f"Non-self-flip pairs (different class): {len(non_self_flip)}")

for b, f in self_flip_pairs:
    d = flip_connections[(b, f)]
    print(f"  Self-flip: {tiling_data[b]['bits']} <-> {tiling_data[f]['bits']}, class={d['c1']}, H={d['H1']}")

# =====================================================================
# TRANSFER MATRIX UNDER FULL FLIP
# =====================================================================
print()
print("=" * 70)
print("TRANSFER MATRIX UNDER FULL FLIP (T -> T^op via non-path reversal)")
print("=" * 70)
print()
print("Flip reverses all non-path arcs, which is NOT the same as T^op")
print("(T^op reverses ALL arcs including path arcs).")
print()

for bits in range(min(num_tilings, 4)):
    flip_bits = ((1 << m) - 1) ^ bits
    M_orig = tiling_data[bits]['M']
    M_flip = tiling_data[flip_bits]['M']
    H_orig = tiling_data[bits]['H']
    H_flip = tiling_data[flip_bits]['H']

    print(f"Tiling {tiling_data[bits]['bits']} (H={H_orig}) vs flip {tiling_data[flip_bits]['bits']} (H={H_flip}):")
    print(f"  M_orig = {M_orig.tolist()}")
    print(f"  M_flip = {M_flip.tolist()}")
    print(f"  M_flip = M_orig? {np.allclose(M_orig, M_flip)}")
    print(f"  M_flip = M_orig^T? {np.allclose(M_flip, M_orig.T)}")
    print(f"  M_flip = -M_orig? {np.allclose(M_flip, -M_orig)}")

    # Check: is M_flip a permutation of M_orig?
    # Under relabeling v -> n-1-v?
    P = np.zeros((n,n))
    for i in range(n):
        P[i][n-1-i] = 1
    M_relabeled = P @ M_orig @ P.T
    print(f"  M_flip = P*M_orig*P^T (v->n-1-v)? {np.allclose(M_flip, M_relabeled)}")

print()
print("=" * 70)
print("DONE — Tiling skeleton analysis")
print("=" * 70)
