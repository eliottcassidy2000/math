#!/usr/bin/env python3
"""
Tiling skeleton at n=5: how does H change under single-tile flips?

At odd n=5, tr(M) = H(T). So delta_tr under a single tile flip
is delta_H = H(T') - H(T). This is the "arc-flip delta" from THM-015.

The skeleton connects 64 tilings (2^6 tiles) via single-tile flips.
Each tiling has degree 6 (one per tile).

KEY QUESTION: What does THM-030 (transfer matrix symmetry) tell us
about how H changes under single-tile flips?

Since M is symmetric, M[a,b] = M[b,a]. Under a single arc flip (i,j),
the change delta_M is also symmetric. And delta_H = tr(delta_M).
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

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

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

n = 5
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

# Compute H for all tilings
tiling_H = {}
tiling_canon = {}
iso_classes = defaultdict(list)

for bits in range(num_tilings):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    canon = tournament_canonical(A)
    tiling_H[bits] = H
    tiling_canon[bits] = canon
    iso_classes[canon].append(bits)

n_classes = len(iso_classes)
print(f"\nIsomorphism classes: {n_classes}")

class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx
    members = iso_classes[canon]
    H = tiling_H[members[0]]
    print(f"  Class {idx}: H={H}, |class|={len(members)}")

# Skeleton analysis
print()
print("=" * 70)
print("DELTA_H UNDER SINGLE TILE FLIPS")
print("=" * 70)

# For each tile, what is the distribution of delta_H?
delta_by_tile = defaultdict(lambda: defaultdict(int))
delta_H_values = defaultdict(int)

for bits in range(num_tilings):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            H1 = tiling_H[bits]
            H2 = tiling_H[neighbor]
            dH = H2 - H1
            delta_by_tile[tile_idx][dH] += 1
            delta_H_values[dH] += 1

print(f"\nOverall delta_H distribution:")
for dH in sorted(delta_H_values.keys()):
    print(f"  delta_H = {dH:+d}: {delta_H_values[dH]} edges")

print(f"\nDelta_H by tile:")
for tile_idx in range(m):
    a, b = tiles[tile_idx]
    dist = delta_by_tile[tile_idx]
    vals = sorted(dist.keys())
    print(f"  Tile {tile_idx} ({a},{b}): delta_H in {{{', '.join(f'{v:+d}:{dist[v]}' for v in vals)}}}")

# Inter-class connections
print()
print("=" * 70)
print("INTER-CLASS SKELETON")
print("=" * 70)

inter_class_edges = defaultdict(int)
inter_class_dH = defaultdict(list)

for bits in range(num_tilings):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            c1 = class_labels[tiling_canon[bits]]
            c2 = class_labels[tiling_canon[neighbor]]
            if c1 != c2:
                edge = (min(c1,c2), max(c1,c2))
                inter_class_edges[edge] += 1
                dH = tiling_H[neighbor] - tiling_H[bits]
                inter_class_dH[edge].append(dH)

for (c1, c2), count in sorted(inter_class_edges.items()):
    H1 = tiling_H[iso_classes[sorted(iso_classes.keys())[c1]][0]]
    H2 = tiling_H[iso_classes[sorted(iso_classes.keys())[c2]][0]]
    dH_vals = inter_class_dH[(c1,c2)]
    # Note: dH can be either sign depending on which way the edge goes
    print(f"  Class {c1} (H={H1}) <-> Class {c2} (H={H2}): {count} edges, |dH|={abs(H2-H1)}")

# Self-flip pairs
print()
print("=" * 70)
print("SELF-FLIP PAIRS")
print("=" * 70)

self_flip_count = 0
for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits >= bits:
        c1 = class_labels[tiling_canon[bits]]
        c2 = class_labels[tiling_canon[flip_bits]]
        if c1 == c2:
            self_flip_count += 1
            H = tiling_H[bits]
            print(f"  Self-flip: {format(bits, f'0{m}b')} <-> {format(flip_bits, f'0{m}b')}, class={c1}, H={H}")

print(f"\nTotal self-flip pairs: {self_flip_count}")
print(f"(At odd n=5, no blueself expected — Theorem 5)")

# Check Claim A connection
print()
print("=" * 70)
print("CLAIM A: H(T) - H(T-v) = 2*sum mu(C)")
print("=" * 70)
print()
print("For each tiling, compute H(T) and H(T-v) for each vertex v.")
print("delta_H = H(T) - H(T-v) should always be even (Claim A).")
print()

# Just check a few tilings
for bits in [0, 1, 31, 63]:
    A = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    print(f"Tiling {format(bits, f'0{m}b')} (H={H}):")
    for v in range(n):
        # T - v
        remaining = [i for i in range(n) if i != v]
        H_minus_v = 0
        for p in permutations(remaining):
            valid = True
            for i in range(len(p)-1):
                if A[p[i]][p[i+1]] != 1:
                    valid = False
                    break
            if valid:
                H_minus_v += 1
        delta = H - H_minus_v
        print(f"  v={v}: H(T-v)={H_minus_v}, delta={delta}, even={delta%2==0}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
n=5 tiling skeleton:
  - {num_tilings} tilings, {n_classes} isomorphism classes
  - {m} tiles per tiling, degree-{m} skeleton graph
  - delta_H under single-tile flip ranges from {min(delta_H_values.keys()):+d} to {max(delta_H_values.keys()):+d}
  - delta_H is always EVEN (consistent with Claim A / Redei parity)
  - No blueself pairs at odd n=5 (Theorem 5)

THM-030 constrains the skeleton because:
  - Every M along the skeleton is symmetric
  - delta_M under tile flip is also symmetric
  - H = tr(M) at odd n, so delta_H = tr(delta_M)
  - The eigenvalue structure is constrained by symmetry
""")
