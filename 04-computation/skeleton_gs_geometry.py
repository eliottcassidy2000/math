#!/usr/bin/env python3
"""
GS (Grid-Symmetric) flip graph = the "blue line skeleton".

GS tilings are tilings invariant under the grid transpose.
At odd n, ALL GS flip pairs cross classes (100% cross-class = "tippy").
The skeleton connects SC tournament classes.

Then: how do NON-SC classes relate to this skeleton?
  - Each tiling T in an NSC class has flip(T) in the partner class
  - But some tilings of NSC class map CLOSER to certain SC classes

The user says: "most paired iso classes connect via black line to only one side"
This means: NSC class i's tilings, when projected onto the SC skeleton,
tend to cluster near one SC class rather than being evenly spread.

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

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

def get_tile_indices(n):
    """Get the staircase tile indices for the tiling representation."""
    tiles = []
    for i in range(n):
        for j in range(i+1, n):
            tiles.append((i, j))
    return tiles

def transpose_index(i, j, n):
    """Grid transpose: (i,j) -> (n-1-j, n-1-i) in the staircase."""
    return (n-1-j, n-1-i)

def is_gs_tiling(bits, n):
    """Check if a tiling is grid-symmetric."""
    tiles = get_tile_indices(n)
    m = len(tiles)
    tile_to_idx = {t: idx for idx, t in enumerate(tiles)}

    for idx, (i, j) in enumerate(tiles):
        ti, tj = transpose_index(i, j, n)
        # Make sure (ti, tj) is in upper-triangular form
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) not in tile_to_idx:
            continue
        tidx = tile_to_idx[(ti, tj)]
        # GS: bit at idx must equal bit at tidx
        if ((bits >> idx) & 1) != ((bits >> tidx) & 1):
            return False
    return True

n = 5
m = n*(n-1)//2
print(f"GS SKELETON GEOMETRY at n={n}")
print(f"{'='*60}")

# Find all GS tilings
gs_tilings = [bits for bits in range(2**m) if is_gs_tiling(bits, n)]
print(f"  {len(gs_tilings)} GS tilings out of {2**m} total")

# Build iso class database
canon_db = {}
class_list = []
for bits in range(2**m):
    A = tournament_from_bits(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({'canon': c, 'rep': A, 'tilings': set(), 'sc': is_SC(A, n),
                          'scores': score_seq(A, n), 'gs_tilings': set()})
    class_list[canon_db[c]]['tilings'].add(bits)

for bits in gs_tilings:
    A = tournament_from_bits(n, bits)
    c_idx = canon_db[canonical(A, n)]
    class_list[c_idx]['gs_tilings'].add(bits)

num_classes = len(class_list)
num_sc = sum(1 for c in class_list if c['sc'])
print(f"  {num_classes} iso classes, {num_sc} SC")

# GS tilings per class
print(f"\n  GS tilings per class:")
for i, c in enumerate(class_list):
    tag = "SC" if c['sc'] else "NSC"
    print(f"    Class {i:2d} ({tag}, scores={c['scores']}): {len(c['gs_tilings'])} GS, {len(c['tilings'])} total")

# Build GS flip graph (blue line skeleton)
print(f"\n{'='*60}")
print(f"BLUE LINE SKELETON (GS flip graph)")
print(f"{'='*60}")

gs_edges = defaultdict(int)  # (class_from, class_to) -> count
for bits in gs_tilings:
    A = tournament_from_bits(n, bits)
    c_from = canon_db[canonical(A, n)]

    flipped = bits ^ ((1 << m) - 1)
    A_flip = tournament_from_bits(n, flipped)
    c_to = canon_db[canonical(A_flip, n)]

    if c_from != c_to:
        edge = (min(c_from, c_to), max(c_from, c_to))
        gs_edges[edge] += 1

print(f"  GS flip edges (blue lines):")
for (i, j), count in sorted(gs_edges.items()):
    print(f"    {i} ({class_list[i]['scores']}) <-> {j} ({class_list[j]['scores']}): weight {count}")

# Now: how do NSC classes connect to the skeleton?
# For each NSC class, find which SC class its tilings are closest to
# "Closest" in Hamming distance
print(f"\n{'='*60}")
print(f"NSC CLASS CONNECTIONS TO SKELETON")
print(f"{'='*60}")

sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]

for i in nsc_indices:
    A_op = converse(class_list[i]['rep'], n)
    partner = canon_db[canonical(A_op, n)]

    # For each tiling in this NSC class, find min Hamming distance to each SC class
    min_dists = {j: float('inf') for j in sc_indices}
    avg_dists = {j: 0 for j in sc_indices}
    counts = {j: 0 for j in sc_indices}

    for bits_i in class_list[i]['tilings']:
        for j in sc_indices:
            for bits_j in class_list[j]['tilings']:
                dist = bin(bits_i ^ bits_j).count('1')
                min_dists[j] = min(min_dists[j], dist)
                avg_dists[j] += dist
                counts[j] += 1

    for j in sc_indices:
        if counts[j] > 0:
            avg_dists[j] /= counts[j]

    # Find closest SC class
    closest = min(sc_indices, key=lambda j: min_dists[j])
    second = sorted(sc_indices, key=lambda j: min_dists[j])[1]

    print(f"\n  NSC class {i} (scores={class_list[i]['scores']}, partner={partner}):")
    for j in sorted(sc_indices, key=lambda j: min_dists[j]):
        marker = " <-- CLOSEST" if j == closest else ""
        print(f"    -> SC {j:2d} (scores={class_list[j]['scores']}): min_dist={min_dists[j]}, avg={avg_dists[j]:.1f}{marker}")

    # Flip connections: where do NSC tilings go?
    flip_targets = defaultdict(int)
    for bits_i in class_list[i]['tilings']:
        flipped = bits_i ^ ((1 << m) - 1)
        A_flip = tournament_from_bits(n, flipped)
        target = canon_db[canonical(A_flip, n)]
        flip_targets[target] += 1

    print(f"    Flip targets: {dict(flip_targets)}")

print("\nDONE")
