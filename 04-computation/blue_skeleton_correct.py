#!/usr/bin/env python3
"""
CORRECT blue line skeleton using tiling flip (complement non-backbone bits).

Tiling model:
- Backbone = standard path 0->1->2->...->n-1 (A[i][i+1] = 1 always)
- Tiling bits = A[i][j] for j-i >= 2 (C(n-1,2) bits)
- Flip = complement tiling bits, keep backbone
- GS = grid-symmetric: tiling bits invariant under staircase transpose

The blue line skeleton = GS flip graph connecting SC tournament classes.
At odd n, ALL GS flip pairs cross classes (the "tippy" property).

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    """Build tournament from tiling bits (backbone = 0->1->...->n-1)."""
    A = [[0]*n for _ in range(n)]
    # Set backbone: A[i][i+1] = 1
    for i in range(n-1):
        A[i][i+1] = 1

    # Set non-backbone from tiling bits
    idx = 0
    for i in range(n):
        for j in range(i+2, n):  # j-i >= 2
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    """Number of non-backbone edges."""
    return n*(n-1)//2 - (n-1)  # = C(n,2) - (n-1) = C(n-1,2)

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

def tiling_transpose_pairs(n):
    """Compute transpose-paired tiling bit indices for GS check.

    Non-backbone edges (i,j) with j-i >= 2.
    Transpose: (i,j) -> (n-1-j, n-1-i).
    """
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
        if tj - ti < 2:
            # Transpose lands on a backbone edge — this is a FIXED tile
            # Actually transpose of non-backbone should be non-backbone for the staircase
            # Wait: (i,j) -> (n-1-j, n-1-i). If n-1-j and n-1-i differ by 1, it's backbone.
            # j-i >= 2, n-1-i - (n-1-j) = j-i >= 2. So it's ALWAYS non-backbone!
            pass
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)

    return pairs, fixed

def is_gs_tiling(tiling_bits, n, pairs, fixed):
    for idx1, idx2 in pairs:
        if ((tiling_bits >> idx1) & 1) != ((tiling_bits >> idx2) & 1):
            return False
    return True

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 5
m = num_tiling_bits(n)
print(f"BLUE LINE SKELETON at n={n}")
print(f"  Non-backbone edges: {m}")
print(f"  Total tilings: {2**m}")

pairs, fixed = tiling_transpose_pairs(n)
print(f"  Transpose pairs: {pairs}")
print(f"  Fixed tiles: {fixed}")
print(f"  GS degrees of freedom: {len(pairs) + len(fixed)}")
print(f"  Expected GS count: {2**(len(pairs)+len(fixed))}")

# Find all GS tilings
gs_tilings = [bits for bits in range(2**m) if is_gs_tiling(bits, n, pairs, fixed)]
print(f"  Actual GS count: {len(gs_tilings)}")

# Build iso class database from tilings
canon_db = {}
class_list = []

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({
            'canon': c, 'rep': A, 'tilings': set(), 'sc': is_SC(A, n),
            'scores': score_seq(A, n), 'gs_tilings': set()
        })
    class_list[canon_db[c]]['tilings'].add(bits)

for bits in gs_tilings:
    A = tournament_from_tiling(n, bits)
    c_idx = canon_db[canonical(A, n)]
    class_list[c_idx]['gs_tilings'].add(bits)

num_classes = len(class_list)
num_sc = sum(1 for c in class_list if c['sc'])
print(f"\n  {num_classes} iso classes, {num_sc} SC")

# GS tilings per class
print(f"\n  Class details:")
for i, c in enumerate(class_list):
    tag = "SC" if c['sc'] else "NSC"
    gs_count = len(c['gs_tilings'])
    if gs_count > 0 or not c['sc']:
        print(f"    Class {i:2d} ({tag}, scores={c['scores']}): {gs_count} GS, {len(c['tilings'])} total")

# Build GS flip graph
print(f"\n{'='*60}")
print(f"GS FLIP GRAPH (BLUE LINE SKELETON)")
print(f"{'='*60}")

same_class = 0
cross_class = 0
gs_edges = defaultdict(int)

for bits in gs_tilings:
    A_from = tournament_from_tiling(n, bits)
    c_from = canon_db[canonical(A_from, n)]

    flipped = flip_tiling(bits, m)
    A_to = tournament_from_tiling(n, flipped)
    c_to = canon_db[canonical(A_to, n)]

    if c_from == c_to:
        same_class += 1
    else:
        cross_class += 1
        edge = (min(c_from, c_to), max(c_from, c_to))
        gs_edges[edge] += 1

print(f"  Same-class flips: {same_class}")
print(f"  Cross-class flips: {cross_class}")
print(f"  Fraction cross-class: {cross_class/(same_class+cross_class):.2%}")

print(f"\n  Blue line edges:")
for (i, j), count in sorted(gs_edges.items()):
    print(f"    {i} ({class_list[i]['scores']}) <-[{count}]-> {j} ({class_list[j]['scores']})")

# Now: full flip scatter matrix for ALL tilings (not just GS)
print(f"\n{'='*60}")
print(f"FULL TILING FLIP SCATTER (where do NSC tilings go?)")
print(f"{'='*60}")

nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
sc_indices = [i for i, c in enumerate(class_list) if c['sc']]

for i in nsc_indices:
    flip_targets = defaultdict(int)
    for bits in class_list[i]['tilings']:
        flipped = flip_tiling(bits, m)
        A_to = tournament_from_tiling(n, flipped)
        c_to = canon_db[canonical(A_to, n)]
        flip_targets[c_to] += 1

    A_op = converse(class_list[i]['rep'], n)
    partner = canon_db[canonical(A_op, n)]

    sc_targets = {k: v for k, v in flip_targets.items() if class_list[k]['sc']}
    nsc_targets = {k: v for k, v in flip_targets.items() if not class_list[k]['sc']}

    print(f"\n  NSC class {i} (scores={class_list[i]['scores']}, partner={partner}):")
    print(f"    -> SC targets: {dict(sc_targets)}")
    print(f"    -> NSC targets: {dict(nsc_targets)}")
    total_to_sc = sum(sc_targets.values())
    total_to_nsc = sum(nsc_targets.values())
    print(f"    SC: {total_to_sc}, NSC: {total_to_nsc} out of {len(class_list[i]['tilings'])}")

print("\nDONE")
