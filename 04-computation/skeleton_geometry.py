#!/usr/bin/env python3
"""
GEOMETRIC EXPLORATION: The blue line skeleton and its asymmetry.

The GS (grid-symmetric) cross-class flip graph is the "blue line skeleton":
- Vertices = self-converse (SC) tournament iso classes
- Edges = GS flip pairs between different classes
- At odd n, 100% of GS flip pairs cross classes (the skeleton is "tippy")
- Non-SC classes connect to this skeleton via "black lines"

Key observation from user: "most paired iso classes connect via black line
to only one side" — meaning the non-SC classes preferentially connect to
one SC class more than its flip partner.

This script explores:
1. The flip scatter matrix F[i][j] at n=5,7
2. How non-SC classes connect to the SC skeleton
3. The "tippiness" / asymmetry of this connection
4. Mobius strip structure: does the skeleton have non-orientable topology?

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

def flip_tiling(bits, m):
    """Flip = complement all bits."""
    return bits ^ ((1 << m) - 1)

def are_isomorphic_cached(A, B, n, cache_A=None):
    if score_seq(A, n) != score_seq(B, n):
        return False
    cA = cache_A if cache_A else canonical(A, n)
    return cA == canonical(B, n)

n = 5
m = n*(n-1)//2
print(f"SKELETON GEOMETRY at n={n}")
print(f"{'='*60}")

# Build iso class database
canon_db = {}
class_list = []
for bits in range(2**m):
    A = tournament_from_bits(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({'canon': c, 'rep': A, 'tilings': [], 'sc': is_SC(A, n),
                          'scores': score_seq(A, n)})
    class_list[canon_db[c]]['tilings'].append(bits)

num_classes = len(class_list)
num_sc = sum(1 for c in class_list if c['sc'])
print(f"  {num_classes} iso classes, {num_sc} self-converse")

# Build flip scatter matrix
F = [[0]*num_classes for _ in range(num_classes)]
for bits in range(2**m):
    A = tournament_from_bits(n, bits)
    c_from = canon_db[canonical(A, n)]

    flipped = flip_tiling(bits, m)
    A_flip = tournament_from_bits(n, flipped)
    c_to = canon_db[canonical(A_flip, n)]

    F[c_from][c_to] += 1

# Identify SC classes and their flip connections
print(f"\n  SC classes and their flip graph neighbors:")
sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]

for i in sc_indices:
    neighbors = [(j, F[i][j]) for j in range(num_classes) if F[i][j] > 0 and j != i]
    self_flip = F[i][i]
    sc_nbrs = [(j, w) for j, w in neighbors if class_list[j]['sc']]
    nsc_nbrs = [(j, w) for j, w in neighbors if not class_list[j]['sc']]
    print(f"  Class {i} (SC, scores={class_list[i]['scores']}, {len(class_list[i]['tilings'])} tilings):")
    print(f"    self-flip: {self_flip}, SC neighbors: {sc_nbrs}, NSC neighbors: {nsc_nbrs}")

# For non-SC classes, find their CONVERSE partner
print(f"\n  Non-SC (paired) classes:")
for i in nsc_indices:
    # Find converse class
    A_op = converse(class_list[i]['rep'], n)
    c_op = canon_db[canonical(A_op, n)]

    # How does this class connect to SC classes via flip?
    sc_connections = [(j, F[i][j]) for j in sc_indices if F[i][j] > 0]

    print(f"  Class {i} (NSC, scores={class_list[i]['scores']}, partner={c_op}):")
    print(f"    flip connections to SC: {sc_connections}")

    # Does the partner connect to the SAME SC classes?
    partner_sc = [(j, F[c_op][j]) for j in sc_indices if F[c_op][j] > 0]
    print(f"    partner's SC connections: {partner_sc}")

    # Asymmetry: does this class connect more to one side than another?
    total_to_sc = sum(w for _, w in sc_connections)
    total_partner_to_sc = sum(w for _, w in partner_sc)
    print(f"    total flip weight to SC: {total_to_sc} vs partner: {total_partner_to_sc}")

# Now analyze the "sidedness" — which SC class does each NSC class prefer?
print(f"\n{'='*60}")
print(f"ASYMMETRY ANALYSIS: Do NSC classes 'lean' toward one SC class?")
print(f"{'='*60}")

for i in nsc_indices:
    A_op = converse(class_list[i]['rep'], n)
    c_op = canon_db[canonical(A_op, n)]

    # Hamming distance from each tiling of class i to each SC class
    sc_distances = {j: [] for j in sc_indices}
    for bits in class_list[i]['tilings']:
        for j in sc_indices:
            for sc_bits in class_list[j]['tilings']:
                dist = bin(bits ^ sc_bits).count('1')
                sc_distances[j].append(dist)

    for j in sc_indices:
        if sc_distances[j]:
            avg_dist = sum(sc_distances[j]) / len(sc_distances[j])
            min_dist = min(sc_distances[j])
            # print(f"  Class {i} -> SC class {j}: avg_dist={avg_dist:.2f}, min_dist={min_dist}")

# Degree centrality in the skeleton
print(f"\n{'='*60}")
print(f"SKELETON DEGREE CENTRALITY")
print(f"{'='*60}")
for i in sc_indices:
    sc_degree = sum(F[i][j] for j in sc_indices if j != i)
    nsc_degree = sum(F[i][j] for j in nsc_indices)
    total = sum(F[i][j] for j in range(num_classes))
    print(f"  SC class {i} (scores={class_list[i]['scores']}): SC-edges={sc_degree}, NSC-edges={nsc_degree}, total={total}")

print("\nDONE")
