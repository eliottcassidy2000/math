#!/usr/bin/env python3
"""
petersen_alpha2_87.py — opus-2026-03-14-S87

THE PETERSEN GRAPH AND THE α₂=4 PATTERNS

The 10 partitions of {0,...,5} into two triples are in bijection with
the EDGES of K₆ (or equivalently, 2-element subsets of {0,...,5} via
the complement triple). The "both-cyclic" patterns form subgraphs.

Actually, the 10 partitions correspond to PERFECT MATCHINGS of K₆.
Wait no — 10 = C(6,3)/2 partitions into {3,3}.

The key observation: the 15 α₂=4 patterns each use 4 of the 10 partitions.
What graph do these form?

Also exploring: the connection to the Fano plane and Steiner systems.
"""

from itertools import combinations
from collections import Counter, defaultdict

# ── Setup: partitions of {0,...,5} into two triples ──

verts = list(range(6))
partitions = []
for triple in combinations(verts, 3):
    comp = tuple(v for v in verts if v not in triple)
    if triple < comp:
        partitions.append((triple, comp))

print("=" * 70)
print("THE 10 PARTITIONS AND THEIR STRUCTURE")
print("=" * 70)
print()

# Each partition corresponds to a perfect matching of K₆ by identifying
# each triple with its set of 3 arcs (3-cycle or transitive).
# Actually, each partition splits the 15 arcs into:
#   3 within A, 3 within B, 9 between A and B.

# The 10 partitions form vertices of a graph where two partitions
# are "adjacent" if they share NO common element in either triple.
# Actually that never happens since any two 3-element subsets of
# a 6-element set must overlap.

# Let's define: two partitions CONFLICT if they share at least 2 vertices
# in one of their triples (strong overlap).

# More naturally: two partitions are COMPATIBLE if one can find a
# perfect matching of K₆ containing both partition edges.

# The INTERSECTION PATTERN:
# Partition P_i = {A_i, B_i}. Define overlap between P_i and P_j as
# (|A_i ∩ A_j|, |A_i ∩ B_j|) (and symmetric variants).

print("Partition intersection patterns:")
for i in range(10):
    for j in range(i+1, 10):
        Ai, Bi = set(partitions[i][0]), set(partitions[i][1])
        Aj, Bj = set(partitions[j][0]), set(partitions[j][1])
        # The overlap type is determined by |A_i ∩ A_j|
        aa = len(Ai & Aj)
        ab = len(Ai & Bj)
        # aa + ab = 3 (since |A_i| = 3 and A_j ∪ B_j = everything)
        if j <= i + 3:  # just print a few
            pass
    break

# Classify partition pairs by overlap type
overlap_types = Counter()
for i in range(10):
    for j in range(i+1, 10):
        Ai = set(partitions[i][0])
        Aj = set(partitions[j][0])
        aa = len(Ai & Aj)
        overlap_types[aa] += 1

print(f"Overlap types (|A_i ∩ A_j|) distribution:")
for aa, count in sorted(overlap_types.items()):
    # aa can be 0, 1, 2, 3
    # If aa=0: A_i and A_j are complementary triples, so A_i = B_j
    # If aa=1: they share 1 vertex in one triple
    # If aa=2: they share 2 vertices
    # If aa=3: they're the same partition
    print(f"  |A_i ∩ A_j| = {aa}: {count} pairs")

# ══════════════════════════════════════════════════════════════════
# THE α₂=4 PATTERNS — STRUCTURE
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE 15 α₂=4 PATTERNS")
print("=" * 70)

# From previous computation: there are 15 distinct 4-element subsets
# that appear as α₂=4 patterns. Let's find them.

def all_tournaments_n6():
    edges = [(i,j) for i in range(6) for j in range(i+1,6)]
    m = len(edges)
    for bits in range(1 << m):
        adj = [[0]*6 for _ in range(6)]
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def is_cyclic_triple(adj, triple):
    a, b, c = triple
    return ((adj[a][b] and adj[b][c] and adj[c][a]) or
            (adj[a][c] and adj[c][b] and adj[b][a]))

alpha4_patterns = set()
for adj in all_tournaments_n6():
    pattern = frozenset(i for i, (A, B) in enumerate(partitions)
                        if is_cyclic_triple(adj, A) and is_cyclic_triple(adj, B))
    if len(pattern) == 4:
        alpha4_patterns.add(pattern)

print(f"\n{len(alpha4_patterns)} distinct α₂=4 patterns:")
for p in sorted(alpha4_patterns, key=lambda x: sorted(x)):
    indices = sorted(p)
    parts = [(partitions[i][0], partitions[i][1]) for i in indices]
    # What's the union of all triples involved?
    all_triples = set()
    for A, B in parts:
        all_triples.add(A)
        all_triples.add(B)
    print(f"  {indices}: triples = {[list(t) for t in sorted(all_triples)]}")

# ══════════════════════════════════════════════════════════════════
# THE COMPLEMENT GRAPH OF PARTITIONS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE PARTITION COMPLEMENT GRAPH")
print("=" * 70)
print()
print("Two partitions are COMPLEMENT-ADJACENT if they can be simultaneously")
print("both-cyclic in some tournament (compatible forcing).")
print()

# Build adjacency matrix: P_i ~ P_j if they appear together in some both-cyclic pattern
co_occurrence = [[0]*10 for _ in range(10)]
pattern_count = Counter()
for adj in all_tournaments_n6():
    both_cyclic = frozenset(i for i, (A, B) in enumerate(partitions)
                            if is_cyclic_triple(adj, A) and is_cyclic_triple(adj, B))
    pattern_count[both_cyclic] += 1
    for i in both_cyclic:
        for j in both_cyclic:
            if i != j:
                co_occurrence[i][j] += 1

print("Co-occurrence matrix (# tournaments where both P_i and P_j are both-cyclic):")
for i in range(10):
    row = [co_occurrence[i][j] for j in range(10)]
    print(f"  P_{i}: {row}")

# Build the "can coexist" graph
can_coexist = [[1 if co_occurrence[i][j] > 0 else 0 for j in range(10)] for i in range(10)]
edge_count = sum(can_coexist[i][j] for i in range(10) for j in range(i+1, 10))
print(f"\nCoexistence graph: {edge_count} edges out of C(10,2) = 45")
print("Degree sequence:", [sum(can_coexist[i][j] for j in range(10) if j != i) for i in range(10)])

# ══════════════════════════════════════════════════════════════════
# THE FORCING PROOF: 3 → 4
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE FORCING PROOF: WHY 3 BOTH-CYCLIC → 4")
print("=" * 70)
print()

# For each triple of partitions that CAN be simultaneously both-cyclic,
# identify which 4th partition is FORCED.
forced = defaultdict(set)
for adj in all_tournaments_n6():
    both_cyclic = frozenset(i for i, (A, B) in enumerate(partitions)
                            if is_cyclic_triple(adj, A) and is_cyclic_triple(adj, B))
    if len(both_cyclic) >= 3:
        for triple in combinations(sorted(both_cyclic), 3):
            extras = both_cyclic - set(triple)
            forced[frozenset(triple)].update(extras)

print(f"Triples that force a 4th (out of C(10,3) = 120):")
forcing_triples = {t: extras for t, extras in forced.items() if extras}
non_forcing = {t: extras for t, extras in forced.items() if not extras}

print(f"  Triples that ALWAYS force at least one extra: {len(forcing_triples)}")
print(f"  Triples that appear but force 0 extras: {len(non_forcing)}")

# Are ALL 3-subsets that can appear forcing?
# The key: we showed α₂=3 is impossible, so whenever 3 appear, 4 must.
# Let's check: do all feasible 3-subsets force exactly 1 extra?

for triple, extras in sorted(forcing_triples.items(), key=lambda x: sorted(x[0])):
    t = sorted(triple)
    e = sorted(extras)
    if len(t) <= 3:  # just first few
        print(f"  {t} → forces {e}")
    if len(e) == 1:
        pass  # exactly 1 forced
    else:
        print(f"  NOTE: {t} has {len(e)} possible extras (but always at least 1)")

# How many force exactly 1?
exactly_1 = sum(1 for t, e in forcing_triples.items() if len(e) == 1)
print(f"\n  Force exactly 1 extra: {exactly_1}")
print(f"  Force 2+ extras: {len(forcing_triples) - exactly_1}")

# ══════════════════════════════════════════════════════════════════
# THE PETERSEN GRAPH CONNECTION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PETERSEN GRAPH AND S₆ SYMMETRY")
print("=" * 70)
print()
print("The 10 partitions of {0,...,5} into {3,3} are the vertices of")
print("the Petersen graph's complement! Two partitions are Petersen-adjacent")
print("iff they have NO common triple (i.e., |A_i ∩ A_j| ∈ {1,2}).")
print()
print("Actually, the PETERSEN GRAPH has 10 vertices and 15 edges.")
print("It's the Kneser graph K(5,2) — but our 10 partitions are different.")
print()
print("Our 10 partitions form the vertices of the COMPLETE MULTIPARTITE")
print("graph K_{1,1,...,1} = K₁₀... let me think about this differently.")
print()

# The 10 partitions have a natural identification:
# Each partition {A, B} with A = {a,b,c} corresponds to the 3-element subset A
# (since B is determined). But we used canonical ordering A < B.
# So the 10 partitions biject with C(6,3)/2 = 10 unordered pairs {A, 6\A}.

# The PAIRING graph: partition i ~ partition j if they share a triple.
# Two partitions share a triple iff A_i = A_j or A_i = B_j.
# A_i = A_j means same partition. A_i = B_j means they're COMPLEMENTARY.
# Since each A has a unique complement B, each partition has exactly 1
# complementary partner. But since we use A < B, these are the SAME partition.
# Wait no, A and B are DIFFERENT subsets but the SAME partition.

# So no two DISTINCT partitions share a triple!
# The overlap is in PAIRS of vertices, not triples.

# The natural adjacency: Kneser-like
# Define P_i ~ P_j iff A_i ∩ A_j = ∅.
# But |A_i| = 3, |{0,...,5}| = 6, so A_i ∩ A_j = ∅ iff A_i = B_j.
# But that means the partitions are "opposite" — same partition!
# So in our canonical form, no two distinct partitions have A_i ∩ A_j = ∅.

# Actually wait: A_i and A_j are 3-subsets of {0,...,5}.
# A_i ∩ A_j can have 0, 1, 2, or 3 elements.
# If A_i ∩ A_j = ∅ then A_j = B_i, so partition j = partition i (same!)
# For DISTINCT partitions, |A_i ∩ A_j| ∈ {1, 2}.

# So the intersection types are:
# |A_i ∩ A_j| = 1: they share 1 vertex in each "A" triple
# |A_i ∩ A_j| = 2: they share 2 vertices in one "A" triple

# Count:
type_1 = 0
type_2 = 0
for i in range(10):
    for j in range(i+1, 10):
        aa = len(set(partitions[i][0]) & set(partitions[j][0]))
        if aa == 1:
            type_1 += 1
        elif aa == 2:
            type_2 += 1

print(f"Intersection types among 10 partitions:")
print(f"  |A_i ∩ A_j| = 1: {type_1} pairs")
print(f"  |A_i ∩ A_j| = 2: {type_2} pairs")
print(f"  Total: {type_1 + type_2} = C(10,2) = {10*9//2}")
print()

# The 15 α₂=4 patterns: are they the 15 MATCHINGS of the Petersen graph?
# The Petersen graph has 2000 perfect matchings... too many.
# But it has 15 edges. And we have 15 α₂=4 patterns.

# Check: do the 15 patterns correspond to the 15 edges of a specific graph on 10 vertices?
# A graph on 10 vertices with 15 edges = K₆ (if we identify vertices with partitions somehow).
# Actually C(6,2) = 15 and we have 15 patterns. Coincidence?

# Build the "co-4" graph: vertices = the 15 α₂=4 patterns
# Actually, let's see if the 15 patterns form the edges of the Petersen graph complement.

# The Petersen graph complement has 10 vertices and C(10,2)-15 = 30 edges.
# So Petersen complement has 30 edges, not 15.

# Let's just examine the structure of the 15 patterns directly.
print("The 15 α₂=4 patterns as 4-subsets of 10 partitions:")
patterns_list = sorted(alpha4_patterns, key=lambda x: sorted(x))
for k, p in enumerate(patterns_list):
    print(f"  Pattern {k}: {sorted(p)}")

# Check: is each partition in the same number of patterns?
partition_appearances = Counter()
for p in patterns_list:
    for i in p:
        partition_appearances[i] += 1
print(f"\nPartition appearance counts in α₂=4 patterns:")
for i in range(10):
    print(f"  P_{i}: appears in {partition_appearances[i]} patterns")

# Each in 6! So it's a 6-regular design.
# 15 patterns × 4 elements / 10 partitions = 6 appearances each. ✓

# Do patterns partition into parallel classes?
# 15 / 5 = 3 parallel classes of 5 if they partition the 10 partitions.
# But patterns have 4 elements each and there are 10 partitions, so
# 5 disjoint patterns would cover 5×4 = 20 slots, but only 10 partitions...
# Not a partition.

# Check: how many patterns are pairwise disjoint?
max_disjoint = 0
for i in range(15):
    for j in range(i+1, 15):
        if not (patterns_list[i] & patterns_list[j]):
            max_disjoint += 1
print(f"\nPairwise disjoint pattern pairs: {max_disjoint}")

# Intersection sizes between patterns
inter_sizes = Counter()
for i in range(15):
    for j in range(i+1, 15):
        inter_sizes[len(patterns_list[i] & patterns_list[j])] += 1
print(f"\nIntersection size distribution between α₂=4 patterns:")
for s, c in sorted(inter_sizes.items()):
    print(f"  |P_i ∩ P_j| = {s}: {c} pairs")

# ══════════════════════════════════════════════════════════════════
# SYNTHESIS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print()
print("The 15 α₂=4 patterns form a (10, 4, 6)-design:")
print("  10 points (partitions), 15 blocks of size 4, each point in 6 blocks.")
print("  This is a 2-design (BIBD) if every pair appears in λ blocks.")
print()

# Check: is it a 2-design?
pair_appearances = Counter()
for p in patterns_list:
    for i, j in combinations(sorted(p), 2):
        pair_appearances[(i,j)] += 1

lambda_vals = set(pair_appearances.values())
print(f"  Pair appearances λ: {lambda_vals}")
if len(lambda_vals) == 1:
    lam = lambda_vals.pop()
    print(f"  ✓ It's a 2-(10, 4, {lam}) design!")
    # Check: b×k×(k-1) = v×(v-1)×λ → 15×4×3 = 10×9×λ → 180 = 90λ → λ=2
    print(f"  Fisher's inequality: b ≥ v: 15 ≥ 10 ✓")
else:
    print(f"  NOT a 2-design (λ varies): {dict(pair_appearances)}")
    # Show the distribution
    for lam, count in sorted(Counter(pair_appearances.values()).items()):
        print(f"    λ = {lam}: {count} pairs")
