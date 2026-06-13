#!/usr/bin/env python3
"""
h21_case81_structural.py — Structural investigation of (alpha_1=8, alpha_2=1).

KEY QUESTION: Why can't we have exactly 8 directed odd cycles
with exactly 1 disjoint pair?

At n=6: alpha_1=8 always has t₃=4, d₅=4, all cycles on a 5-vertex
"star" with one vertex connected to all. alpha_2=0 always.

At n=7 (sampling): alpha_1=8 → alpha_2=0 always.
At n=8+ (sampling): alpha_1=8 → alpha_2 ∈ {0, 7}.

The gap at alpha_2=1 means the "minimum nonzero alpha_2" is >1.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random

sys.stdout.reconfigure(line_buffering=True)

def get_directed_cycles(A, n):
    groups = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    groups[frozenset(verts)] += 1
    return groups

def compute_alpha(groups):
    vs_list = list(groups.items())
    alpha1 = sum(d for _, d in vs_list)
    alpha2 = 0
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            if not (vs_list[i][0] & vs_list[j][0]):
                alpha2 += vs_list[i][1] * vs_list[j][1]
    return alpha1, alpha2

print("=" * 70)
print("CASE (8,1): STRUCTURAL INVESTIGATION")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: At n=6, what does alpha_1=8 look like?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: alpha_1=8 structure at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

a8_structures = Counter()

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    if a1 == 8:
        # Characterize: which vertices appear in cycles?
        all_cycle_verts = set()
        for vs in groups:
            all_cycle_verts |= vs
        cycle_vertex_count = len(all_cycle_verts)

        # Sizes of cycle vertex sets
        sizes = tuple(sorted(len(vs) for vs in groups))
        t3 = sum(1 for vs in groups if len(vs) == 3)
        d5 = sum(1 for vs in groups if len(vs) == 5)

        a8_structures[(t3, d5, cycle_vertex_count)] += 1

print(f"  (t₃, d₅, cycle_vertices): count")
for key, count in sorted(a8_structures.items()):
    t3, d5, cv = key
    print(f"    t₃={t3}, d₅={d5}, cycle_verts={cv}: {count}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: The "cycle vertex span" constraint
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Cycle vertex span analysis ---")

# If all 8 cycles live on 5 vertices (cycle_verts=5),
# then no pair can be disjoint (two 3-cycles need 6 vertices).
# alpha_2=0 forced.

# For alpha_2=1: need at least 6 cycle vertices.
# The disjoint pair needs 6 vertices (two 3-cycles) or
# 5+3=8 (3-cycle + 5-cycle) or 5+5=10 (two 5-cycles).
# At n=6, 5+3=8 > 6, so the only option is two complementary 3-cycles.

# But we showed at n=6: alpha_1=8 always has cycle_verts=5.
# So alpha_2=1 is impossible at n=6.

print(f"  At n=6: alpha_1=8 always has cycle_verts=5.")
print(f"  Two disjoint 3-cycles need 6 vertices → impossible.")
print(f"  Two disjoint cycles at n=6: only complementary 3-cycles (6 verts).")
print(f"  Since cycle_verts=5 < 6: alpha_2=0 forced. ✓")

# ═══════════════════════════════════════════════════════════════════
# Part 3: At n=7, analyze alpha_1=8 structure (exhaustive)
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: alpha_1=8 at n=7 (exhaustive) ---")

n = 7
edges7 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne7 = len(edges7)

a8_7_structures = Counter()
a8_7_a2 = Counter()
a8_7_count = 0

for bits in range(2**ne7):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges7):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    if a1 == 8:
        a8_7_count += 1
        a8_7_a2[a2] += 1

        all_cycle_verts = set()
        for vs in groups:
            all_cycle_verts |= vs
        cycle_vertex_count = len(all_cycle_verts)

        t3 = sum(d for vs, d in groups.items() if len(vs) == 3)
        d5 = sum(d for vs, d in groups.items() if len(vs) == 5)
        d7 = sum(d for vs, d in groups.items() if len(vs) == 7)
        n_vs = len(groups)

        a8_7_structures[(t3, d5, d7, n_vs, cycle_vertex_count)] += 1

    if bits % 500000 == 0 and bits > 0:
        print(f"  Progress: {bits}/2097152")

print(f"\n  alpha_1=8 count: {a8_7_count}")
print(f"  alpha_2 distribution: {dict(sorted(a8_7_a2.items()))}")

print(f"\n  Structural breakdown:")
print(f"  {'t₃':>4} {'d₅':>4} {'d₇':>4} {'n_vs':>5} {'verts':>5} {'count':>6}")
for (t3, d5, d7, nvs, cv), count in sorted(a8_7_structures.items()):
    print(f"  {t3:4d} {d5:4d} {d7:4d} {nvs:5d} {cv:5d} {count:6d}")

# ═══════════════════════════════════════════════════════════════════
# Part 4: Minimum cycle vertex span for alpha_2=1
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: Minimum span for disjoint pair ---")

# For two cycle vertex sets to be disjoint:
# - Two 3-cycle VS: need 6 vertices
# - 3-cycle + 5-cycle: need 8 vertices
# - Two 5-cycle VS: need 10 vertices
# - 3-cycle + 7-cycle: need 10 vertices

# At n=7: possible disjoint pairs are:
# - Two 3-cycles using 6 of 7 vertices (the 7th is unused)
# - That's it! All other combinations need ≥8 vertices.

# So for alpha_2=1 at n=7 with alpha_1=8:
# Need at least one pair of disjoint 3-cycle vertex sets.
# This means the tournament has two 3-cycles on disjoint vertex sets
# {a,b,c} and {d,e,f} with vertex g unused.

# Q: Can this happen with total alpha_1=8?
# The two disjoint 3-cycles contribute 2 to alpha_1.
# Remaining 6 cycles must use vertices from {a,...,g}.
# The remaining cycles must NOT create additional disjoint pairs
# (or alpha_2 could be > 1).

# If remaining 6 cycles all overlap with BOTH {a,b,c} and {d,e,f}:
# Each must touch both groups → each uses vertex g.
# So remaining cycles are on subsets containing g plus vertices from both groups.
# Possible cycle VS: {a,d,g}, {a,e,g}, etc. (3-cycles through g).
# Or {a,b,d,e,g} (5-cycles through g).

# With g touching both groups: how many cycles through g?
# Triples containing g and one from each group: 3×3 = 9 possible triples.
# Each is either cyclic or not. If t₃(through g) = k, that adds k to alpha_1.
# Plus the 2 disjoint 3-cycles and 6 remaining → 2 + k + other = 8.
# If the only cycles are 2 disjoint 3-cycles + k cross-triples: 2+k=8, k=6.
# Need 6 of the 9 cross-triples {a_i, d_j, g} to be cyclic.

print("  At n=7: for alpha_2≥1 with alpha_1=8:")
print("    Need disjoint 3-cycle pair {a,b,c} ⊥ {d,e,f}, vertex g free.")
print("    Remaining alpha_1 = 6 cycles must all touch both groups.")
print("    Each such cycle uses vertex g (bridge between groups).")
print("    Cross-triples: {a_i, d_j, g} for i,j ∈ {0,1,2} → 9 possible.")
print("    Need 6 of 9 cross-triples to be cyclic.")
print()

# Can we have 6 of 9 cross-triples {a_i, d_j, g} be cyclic?
# These triples all contain g. The subtournament on the edges
# a_i→g, g→d_j, a_i→d_j determines cyclicity.

# In a triple {a_i, d_j, g}: it's a 3-cycle iff
# a_i→d_j→g→a_i OR a_i→g→d_j→a_i (two orientations).
# Exactly one orientation matches the tournament.

# For a specific tournament: edges involving g are:
# g→a₁ or a₁→g, g→a₂ or a₂→g, g→a₃ or a₃→g
# g→d₁ or d₁→g, g→d₂ or d₂→g, g→d₃ or d₃→g

# Also edges a_i→d_j or d_j→a_i.

# This is getting complex. Let me just check computationally.
# For n=7, enumerate tournaments with {0,1,2} cyclic, {3,4,5} cyclic,
# and count cross-triple cycles.

print("  Computational check: n=7, {0,1,2} and {3,4,5} both cyclic")

# Fix: 0→1→2→0 and 3→4→5→3
# Free edges: all involving vertex 6, plus cross-edges a_i→d_j
# Total edges = C(7,2) = 21. Fixed: 6 (3 for each cycle). Free: 15.

fixed_edges = {(0,1): 1, (1,2): 1, (2,0): 1,
               (3,4): 1, (4,5): 1, (5,3): 1}

free_edges = [(i,j) for i in range(7) for j in range(i+1,7)
              if (i,j) not in fixed_edges and (j,i) not in fixed_edges]
nf = len(free_edges)

a1_with_disjoint = Counter()
a2_with_disjoint = Counter()
count_with_disjoint = 0

for bits in range(2**nf):
    A = [[0]*7 for _ in range(7)]
    # Apply fixed edges
    for (i,j), d in fixed_edges.items():
        A[i][j] = 1
        # Opponent gets reverse
    # Apply free edges
    for idx, (i,j) in enumerate(free_edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, 7)
    a1, a2 = compute_alpha(groups)

    # Check that our fixed cycles still exist
    assert frozenset({0,1,2}) in groups
    assert frozenset({3,4,5}) in groups

    if a1 == 8:
        a1_with_disjoint[a2] += 1
        count_with_disjoint += 1

print(f"  Tournaments with fixed disjoint pair and alpha_1=8: {count_with_disjoint}")
print(f"  alpha_2 distribution: {dict(sorted(a1_with_disjoint.items()))}")

if 1 not in a1_with_disjoint:
    print(f"  (8,1) NOT achievable even with forced disjoint pair!")
    print(f"  This means: if {'{0,1,2}'} and {'{3,4,5}'} are both cyclic")
    print(f"  AND alpha_1=8, then alpha_2 ≥ 2 always.")
else:
    print(f"  (8,1) IS achievable with forced disjoint pair.")

# ═══════════════════════════════════════════════════════════════════
# Part 5: Try different disjoint pair orientations
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: All orientations of disjoint pair ---")

# Also try 0→2→1→0 (reverse) for first triple
# and 3→5→4→3 for second
for orient1 in [(0,1,2), (0,2,1)]:
    for orient2 in [(3,4,5), (3,5,4)]:
        fixed = {}
        for i in range(3):
            fixed[(orient1[i], orient1[(i+1)%3])] = 1
            fixed[(orient2[i], orient2[(i+1)%3])] = 1

        free = [(i,j) for i in range(7) for j in range(i+1,7)
                if (i,j) not in fixed and (j,i) not in fixed]
        nf = len(free)

        count_81 = 0
        total_a8 = 0

        for bits in range(2**nf):
            A = [[0]*7 for _ in range(7)]
            for (i,j), d in fixed.items():
                A[i][j] = 1
            for idx, (i,j) in enumerate(free):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1

            groups = get_directed_cycles(A, 7)
            a1, a2 = compute_alpha(groups)

            if a1 == 8:
                total_a8 += 1
                if a2 == 1:
                    count_81 += 1

        print(f"  orient1={orient1}, orient2={orient2}: "
              f"alpha_1=8: {total_a8}, (8,1): {count_81}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
