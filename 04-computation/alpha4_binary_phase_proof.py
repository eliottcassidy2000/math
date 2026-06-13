#!/usr/bin/env python3
"""
alpha4_binary_phase_proof.py — WHY does alpha_1=4 force alpha_2 ∈ {0,4}?

From kind-pasteur's Phase Transition Table:
  alpha_1=4 means exactly 4 directed odd cycles (as vertex sets, with multiplicity)
  alpha_2 counts pairwise-disjoint pairs among them
  Observed: alpha_2 ∈ {0,4} only — no values 1,2,3,5,6

Hypothesis: When a tournament has 4 odd-cycle vertex sets,
the intersection structure is either "all pairwise intersecting" (α₂=0)
or "C₄/K_{2,2} topology" (α₂=4 = two disjoint pairs each internally overlapping).

This script investigates the structural reason.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def get_directed_cycles(A, n):
    """Get all directed odd cycles as (vertex_set, count) pairs."""
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
    """Compute alpha_1 (= beta_1) and alpha_2 (= beta_2)."""
    vs_list = list(groups.items())
    n_vs = len(vs_list)

    # alpha_1 = sum of all d(S) = total directed cycles
    # Actually: alpha_1 = number of cycle vertex sets (with multiplicity)
    # Wait - need to clarify. In the beta formulation:
    # beta_1 = sum over vertex sets S of d(S), where d(S) = #directed cycles on S
    # beta_2 = sum over disjoint pairs (S,S') of d(S)*d(S')

    alpha1 = sum(d for _, d in vs_list)

    alpha2 = 0
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if not (vs_list[i][0] & vs_list[j][0]):  # disjoint
                alpha2 += vs_list[i][1] * vs_list[j][1]

    return alpha1, alpha2

print("=" * 70)
print("ALPHA_1=4 BINARY PHASE THEOREM INVESTIGATION")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: At n=6, analyze all tournaments with alpha_1=4
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: alpha_1=4 structure at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

alpha4_cases = []  # (bits, groups_dict, alpha2)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    if a1 == 4:
        alpha4_cases.append((bits, dict(groups), a2))

print(f"  Tournaments with alpha_1=4: {len(alpha4_cases)}")
a2_dist = Counter(a2 for _, _, a2 in alpha4_cases)
print(f"  alpha_2 distribution: {dict(sorted(a2_dist.items()))}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: Analyze the cycle vertex set intersection patterns
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: Intersection patterns ---")

# Classify by intersection graph type
patterns = defaultdict(list)

for bits, groups, a2 in alpha4_cases:
    vs_list = list(groups.keys())
    ds = [groups[vs] for vs in vs_list]

    # How many vertex sets, and what are their sizes?
    sizes = sorted([len(vs) for vs in vs_list])

    # Build intersection graph
    n_vs = len(vs_list)
    int_edges = []
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if vs_list[i] & vs_list[j]:
                int_edges.append((i,j))

    # Key: (sizes, multiplicities, #intersecting_pairs, alpha_2)
    mults = sorted(ds)
    key = (tuple(sizes), tuple(mults), len(int_edges), a2)
    patterns[key].append(bits)

print(f"  Distinct structural patterns: {len(patterns)}")
for key, cases in sorted(patterns.items()):
    sizes, mults, int_pairs, a2 = key
    n_vs = len(sizes)
    total_pairs = n_vs * (n_vs - 1) // 2
    disj_pairs = total_pairs - int_pairs
    print(f"    sizes={sizes}, d={mults}, "
          f"int_pairs={int_pairs}/{total_pairs}, "
          f"alpha_2={a2}, count={len(cases)}")

    # Show first example in detail
    if len(cases) > 0:
        bits = cases[0]
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
        groups = get_directed_cycles(A, n)
        print(f"      Example (bits={bits}):")
        for vs, d in sorted(groups.items(), key=lambda x: (len(x[0]), x[0])):
            print(f"        {sorted(vs)} d={d}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: Why can't we have exactly 1, 2, or 3 disjoint pairs?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: Why alpha_2 ∈ {0,4} only? ---")

# Key question: if we have 4 cycle vertex sets with alpha_1=4,
# that means each vertex set has d=1 (one directed cycle each),
# OR some have d>1 and fewer vertex sets exist.

# Case A: 4 vertex sets, each with d=1
# Then alpha_2 = #disjoint pairs
# With 4 sets, C(4,2)=6 possible pairs
# The intersection graph is a graph on 4 vertices
# alpha_2 = 6 - #edges in intersection graph

# For alpha_2=0: all 6 pairs intersect → K₄ intersection graph (complete)
# For alpha_2=4: 2 pairs intersect, 4 disjoint → intersection graph has 2 edges
# That's a matching (2 disjoint edges) → K_{2,2} bipartite → "two disjoint pairs internally overlapping"

# For alpha_2=1: 5 pairs intersect → K₄ minus 1 edge → almost complete
# For alpha_2=2: 4 pairs intersect → C₄ or K₃∪K₁ etc
# For alpha_2=3: 3 pairs intersect → matching or triangle

# Case B: fewer vertex sets with higher d
# 2 vertex sets with d=2 → alpha_1 = 2+2 = 4
# 1 vertex set with d=2, 2 with d=1 → alpha_1 = 2+1+1 = 4
# 1 vertex set with d=3, 1 with d=1 → alpha_1 = 3+1 = 4
# 1 vertex set with d=4 → alpha_1 = 4

# Let's check which cases actually occur
print("  Case analysis of alpha_1=4 decomposition:")
decomp_types = Counter()
for bits, groups, a2 in alpha4_cases:
    ds = tuple(sorted(groups.values(), reverse=True))
    decomp_types[(ds, a2)] += 1

for (ds, a2), count in sorted(decomp_types.items()):
    print(f"    d-vector={ds}, alpha_2={a2}: {count} tournaments")

# ═══════════════════════════════════════════════════════════════════
# Part 4: Focus on 4 vertex sets with d=1 each (the interesting case)
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: 4 vertex sets with d=(1,1,1,1) ---")

case_4111 = [(bits, groups, a2) for bits, groups, a2 in alpha4_cases
             if len(groups) == 4 and all(d == 1 for d in groups.values())]

print(f"  Count: {len(case_4111)}")
a2_dist_4111 = Counter(a2 for _, _, a2 in case_4111)
print(f"  alpha_2 distribution: {dict(sorted(a2_dist_4111.items()))}")

# For each, what is the intersection graph?
print("\n  Intersection graph types:")
int_graph_types = Counter()
for bits, groups, a2 in case_4111:
    vs_list = list(groups.keys())
    # Intersection graph: which pairs of vertex sets share a vertex
    int_adj = set()
    for i in range(4):
        for j in range(i+1, 4):
            if vs_list[i] & vs_list[j]:
                int_adj.add((i,j))
    int_graph_types[(frozenset(int_adj), a2)] += 1

for (edges_set, a2), count in sorted(int_graph_types.items(), key=lambda x: x[1], reverse=True):
    n_edges = len(edges_set)
    print(f"    {n_edges} intersecting pairs (alpha_2={6-n_edges}={a2}): {count} cases")

# ═══════════════════════════════════════════════════════════════════
# Part 5: The SPLICING argument for why 1,2,3,5 disjoint pairs are impossible
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: Why certain intersection graphs are impossible ---")

# If exactly 4 vertex sets exist with d=1 each, and they're all 3-cycles (most common),
# let's see what constraint the tournament structure imposes.

# Key insight from Splicing Lemma: two odd cycles sharing exactly 1 vertex
# produce a 3rd cycle. If we start with 4 3-cycles and they share vertices,
# the splicing creates additional cycles, potentially violating alpha_1=4.

# Let's trace this for a concrete example
print("  Checking: do intersecting 3-cycles always splice to create more?")

for bits, groups, a2 in case_4111[:5]:  # first 5 examples
    vs_list = list(groups.keys())
    sizes = [len(vs) for vs in vs_list]

    # Check pairwise intersections
    for i in range(4):
        for j in range(i+1, 4):
            overlap = vs_list[i] & vs_list[j]
            if overlap and len(overlap) == 1:
                # Splicing should create a 3rd cycle!
                union_size = len(vs_list[i] | vs_list[j])
                print(f"    bits={bits}: {sorted(vs_list[i])} ∩ {sorted(vs_list[j])} = {sorted(overlap)}"
                      f" → splice creates {union_size}-cycle on {sorted(vs_list[i] | vs_list[j])}")

# ═══════════════════════════════════════════════════════════════════
# Part 6: Check at n=7 (sample)
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 6: alpha_1=4 at n=7 (sampling) ---")

import random
random.seed(42)

n = 7
edges7 = [(i,j) for i in range(n) for j in range(i+1,n)]
ne7 = len(edges7)

a4_count = 0
a2_vals_7 = Counter()
decomp_7 = Counter()

SAMPLES = 200000
for _ in range(SAMPLES):
    bits = random.randint(0, 2**ne7 - 1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges7):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    if a1 == 4:
        a4_count += 1
        a2_vals_7[a2] += 1
        ds = tuple(sorted(groups.values(), reverse=True))
        decomp_7[(ds, a2)] += 1

print(f"  Sampled {SAMPLES} tournaments, {a4_count} with alpha_1=4")
print(f"  alpha_2 distribution: {dict(sorted(a2_vals_7.items()))}")
print(f"  Decomposition types:")
for (ds, a2), count in sorted(decomp_7.items()):
    print(f"    d={ds}, alpha_2={a2}: {count}")

# ═══════════════════════════════════════════════════════════════════
# Part 7: The KEY structural theorem
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 7: Structural theorem ---")
print("""
CLAIM: For alpha_1=4 with 4 vertex sets each having d=1,
the intersection graph is either K_4 (all intersect) or
the complement of K_4 (all disjoint) or has exactly 2 edges
forming a matching.

REASON: If two 3-cycles share exactly 1 vertex, Splicing creates
a 5th cycle → alpha_1 ≥ 5. So any intersection must involve
sharing ≥2 vertices. But two 3-cycles sharing 2 vertices must be
on vertex sets like {a,b,c} and {a,b,d}, sharing {a,b}.
These DON'T splice (they share 2 vertices, not 1).

The constraint is: no two of the 4 cycles can share EXACTLY 1 vertex
(otherwise splicing creates a 5th cycle, violating alpha_1=4).

With this constraint, each pair either:
- Shares 0 vertices (disjoint)
- Shares 2 vertices (heavily overlapping)

4 triangles on 6 vertices with this constraint:
If all 4 pairwise share ≥2 vertices → all on ≤4 vertices → K_4 type
If some pairs disjoint → disjoint pairs = alpha_2
""")

# Verify: for all alpha_1=4 cases at n=6, no two cycles share exactly 1 vertex
print("  VERIFICATION: Do any intersecting cycle pairs share exactly 1 vertex?")
single_overlap_count = 0
for bits, groups, a2 in alpha4_cases:
    vs_list = list(groups.keys())
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            overlap = vs_list[i] & vs_list[j]
            if len(overlap) == 1:
                single_overlap_count += 1

print(f"  Single-vertex overlaps found: {single_overlap_count}")
if single_overlap_count == 0:
    print("  CONFIRMED: No alpha_1=4 tournament has two cycles sharing exactly 1 vertex!")
    print("  This is BECAUSE sharing 1 vertex → Splicing → 5th cycle → alpha_1≥5")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
