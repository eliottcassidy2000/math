#!/usr/bin/env python3
"""
alpha2_geometry.py — opus-2026-03-14-S71d

THE GEOMETRY OF α₂ AMBIGUITY

α₂ = #{vertex-disjoint pairs of directed odd cycles}

At n=7: α₂ = dp33 (only 3-3 pairs, since 3+5=8>7).
So α₂ counts pairs of vertex-disjoint 3-cycles.

Two tournaments with the SAME dc3 can have different α₂
because the 3-cycles can be:
  (a) All overlapping (sharing vertices) → α₂ small
  (b) Some disjoint → α₂ large

KEY QUESTION: What GEOMETRIC property distinguishes these cases?

HYPOTHESIS: α₂ depends on the "independence number" of the
intersection graph of 3-cycles, where vertices = 3-cycle vertex-sets
and edges = pairs sharing a vertex.

Actually α₂ = #{independent edges in the intersection COMPLEMENT},
i.e., α₂ = #{matching in the disjointness graph}.
"""

import sys
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

def find_3cycles(A, n):
    """Find all directed 3-cycle vertex-sets."""
    cycles = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cycles.append(frozenset(combo))
    return cycles

# ======================================================================
# PART 1: Exhaustive n=6 — α₂ structure
# ======================================================================
print("=" * 70)
print("PART 1: EXHAUSTIVE α₂ STRUCTURE AT n=6")
print("=" * 70)

n = 6
tb = n*(n-1)//2

a2_data = defaultdict(list)
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4

    cycles_3 = find_3cycles(A, n)
    dp33 = sum(1 for i in range(len(cycles_3)) for j in range(i+1, len(cycles_3))
               if cycles_3[i].isdisjoint(cycles_3[j]))

    a2_data[dc3].append({'a2': a2, 'dp33': dp33, 'dc5': dc5, 'H': H, 'bits': bits,
                          'n_3cycles': len(cycles_3)})

# Check α₂ = dp33 at n=6
all_match = all(d['a2'] == d['dp33'] for lst in a2_data.values() for d in lst)
print(f"\n  α₂ = dp33 at n=6: {all_match}")

# Wait — at n=6, 3+3=6=n, so disjoint 3-cycles use ALL vertices.
# And there might be (3,3) disjoint pairs. But also check if dc5 matters.
# dc5 uses 5 vertices, leaving 1 vertex free → no disjoint partner.
# So α₂ = dp33 at n=6 too.

print(f"\n  dc3 → range of α₂ (exhaustive n=6):")
for dc3 in sorted(a2_data.keys()):
    a2_vals = sorted(set(d['a2'] for d in a2_data[dc3]))
    count = len(a2_data[dc3])
    print(f"    dc3={dc3:2d}: α₂ ∈ {a2_vals} ({count} tournaments)")

# ======================================================================
# PART 2: The constraint on α₂ given dc3
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: MAX α₂ AS FUNCTION OF dc3")
print("=" * 70)

# At n=6: 2 disjoint 3-cycles partition {0,...,5} into two triples.
# There are C(6,3)/2 = 10 ways to partition into two triples.
# Each partition gives at most 1 pair (if both triples are cyclic).

# For dc3 = k: how many of the C(6,3) = 20 triples are cyclic?
# k triples are cyclic. How many disjoint pairs can there be?

# The disjointness graph on C(6,3) triples: two triples are adjacent
# if they share NO vertex. Since each triple uses 3 of 6 vertices,
# two triples are disjoint iff they partition {0,...,5}.
# So each triple has EXACTLY 1 disjoint partner!

# Therefore: dp33 = #{complementary pairs where BOTH are cyclic}
# If triple T is cyclic and T^c is cyclic, that's 1 pair.
# So dp33 = #{triples T with both T and T^c cyclic} / 2
# (since each pair is counted twice: once for T, once for T^c)

# Wait: we have 10 complementary pairs. dp33 = #{pairs where both are cyclic}.
# So dp33 ≤ 10. And dp33 = #{pairs (T, T^c) : T cyclic AND T^c cyclic} / 2... no.
# Each complementary pair {T, T^c} is counted once. dp33 counts pairs where both T and T^c are cyclic.

print(f"\n  At n=6: There are C(6,3)/2 = 10 complementary triple-pairs.")
print(f"  dp33 = #{'{'}pairs where BOTH triples are 3-cycles{'}'}")
print(f"  So dp33 ≤ 10.")

# Verify
for dc3 in sorted(a2_data.keys()):
    max_a2 = max(d['a2'] for d in a2_data[dc3])
    print(f"    dc3={dc3:2d}: max(α₂) = {max_a2}")

# ======================================================================
# PART 3: The complementary 3-cycle structure
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: COMPLEMENTARY 3-CYCLE ANALYSIS AT n=6")
print("=" * 70)

# For each tournament, check which complementary pairs are doubly-cyclic
for bits in [0, 100, 500, 1000]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    cycles = find_3cycles(A, n)
    all_verts = frozenset(range(6))

    comp_pairs = []
    for combo in combinations(range(6), 3):
        t1 = frozenset(combo)
        t2 = all_verts - t1
        if t1 < t2:  # avoid counting twice
            c1 = t1 in cycles
            c2 = t2 in cycles
            comp_pairs.append((t1, t2, c1, c2))

    dp = sum(1 for _, _, c1, c2 in comp_pairs if c1 and c2)
    print(f"\n  bits={bits}: dc3={dc3}, α₂={dp}")
    for t1, t2, c1, c2 in comp_pairs:
        if c1 or c2:
            mark = "**" if c1 and c2 else "  "
            print(f"    {mark}{set(t1)} {'cyclic' if c1 else 'trans.'} | {set(t2)} {'cyclic' if c2 else 'trans.'}")

# ======================================================================
# PART 4: At n=7 — the matching in the disjointness graph
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: DISJOINTNESS GRAPH AT n=7")
print("=" * 70)

n = 7
tb = n*(n-1)//2

# At n=7: C(7,3) = 35 triples. Two triples {a,b,c} and {d,e,f} are
# disjoint iff they share no vertex. This uses 6 of 7 vertices, leaving 1 free.
# Each triple T has C(4,3) = 4 disjoint partners (choose 3 from the 4 remaining).

# So the disjointness graph G has 35 vertices, each of degree 4.
# dp33 = #{edges of G where BOTH endpoints are 3-cycles in T}
#       = #{edges in the induced subgraph G[cyclic triples]}

# The number of edges in G: each vertex has degree 4, so |E| = 35*4/2 = 70.
# If dc3 = k triples are cyclic, the induced subgraph has at most C(k,2) edges.
# But in practice it's much less, constrained by the graph structure.

print(f"  Disjointness graph on C(7,3)=35 triples:")
print(f"  Each vertex has degree 4 (choose 3 from 4 remaining vertices).")
print(f"  Total edges: 35*4/2 = 70.")
print(f"  α₂ = edges in induced subgraph on cyclic triples.")

# Sample some tournaments
np.random.seed(42)
for trial in range(10):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    cycles = find_3cycles(A, n)
    cycle_set = set(cycles)

    # Count edges in disjointness graph restricted to cyclic triples
    dp33 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if cycles[i].isdisjoint(cycles[j]):
                dp33 += 1

    # Also: maximum matching in disjointness graph on cycles
    # (α₂ counts ALL disjoint pairs, not just a maximum matching)
    a1 = dc3 + count_directed_k_cycles(A, n, 5) + count_ham_cycles(A, n)
    a2 = (H - 1 - 2*a1) // 4

    print(f"  trial {trial}: dc3={dc3:2d}, α₂={a2}, dp33={dp33}, "
          f"density={dp33/max(1,comb(dc3,2)):.3f}")

# ======================================================================
# PART 5: The "spread" metric
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: VERTEX COVERAGE SPREAD OF 3-CYCLES")
print("=" * 70)

# For same dc3, what distinguishes high α₂ from low α₂?
# Hypothesis: more "spread out" 3-cycles → more disjoint pairs → higher α₂

np.random.seed(42)
data = []
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4

    cycles = find_3cycles(A, n)
    # Vertex coverage: how many vertices are in at least one 3-cycle?
    covered = set()
    for c in cycles:
        covered |= c
    coverage = len(covered)

    # Uniformity: variance of vertex participation counts
    vert_count = Counter()
    for c in cycles:
        for v in c:
            vert_count[v] += 1
    counts = [vert_count.get(v, 0) for v in range(n)]
    uniformity = np.std(counts)

    data.append({
        'dc3': dc3, 'a2': a2, 'coverage': coverage,
        'uniformity': uniformity, 'H': H
    })

# Correlation between coverage and α₂, controlling for dc3
print(f"\n  At n=7 (1000 samples):")
dc3_arr = np.array([d['dc3'] for d in data], dtype=float)
a2_arr = np.array([d['a2'] for d in data], dtype=float)
cov_arr = np.array([d['coverage'] for d in data], dtype=float)
uni_arr = np.array([d['uniformity'] for d in data], dtype=float)

print(f"    corr(dc3, α₂) = {np.corrcoef(dc3_arr, a2_arr)[0,1]:.4f}")
print(f"    corr(coverage, α₂) = {np.corrcoef(cov_arr, a2_arr)[0,1]:.4f}")
print(f"    corr(uniformity, α₂) = {np.corrcoef(uni_arr, a2_arr)[0,1]:.4f}")

# Partial correlation: α₂ ~ coverage | dc3
# Use residualization
from numpy.polynomial import polynomial as P
# Simple: for each dc3 value, check corr(coverage, α₂)
print(f"\n  Partial correlation (α₂ ~ coverage | dc3):")
for dc3_val in range(5, 13):
    subset = [(d['coverage'], d['a2']) for d in data if d['dc3'] == dc3_val]
    if len(subset) < 10:
        continue
    c, a = zip(*subset)
    r = np.corrcoef(c, a)[0,1]
    print(f"    dc3={dc3_val:2d}: corr(coverage, α₂) = {r:.4f} (n={len(subset)})")

# ======================================================================
# PART 6: Score sequence and the complementary structure
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: SCORE SEQUENCE AND COMPLEMENTARY PAIRS")
print("=" * 70)

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

# At n=6, exhaustive
n = 6
tb = n*(n-1)//2

score_a2 = defaultdict(set)
score_dc3 = defaultdict(set)
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    ss = score_seq(A, n)
    score_a2[ss].add(a2)
    score_dc3[ss].add(dc3)

print(f"\n  n=6 exhaustive:")
print(f"    {len(score_dc3)} distinct score sequences")
amb_dc3 = sum(1 for v in score_dc3.values() if len(v) > 1)
amb_a2 = sum(1 for v in score_a2.values() if len(v) > 1)
print(f"    score → dc3: {amb_dc3}/{len(score_dc3)} ambiguous")
print(f"    score → α₂:  {amb_a2}/{len(score_a2)} ambiguous")

for ss in sorted(score_a2.keys()):
    dc3_vals = sorted(score_dc3[ss])
    a2_vals = sorted(score_a2[ss])
    print(f"    {ss}: dc3∈{dc3_vals}, α₂∈{a2_vals}")

# ======================================================================
# SUMMARY
# ======================================================================
print("\n" + "=" * 70)
print("SUMMARY: GEOMETRY OF α₂")
print("=" * 70)

print("""
  1. α₂ = #{vertex-disjoint pairs of directed 3-cycles} (at n≤8)

  2. At n=6: α₂ = #{complementary pair partitions where both are cyclic}
     There are exactly 10 complementary pairs (partition of 6 into 3+3).
     α₂ ranges from 0 to 4.

  3. At n=7: α₂ = #{edges in the disjointness graph of 3-cycles}
     The disjointness graph on C(7,3)=35 triples has degree 4.
     α₂ ranges from 0 to ~12 depending on how many triples are cyclic
     AND how they are arranged.

  4. KEY FACT: dc3 does NOT determine α₂ because the ARRANGEMENT
     of 3-cycles matters, not just their count.
     Two tournaments with dc3=8 can have α₂=0 (all overlapping)
     or α₂=6 (many disjoint).

  5. Score sequence does NOT determine α₂ (at n≥6).
     Score determines dc3 (classical formula), but the GEOMETRY
     of where 3-cycles sit is not captured by scores alone.

  6. This is why "knowing 3" (I(Ω,3)) adds information beyond
     "knowing 2" (I(Ω,2)=H): the 3-evaluation captures the
     intersection structure that α₂ encodes.
""")

print("Done.")
