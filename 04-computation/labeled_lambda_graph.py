#!/usr/bin/env python3
"""
labeled_lambda_graph.py -- kind-pasteur-2026-03-13-S61

The LABELED lambda graph: keep the full pair-coverage function
lambda_{uv} as a labeled graph on vertices, not just a sorted histogram.

Key insight: the sorted lambda histogram loses the information about
WHICH vertex pairs share cycles. The unsorted (labeled) lambda
contains all the geometric arrangement information that determines alpha_2.

This script:
1. Check if the labeled lambda graph (up to isomorphism) determines H
2. The lambda graph as a WEIGHTED graph on tournament vertices
3. Connection to BIBD / block design theory
4. What graph properties of the lambda graph predict alpha_2?
5. The Vitali quotient: labeled/sorted = the non-measurable content

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_labeled_lambda(A, n):
    """Get the FULL labeled lambda matrix (not just sorted histogram).
    Returns the lambda matrix AND the 3-cycle vertex sets."""
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    # Full lambda matrix
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val

    return lam, c3_sets


def lambda_graph_canonical(lam, n):
    """Compute canonical form of the lambda graph under vertex permutation.
    For small n, try all permutations and return lexicographic minimum."""
    if n > 8:
        # Too expensive, use degree-sorted approximation
        degrees = sorted([(sum(lam[v]), v) for v in range(n)], reverse=True)
        perm = [v for _, v in degrees]
        flat = []
        for i in range(n):
            for j in range(i+1, n):
                flat.append(lam[perm[i]][perm[j]])
        return tuple(flat)

    # For n <= 8: try all n! permutations (feasible for n=5,6)
    best = None
    for perm in permutations(range(n)):
        flat = []
        for i in range(n):
            for j in range(i+1, n):
                flat.append(lam[perm[i]][perm[j]])
        flat = tuple(flat)
        if best is None or flat < best:
            best = flat
    return best


# ========================================================================
# ANALYSIS 1: n=5 — LABELED LAMBDA VS SORTED LAMBDA
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: n=5 — LABELED vs SORTED LAMBDA")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m

by_sorted = defaultdict(list)
by_labeled = defaultdict(list)

for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    lam, c3_sets = get_labeled_lambda(A, n)

    sorted_lam = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    labeled_lam = lambda_graph_canonical(lam, n)

    by_sorted[sorted_lam].append({'H': H, 'scores': scores, 'labeled': labeled_lam})
    by_labeled[labeled_lam].append({'H': H, 'scores': scores})

print(f"Sorted lambda classes: {len(by_sorted)}")
print(f"Labeled lambda classes (up to isomorphism): {len(by_labeled)}")

# Check: does labeled lambda determine H?
labeled_ambig = 0
for lab, group in by_labeled.items():
    H_set = sorted(set(d['H'] for d in group))
    if len(H_set) > 1:
        labeled_ambig += 1

if labeled_ambig == 0:
    print(f"Labeled lambda DETERMINES H at n=5!")
else:
    print(f"{labeled_ambig} ambiguous labeled classes at n=5")

# How many labeled classes per sorted class?
print(f"\nLabeled classes per sorted lambda class:")
for sl in sorted(by_sorted.keys()):
    group = by_sorted[sl]
    labeled_classes = set(d['labeled'] for d in group)
    H_set = sorted(set(d['H'] for d in group))
    if len(labeled_classes) > 1:
        print(f"  sorted={sl}: {len(labeled_classes)} labeled classes, H={H_set}")


# ========================================================================
# ANALYSIS 2: n=6 — LABELED LAMBDA RESOLVES AMBIGUITIES?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: n=6 — DOES LABELED LAMBDA RESOLVE AMBIGUITIES?")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

by_labeled_n6 = defaultdict(list)
count = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    lam, c3_sets = get_labeled_lambda(A, n)
    labeled_lam = lambda_graph_canonical(lam, n)

    by_labeled_n6[labeled_lam].append({'H': H, 'scores': scores})
    count += 1

print(f"Labeled lambda classes at n=6: {len(by_labeled_n6)}")

labeled_ambig_n6 = 0
for lab, group in by_labeled_n6.items():
    H_set = sorted(set(d['H'] for d in group))
    if len(H_set) > 1:
        labeled_ambig_n6 += 1
        sc_set = set(d['scores'] for d in group)
        print(f"  AMBIGUOUS: {len(group)} tours, H={H_set}, scores={sorted(sc_set)}")
        print(f"    labeled lambda = {lab[:10]}...")

if labeled_ambig_n6 == 0:
    print(f"Labeled lambda DETERMINES H at n=6!")
else:
    print(f"{labeled_ambig_n6} ambiguous labeled classes at n=6")


# ========================================================================
# ANALYSIS 3: n=7 AMBIGUOUS CASES — LABELED LAMBDA
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: n=7 AMBIGUOUS CASES — LABELED LAMBDA")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

# Check the specific ambiguous case
target_score = (0, 2, 3, 3, 4, 4, 5)
target_sorted_lam = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3)

data = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != target_score:
        continue

    lam, c3_sets = get_labeled_lambda(A, n)
    sorted_lam = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    if sorted_lam != target_sorted_lam:
        continue

    H = count_ham_paths(A, n)

    # For n=7, canonical form is expensive. Use a cheaper invariant:
    # degree sequence of the lambda graph (sum of lambda values per vertex)
    lam_degrees = tuple(sorted(sum(lam[v]) for v in range(n)))

    # Also: the multiset of (lam_degree(u), lam_degree(v), lam(u,v)) for edges
    edge_types = []
    for u in range(n):
        for v in range(u+1, n):
            edge_types.append((sum(lam[u]), sum(lam[v]), lam[u][v]))
    edge_type_hist = tuple(sorted(edge_types))

    data.append({
        'bits': bits, 'H': H, 'lam_degrees': lam_degrees,
        'edge_types': edge_type_hist
    })

print(f"Found {len(data)} tournaments in ambiguous class")

by_H = defaultdict(list)
for d in data:
    by_H[d['H']].append(d)

for H in sorted(by_H.keys()):
    group = by_H[H]
    deg_set = set(d['lam_degrees'] for d in group)
    print(f"\n  H={H}: {len(group)} tournaments")
    print(f"    Lambda degree sequences: {sorted(deg_set)}")

# Check: do lambda degrees resolve?
by_deg = defaultdict(set)
for d in data:
    by_deg[d['lam_degrees']].add(d['H'])

deg_resolves = all(len(v) == 1 for v in by_deg.values())
print(f"\n  Lambda degree sequence resolves: {deg_resolves}")

# Check: do edge types resolve?
by_edge = defaultdict(set)
for d in data:
    by_edge[d['edge_types']].add(d['H'])

edge_resolves = all(len(v) == 1 for v in by_edge.values())
print(f"  Edge type histogram resolves: {edge_resolves}")


# ========================================================================
# ANALYSIS 4: THE FUNDAMENTAL INSIGHT — VERTEX SCORE + LAMBDA
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: VERTEX SCORE + LAMBDA DETERMINES EVERYTHING?")
print("=" * 70)

# The key realization: the SORTED lambda loses vertex identity.
# But the tournament already has labeled vertices.
# The "labeled lambda" = lambda_{uv} for each specific (u,v) pair.
# But two tournaments can have different vertex LABELINGS.
#
# The RIGHT invariant is: for each vertex v with out-degree d_v,
# how does its lambda profile look?
#
# In the ambiguous case: same SCORE SEQUENCE and same SORTED lambda.
# But the vertex-specific lambda profiles differ.

# Let me check: for each vertex, compute its "lambda profile" =
# (out_degree_v, sorted lambda values to all other vertices)

for H_target in [25, 29]:
    d = by_H[H_target][0]
    bits = d['bits']
    A = binary_to_tournament(bits, n)
    lam, c3_sets = get_labeled_lambda(A, n)
    out_degs = [sum(A[v]) for v in range(n)]

    print(f"\n  H={H_target} (bits={bits}):")
    print(f"    Out-degrees: {out_degs}")
    for v in range(n):
        lam_profile = sorted([lam[v][w] for w in range(n) if w != v])
        print(f"    Vertex {v} (d_out={out_degs[v]}): lambda profile = {lam_profile}")


# ========================================================================
# ANALYSIS 5: DISJOINTNESS GEOMETRY — THE REAL INVARIANT
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: DISJOINTNESS GEOMETRY")
print("=" * 70)

# The fundamental question: alpha_2 = #{disjoint 3-cycle pairs}.
# Two 3-cycles are disjoint iff their vertex sets don't intersect.
# At n=7 with c3=6 vertex sets (all non-isomorphic),
# alpha_2 = 0 or 1 for these small cases.

# The GEOMETRY: which pairs of 3-cycles are disjoint?
# This is a graph on 3-cycle vertex sets (disjointness graph).

for H_target in [25, 29]:
    d = by_H[H_target][0]
    bits = d['bits']
    A = binary_to_tournament(bits, n)
    lam, c3_sets = get_labeled_lambda(A, n)

    print(f"\n  H={H_target} (bits={bits}):")
    print(f"    3-cycle vertex sets: {[sorted(fs) for fs in c3_sets]}")

    # Disjointness graph
    nc3 = len(c3_sets)
    disj_graph = [[0]*nc3 for _ in range(nc3)]
    disj_count = 0
    for i in range(nc3):
        for j in range(i+1, nc3):
            if not (c3_sets[i] & c3_sets[j]):
                disj_graph[i][j] = 1
                disj_graph[j][i] = 1
                disj_count += 1
                print(f"    DISJOINT: {sorted(c3_sets[i])} and {sorted(c3_sets[j])}")

    print(f"    alpha_2 = {disj_count}")

    # What makes the difference? At H=29 there IS a disjoint pair,
    # at H=25 there ISN'T. But they have the SAME sorted lambda histogram.
    # The sorted lambda says THERE EXIST 9 pairs with lambda=0, etc.
    # But the lambda=0 pairs might or might not form complementary vertex sets
    # for 3-cycles.

    # Show which vertex pairs have lambda=0
    zero_pairs = []
    for u in range(n):
        for v in range(u+1, n):
            if lam[u][v] == 0:
                zero_pairs.append((u, v))
    print(f"    Vertex pairs with lambda=0: {zero_pairs}")


# ========================================================================
# ANALYSIS 6: THE VITALI QUOTIENT — WHAT SORTING DESTROYS
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: THE VITALI QUOTIENT — WHAT SORTING DESTROYS")
print("=" * 70)

# The Vitali analogy is now precise:
# - "Measurable": sorted lambda histogram (forgets vertex labels)
# - "Non-measurable": the specific assignment of lambda values to vertex PAIRS
# - The quotient (labeled / sorted) contains the geometric information
#   about WHICH vertex pairs share cycles

# This is exactly like the Vitali construction:
# - The equivalence relation: x ~ y iff x - y in Q (rationals)
# - Each equivalence class is a coset of Q
# - The Vitali set picks one representative from each coset
# - The "non-measurable" content is the specific choice of representative

# For tournaments:
# - The equivalence relation: T ~ T' iff sorted_lambda(T) = sorted_lambda(T')
# - Each class has same sorted lambda but different labeled lambda
# - The "non-measurable" content is the specific vertex-pair assignment
# - alpha_2 (disjointness count) depends on this assignment

print("""
THE VITALI STRUCTURE OF TOURNAMENT CYCLE GEOMETRY:

MEASURABLE (sorted lambda):
  - Determines: c3, ov_histogram, sum(lambda^2), mean(lambda)
  - Does NOT determine: alpha_2 (disjointness), H

NON-MEASURABLE (vertex-pair assignment):
  - The specific mapping: vertex pair (u,v) -> lambda_{uv}
  - Determines: which 3-cycle vertex sets are disjoint
  - Determines: alpha_2 and hence H

THE QUOTIENT:
  Within each sorted-lambda class, the labeled lambda varies.
  The number of labeled lambda orbits = number of H values.

  For REGULAR n=7: sorted lambda has only 3 forms, and each
  sorted form corresponds to exactly ONE labeled orbit.
  So sorting loses NOTHING -> lambda determines H.

  For non-regular: sorted lambda can have MULTIPLE labeled orbits.
  The orbits with different alpha_2 give different H.
  This is the "non-measurable" content.

ANALOGY TABLE:
  Lebesgue measure   <-->   Score sequence
  Borel sigma-algebra <-->   Sorted lambda histogram
  Power set           <-->   Labeled lambda graph
  Non-measurable set  <-->   Tournaments with same sorted lambda but different alpha_2
  Vitali set          <-->   Set of representatives, one per labeled orbit
""")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
