#!/usr/bin/env python3
"""
complete_invariant_search.py -- kind-pasteur-2026-03-13-S61

Does the labeled lambda graph (up to isomorphism) determine H at n=7?
If yes, this is a MAJOR theorem: H is determined by the pair-coverage
function lambda_{uv} (the labeled lambda graph).

Also: systematize the {2,1,0} overlap weight connection.

Strategy:
- At n=7, computing full graph isomorphism canonical form for all 2M
  tournaments is expensive. Instead, use INVARIANTS of the labeled
  lambda graph to check if they resolve all ambiguities.

Invariants of the labeled lambda graph (cheapest to most expensive):
  Level 1: Sorted histogram (INSUFFICIENT — proved)
  Level 2: Lambda degree sequence (row sums, sorted)
  Level 3: Lambda 2-step walk counts (paths of length 2)
  Level 4: Lambda spectrum (eigenvalues) — needs float approx
  Level 5: Full canonical form (exact but expensive)

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
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


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val

    return lam, c3_sets


def lambda_degree_sequence(lam, n):
    """Sorted row sums of the lambda matrix."""
    return tuple(sorted(sum(lam[v]) for v in range(n)))


def lambda_walk_invariant(lam, n):
    """Compute invariants from 2-step walks in the lambda graph.
    For each vertex v, compute the sum of lam[v][w]*lam[w][u] over all w,u.
    This gives a "2-neighborhood" invariant."""
    # L^2[v][u] = sum_w lam[v][w] * lam[w][u]
    L2 = [[0]*n for _ in range(n)]
    for v in range(n):
        for u in range(n):
            s = 0
            for w in range(n):
                s += lam[v][w] * lam[w][u]
            L2[v][u] = s

    # Invariant: sorted diagonal of L^2 (= degree-weighted degree)
    diag_L2 = tuple(sorted(L2[v][v] for v in range(n)))

    # Trace of L^2
    trace_L2 = sum(L2[v][v] for v in range(n))

    # Sorted row sums of L^2
    row_sums_L2 = tuple(sorted(sum(L2[v]) for v in range(n)))

    return diag_L2, trace_L2, row_sums_L2


def lambda_combined_invariant(lam, n, scores):
    """Combined invariant: lambda degree sequence + walk invariants + score correlation."""
    deg_seq = lambda_degree_sequence(lam, n)
    diag_L2, trace_L2, row_sums_L2 = lambda_walk_invariant(lam, n)

    # Also: correlation between out-degree and lambda-degree
    # lam_degs are UNsorted here, paired with out-degrees
    lam_degs = [sum(lam[v]) for v in range(n)]

    # Joint distribution: (out_degree, lam_degree) pairs, sorted
    joint = tuple(sorted(zip(scores, lam_degs)))

    return (deg_seq, diag_L2, trace_L2, row_sums_L2, joint)


# ========================================================================
# ANALYSIS 1: n=7 — DOES LAMBDA DEGREE SEQUENCE RESOLVE ALL?
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: n=7 — LAMBDA DEGREE SEQUENCE COMPLETENESS")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

by_score_sorted_lam = defaultdict(list)  # (score, sorted_lam) -> list of (H, deg_seq)
count = 0

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    H = count_ham_paths(A, n)

    lam, c3_sets = get_labeled_lambda(A, n)
    sorted_lam = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    deg_seq = lambda_degree_sequence(lam, n)

    by_score_sorted_lam[(scores, sorted_lam)].append({
        'H': H, 'deg_seq': deg_seq, 'bits': bits
    })

    count += 1
    if count % 500000 == 0:
        print(f"  Processed {count}/{total}...")

# Check: for each ambiguous (score, sorted_lam) class, does deg_seq resolve?
print(f"\nTotal (score, sorted_lam) classes: {len(by_score_sorted_lam)}")

resolved_count = 0
unresolved_count = 0
total_ambig_classes = 0

for (sc, sl), group in by_score_sorted_lam.items():
    H_set = sorted(set(d['H'] for d in group))
    if len(H_set) <= 1:
        continue

    total_ambig_classes += 1

    # Check if deg_seq resolves
    by_deg = defaultdict(set)
    for d in group:
        by_deg[d['deg_seq']].add(d['H'])

    all_resolved = all(len(v) == 1 for v in by_deg.values())
    if all_resolved:
        resolved_count += 1
    else:
        unresolved_count += 1
        print(f"  UNRESOLVED: score={sc}, sl=...({len(sl)} vals), H={H_set}")
        for deg, Hs in sorted(by_deg.items()):
            if len(Hs) > 1:
                print(f"    deg_seq={deg}: H={sorted(Hs)}")

print(f"\nAmbiguous sorted-lambda classes: {total_ambig_classes}")
print(f"Resolved by lambda degree sequence: {resolved_count}")
print(f"Unresolved: {unresolved_count}")


# ========================================================================
# ANALYSIS 2: WALK INVARIANTS FOR REMAINING CASES
# ========================================================================
if unresolved_count > 0:
    print(f"\n{'='*70}")
    print("ANALYSIS 2: WALK INVARIANTS FOR UNRESOLVED CASES")
    print("=" * 70)

    # Re-scan only unresolved classes
    for (sc, sl), group in by_score_sorted_lam.items():
        H_set = sorted(set(d['H'] for d in group))
        if len(H_set) <= 1:
            continue

        by_deg = defaultdict(set)
        for d in group:
            by_deg[d['deg_seq']].add(d['H'])

        if all(len(v) == 1 for v in by_deg.values()):
            continue

        # This class is unresolved. Compute walk invariants.
        walk_data = []
        for d in group:
            A_loc = binary_to_tournament(d['bits'], n)
            lam_loc, _ = get_labeled_lambda(A_loc, n)
            scores_loc = [sum(A_loc[v]) for v in range(n)]
            combined = lambda_combined_invariant(lam_loc, n, scores_loc)
            walk_data.append({
                'H': d['H'], 'combined': combined, 'bits': d['bits']
            })

        by_combined = defaultdict(set)
        for wd in walk_data:
            by_combined[wd['combined']].add(wd['H'])

        all_walk_resolved = all(len(v) == 1 for v in by_combined.values())
        print(f"\n  Score {sc}: walk invariants resolve = {all_walk_resolved}")
        if not all_walk_resolved:
            for comb, Hs in sorted(by_combined.items()):
                if len(Hs) > 1:
                    print(f"    combined invariant: H={sorted(Hs)}")
                    # Show one tournament from each H class
                    for h in sorted(Hs):
                        t = next(wd for wd in walk_data if wd['H'] == h and wd['combined'] == comb)
                        print(f"      H={h}, bits={t['bits']}")

else:
    print(f"\n{'='*70}")
    print("ALL AMBIGUITIES RESOLVED BY LAMBDA DEGREE SEQUENCE!")
    print("=" * 70)
    print("""
THEOREM: At n=7, the Hamiltonian path count H(T) is determined by:
  1. The score sequence (sorted out-degrees)
  2. The lambda degree sequence (sorted row sums of pair-coverage matrix)

Equivalently: the lambda pair-coverage matrix, up to graph isomorphism,
determines H(T) for all n <= 7.

This means the labeled lambda graph is a COMPLETE INVARIANT for H
at n <= 7 (and possibly beyond).
""")


# ========================================================================
# ANALYSIS 3: THE {2,1,0} OVERLAP WEIGHT CONNECTION
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE {2,1,0} OVERLAP WEIGHT CONNECTION")
print("=" * 70)

# The overlap weight between two 3-cycles C_i, C_j is:
# w(C_i, C_j) = |V(C_i) cap V(C_j)| in {0, 1, 2}
# (Not 3, since two distinct 3-cycle vertex sets can share at most 2 vertices)

# The lambda value for a vertex pair (u,v) counts how many 3-cycles contain both.
# Connection: w(C_i, C_j) = #{shared vertices} = sum_{v in C_i cap C_j} 1

# The overlap weight between C_i and C_j can be computed from lambda:
# w = sum over vertex-pairs (u,v) in C_i cap C_j of ... wait, w is just |C_i cap C_j|

# The relationship between lambda and overlap:
# lambda_{u,v} = #{C : u in C and v in C} (pair-coverage)
# ov(C_i, C_j) = |C_i cap C_j| (overlap weight)

# These are DUAL perspectives:
# lambda is "how many cycles does this vertex pair see"
# overlap is "how many vertices do these two cycles share"

# The key: ov(C_i, C_j) = 0 iff C_i and C_j are disjoint.
# lambda_{u,v} = 0 iff no cycle contains both u and v.

# For the ambiguous case: same sorted lambda but different disjointness geometry.
# The overlap weight histogram is determined by sorted lambda (proved earlier).
# But the specific ASSIGNMENT of overlap weights to cycle pairs is NOT.

print("""
THE {2,1,0} OVERLAP HIERARCHY:

Level 0: Score sequence
  -> Determines c3 (number of 3-cycle vertex sets)

Level 1: Overlap weight HISTOGRAM (sorted counts of 0,1,2 overlaps)
  -> Same as sorted lambda histogram (proved: ov0+ov1+ov2 = C(c3,2),
     and ov2 = sum C(lam,2), so ov histogram <-> sorted lambda)
  -> Does NOT determine alpha_2 or H in general

Level 2: Overlap weight MATRIX (which cycle pairs have which overlap)
  -> Labeled lambda graph (up to isomorphism)
  -> DETERMINES H at n=5, n=6 (proved), and n=7 (checking)

Level 3: Full tournament adjacency
  -> Contains all information

The {2,1,0} overlap weights between 3-cycle pairs are a TERNARY code.
The sorted histogram is the "measurable" part.
The labeled assignment is the "non-measurable" part.
The non-measurable part determines alpha_2 and hence H.
""")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
