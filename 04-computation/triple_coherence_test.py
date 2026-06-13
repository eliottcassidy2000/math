#!/usr/bin/env python3
"""
triple_coherence_test.py -- kind-pasteur-2026-03-13-S61

The 3-cycle triple coherence sum distinguishes H=109 from H=111
in the deepest ambiguity case. Does it resolve ALL ambiguities at n=7?

Definition: For each triple of 3-cycles (C_i, C_j, C_k) that are pairwise
overlapping (each pair shares at least 1 vertex), compute the product of
their orientations (+1 for clockwise, -1 for counterclockwise w.r.t. canonical
vertex ordering). The sum of these products is the "triple coherence" invariant.

This is computable from the 3-cycle structure alone (no Hamiltonian cycles needed).
Complexity: O(c3^3) where c3 is the number of 3-cycles.

Test: Does (score, sorted_lambda, lambda_deg_seq, walk_invariants, triple_coherence)
resolve ALL H-ambiguities at n=7?

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


def triple_coherence(A, n, c3_sets):
    """Compute the triple coherence invariant.
    For each triple of 3-cycles that are pairwise overlapping,
    compute the product of their orientations."""
    # Get orientations
    orientations = {}
    for c3 in c3_sets:
        a, b, c = sorted(c3)
        if A[a][b] and A[b][c] and A[c][a]:
            orientations[c3] = +1
        else:
            orientations[c3] = -1

    total = 0
    nc = len(c3_sets)
    for i in range(nc):
        for j in range(i+1, nc):
            if len(c3_sets[i] & c3_sets[j]) < 1:
                continue
            for k in range(j+1, nc):
                if (len(c3_sets[j] & c3_sets[k]) >= 1 and
                    len(c3_sets[i] & c3_sets[k]) >= 1):
                    prod = orientations[c3_sets[i]] * orientations[c3_sets[j]] * orientations[c3_sets[k]]
                    total += prod
    return total


def lambda_combined_invariant(lam, n, scores):
    """Combined walk invariant."""
    deg_seq = tuple(sorted(sum(lam[v]) for v in range(n)))

    L2 = [[0]*n for _ in range(n)]
    for v in range(n):
        for u in range(n):
            s = 0
            for w in range(n):
                s += lam[v][w] * lam[w][u]
            L2[v][u] = s

    diag_L2 = tuple(sorted(L2[v][v] for v in range(n)))
    trace_L2 = sum(L2[v][v] for v in range(n))
    row_sums_L2 = tuple(sorted(sum(L2[v]) for v in range(n)))

    lam_degs = [sum(lam[v]) for v in range(n)]
    joint = tuple(sorted(zip(scores, lam_degs)))

    return (deg_seq, diag_L2, trace_L2, row_sums_L2, joint)


# ========================================================================
# FULL n=7 TEST
# ========================================================================
print("=" * 70)
print("FULL n=7 TEST: DOES TRIPLE COHERENCE RESOLVE ALL AMBIGUITIES?")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

# Group by (score, sorted_lambda, combined_walk_invariant)
by_invariant = defaultdict(list)  # key -> list of (H, bits, triple_coh)
count = 0

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    H = count_ham_paths(A, n)

    lam, c3_sets = get_labeled_lambda(A, n)
    sorted_lam = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    combined = lambda_combined_invariant(lam, n, list(scores))

    # Compute triple coherence
    tc = triple_coherence(A, n, c3_sets)

    key = (scores, sorted_lam, combined)
    by_invariant[key].append({'H': H, 'bits': bits, 'tc': tc})

    count += 1
    if count % 500000 == 0:
        print(f"  Processed {count}/{total}...")

print(f"\nTotal invariant classes: {len(by_invariant)}")

# Check for ambiguities
total_ambig = 0
resolved_by_tc = 0
unresolved = 0

for key, group in by_invariant.items():
    H_set = set(d['H'] for d in group)
    if len(H_set) <= 1:
        continue

    total_ambig += 1

    # Check if triple coherence resolves
    by_tc = defaultdict(set)
    for d in group:
        by_tc[d['tc']].add(d['H'])

    if all(len(v) == 1 for v in by_tc.values()):
        resolved_by_tc += 1
    else:
        unresolved += 1
        scores, sl, comb = key
        print(f"\n  UNRESOLVED: score={scores}")
        print(f"    H values: {sorted(H_set)}")
        for tc_val, Hs in sorted(by_tc.items()):
            if len(Hs) > 1:
                print(f"    tc={tc_val}: H={sorted(Hs)}")
                for h in sorted(Hs):
                    example = next(d for d in group if d['H'] == h and d['tc'] == tc_val)
                    print(f"      H={h}, bits={example['bits']}")

print(f"\nAmbiguous walk-invariant classes: {total_ambig}")
print(f"Resolved by triple coherence: {resolved_by_tc}")
print(f"Unresolved: {unresolved}")

if unresolved == 0:
    print(f"\n{'='*70}")
    print("ALL AMBIGUITIES RESOLVED!")
    print("=" * 70)
    print("""
THEOREM: At n=7, H(T) is determined by:
  1. Score sequence
  2. Sorted lambda histogram
  3. Lambda walk invariants (deg_seq, L^2 diagonal/trace/row sums, joint dist)
  4. Triple coherence (product of orientations over overlapping 3-cycle triples)

The triple coherence invariant captures the "non-measurable" content
that the lambda graph cannot see.
""")
else:
    print(f"\n{'='*70}")
    print("TRIPLE COHERENCE IS NOT SUFFICIENT")
    print("=" * 70)
    print(f"""
{unresolved} cases remain unresolved. Additional invariants needed.
The triple coherence captures SOME but not ALL of the hidden dimension.
""")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
