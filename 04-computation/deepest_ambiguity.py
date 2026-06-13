#!/usr/bin/env python3
"""
deepest_ambiguity.py -- kind-pasteur-2026-03-13-S61

Investigate the 2 hardest cases where even walk invariants of
the labeled lambda graph don't determine H:

1. Score (2,2,3,3,3,3,5): bits 4728 (H=109) vs bits 4658 (H=111)
2. Score (1,3,3,3,3,4,4): bits 9388 (H=109) vs bits 9653 (H=111)

Both have H difference of 2 = one more directed cycle.
Question: are their lambda graphs isomorphic?

If YES: labeled lambda does NOT determine H at n=7.
If NO: our walk invariants just weren't fine enough.

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


def are_lambda_isomorphic(lam1, lam2, n):
    """Check if two lambda matrices are isomorphic (related by vertex permutation).
    For n=7 this checks 7! = 5040 permutations."""
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(i+1, n):
                if lam1[perm[i]][perm[j]] != lam2[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True, perm
    return False, None


n = 7

# ========================================================================
# CASE 1: Score (2,2,3,3,3,3,5): bits 4728 (H=109) vs bits 4658 (H=111)
# ========================================================================
print("=" * 70)
print("CASE 1: Score (2,2,3,3,3,3,5)")
print("=" * 70)

for bits, expected_H in [(4728, 109), (4658, 111)]:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = [sum(A[v]) for v in range(n)]
    lam, c3_sets = get_labeled_lambda(A, n)

    print(f"\n  bits={bits}, H={H} (expected {expected_H})")
    print(f"  Out-degrees: {scores}")
    print(f"  3-cycle vertex sets: {[sorted(fs) for fs in c3_sets]}")

    # Lambda matrix
    print(f"  Lambda matrix:")
    for v in range(n):
        row = [lam[v][w] for w in range(n)]
        print(f"    {row}")

    # Cycle counts
    for k in [3, 5, 7]:
        ck = 0
        for subset in combinations(range(n), k):
            ck += count_directed_ham_cycles_on_subset(A, list(subset))
        print(f"  c{k}_dir = {ck}")

    # Alpha_2
    all_cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            nc = count_directed_ham_cycles_on_subset(A, list(subset))
            for _ in range(nc):
                all_cycles.append(frozenset(subset))

    alpha2 = 0
    nc = len(all_cycles)
    for i in range(nc):
        for j in range(i+1, nc):
            if not (all_cycles[i] & all_cycles[j]):
                alpha2 += 1
    print(f"  Total cycles: {nc}, alpha_2: {alpha2}")

# Lambda isomorphism check
A1 = binary_to_tournament(4728, n)
A2 = binary_to_tournament(4658, n)
lam1, _ = get_labeled_lambda(A1, n)
lam2, _ = get_labeled_lambda(A2, n)

print(f"\nLambda isomorphism check...")
iso, perm = are_lambda_isomorphic(lam1, lam2, n)
if iso:
    print(f"  ISOMORPHIC! Permutation: {perm}")
    print(f"  => LABELED LAMBDA DOES NOT DETERMINE H AT n=7!")
else:
    print(f"  NOT isomorphic!")
    print(f"  => Walk invariants were too weak. Labeled lambda IS finer.")


# ========================================================================
# CASE 2: Score (1,3,3,3,3,4,4): bits 9388 (H=109) vs bits 9653 (H=111)
# ========================================================================
print(f"\n{'='*70}")
print("CASE 2: Score (1,3,3,3,3,4,4)")
print("=" * 70)

for bits, expected_H in [(9388, 109), (9653, 111)]:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = [sum(A[v]) for v in range(n)]
    lam, c3_sets = get_labeled_lambda(A, n)

    print(f"\n  bits={bits}, H={H} (expected {expected_H})")
    print(f"  Out-degrees: {scores}")
    print(f"  3-cycle vertex sets: {[sorted(fs) for fs in c3_sets]}")

    # Lambda matrix
    print(f"  Lambda matrix:")
    for v in range(n):
        row = [lam[v][w] for w in range(n)]
        print(f"    {row}")

    # Cycle counts
    for k in [3, 5, 7]:
        ck = 0
        for subset in combinations(range(n), k):
            ck += count_directed_ham_cycles_on_subset(A, list(subset))
        print(f"  c{k}_dir = {ck}")

    # Alpha_2
    all_cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            nc = count_directed_ham_cycles_on_subset(A, list(subset))
            for _ in range(nc):
                all_cycles.append(frozenset(subset))

    alpha2 = 0
    nc = len(all_cycles)
    for i in range(nc):
        for j in range(i+1, nc):
            if not (all_cycles[i] & all_cycles[j]):
                alpha2 += 1
    print(f"  Total cycles: {nc}, alpha_2: {alpha2}")

# Lambda isomorphism check
A1 = binary_to_tournament(9388, n)
A2 = binary_to_tournament(9653, n)
lam1, _ = get_labeled_lambda(A1, n)
lam2, _ = get_labeled_lambda(A2, n)

print(f"\nLambda isomorphism check...")
iso, perm = are_lambda_isomorphic(lam1, lam2, n)
if iso:
    print(f"  ISOMORPHIC! Permutation: {perm}")
    print(f"  => LABELED LAMBDA DOES NOT DETERMINE H AT n=7!")
else:
    print(f"  NOT isomorphic!")
    print(f"  => Walk invariants were too weak. Labeled lambda IS finer.")


# ========================================================================
# CONCLUSION
# ========================================================================
print(f"\n{'='*70}")
print("CONCLUSION")
print("=" * 70)

if not iso:
    print("""
The labeled lambda graph distinguishes all tested ambiguous cases.
Walk invariants (degree sequence + L^2 trace + row sums + joint distribution)
were insufficient because they are LOW-ORDER approximations of the full
graph isomorphism class.

Current status:
  - Labeled lambda determines H at n=5 (PROVED, exhaustive canonical form)
  - Labeled lambda determines H at n=6 (PROVED, exhaustive canonical form)
  - Labeled lambda LIKELY determines H at n=7 (all specific counterexample
    candidates are NOT lambda-isomorphic)

This suggests:
  CONJECTURE: The labeled lambda graph (pair-coverage matrix up to
  permutation) is a COMPLETE INVARIANT for H(T) at all n.
""")
else:
    print("""
NEGATIVE RESULT: The labeled lambda graph does NOT determine H at n=7.
Two tournaments with isomorphic lambda graphs have different H values.

This means the pair-coverage of 3-cycles (even as a labeled graph)
is NOT sufficient to determine H. Additional cycle information (5-cycle
orientation, Hamiltonian cycle count, or higher-order overlap invariants)
is genuinely independent of the 3-cycle pair-coverage structure.

The Vitali analogy is CONFIRMED at the labeled level:
  The full tournament contains "non-measurable" information
  even beyond the pair-coverage graph.
""")

print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
