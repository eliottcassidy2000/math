#!/usr/bin/env python3
"""
gk_cluster_corr_s112.py — Cluster correlation structure for the cubic proof
kind-pasteur-2026-03-15-S112

Key finding from gk_tridiag_proof_s112.py:
  E[Z_{j}Z_{j+1} * Z_{k}Z_{k+1}] = 4/(n)_4  for NON-ADJACENT pairs
  while E[Z_j Z_{j+1}]^2 = 4/(n(n-1))^2
  Ratio = n(n-1)/((n-2)(n-3)) --- UNIVERSAL (gap-independent)

Hypothesis: E[prod of r non-adjacent pairs] = 2^r / (n)_{2r}
If true, then the "naive formula" with this weight reproduces g_k exactly,
and the degree-3 cap follows from the combinatorial structure.

Testing at n=8 (3 pairs possible) and verifying the hypothesis.
"""

from itertools import permutations, combinations
from fractions import Fraction
from math import comb, factorial
from collections import defaultdict

def compute_all_pair_expectations(n):
    """For n up to 8, compute E[prod Z_j for j in S] for all domino subsets."""
    m = n - 1
    nfact = factorial(n)

    # Precompute Z for all permutations
    z_all = []
    for perm in permutations(range(n)):
        z = []
        for j in range(m):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z.append(x - y)
        z_all.append(z)

    # Adjacent pairs
    adj_pairs = [(j, j+1) for j in range(m - 1)]

    results = {}

    # Generate all subsets of adj_pairs up to size floor(m/2)
    max_r = min(len(adj_pairs), m // 2)

    for r in range(1, max_r + 1):
        for combo in combinations(range(len(adj_pairs)), r):
            chosen_pairs = [adj_pairs[i] for i in combo]

            # Get all positions
            positions = set()
            for p in chosen_pairs:
                positions.update(p)
            S = frozenset(positions)

            if S in results:
                continue

            # Compute expectation
            total = Fraction(0)
            for z in z_all:
                prod = 1
                for j in S:
                    prod *= z[j]
                total += prod

            results[S] = total / nfact

    return results

def classify_by_pairs(S):
    """Given a domino subset, find the pair decomposition and gaps."""
    positions = sorted(S)
    # Find contiguous clusters
    clusters = []
    current = [positions[0]]
    for p in positions[1:]:
        if p == current[-1] + 1:
            current.append(p)
        else:
            clusters.append(current)
            current = [p]
    clusters.append(current)

    cluster_sizes = tuple(len(c) for c in clusters)

    # Compute gaps between clusters
    gaps = []
    for i in range(len(clusters) - 1):
        gap = clusters[i+1][0] - clusters[i][-1] - 1
        gaps.append(gap)

    # Number of pairs in each cluster
    pairs_per_cluster = tuple((len(c) + 1) // 2 for c in clusters)

    return cluster_sizes, gaps, pairs_per_cluster

def falling_factorial(n, k):
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

print("="*70)
print("CLUSTER CORRELATION ANALYSIS")
print("="*70)

for n in range(3, 9):
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")

    results = compute_all_pair_expectations(n)

    # Group by cluster structure
    by_structure = defaultdict(list)
    for S, val in results.items():
        if val == 0:
            continue
        cs, gaps, ppc = classify_by_pairs(S)
        by_structure[(cs, tuple(gaps))].append((S, val))

    for key in sorted(by_structure.keys()):
        cs, gaps = key
        items = by_structure[key]
        # All should have same expectation (translation invariance of the formula)
        vals = [v for _, v in items]
        all_same = all(v == vals[0] for v in vals)

        # Count pairs
        total_pairs = sum((s + 1) // 2 for s in cs)  # approximate
        # Actually, need exact pair count
        # A cluster of size L covers L-1 adjacent pair slots, so contributes ceil(L/2) pairs
        # But that's not right either. Actually the pair count = |S|//2 for even |S|

        total_size = sum(cs)
        k = total_size // 2  # number of pairs

        # Prediction: 2^r / (n)_{2r} where r = number of INDEPENDENT pairs
        # For separated clusters: r = number of clusters (each contributes 1 pair? No)
        # Let me just check against 2^k / (n)_{2k}
        pred_simple = Fraction(2**k, falling_factorial(n, 2*k)) if 2*k <= n else None

        print(f"\n  cluster sizes={cs}, gaps={gaps}: {len(items)} subsets, "
              f"all same={all_same}, val={vals[0]}")
        if pred_simple is not None:
            ratio = vals[0] / pred_simple
            print(f"    vs 2^{k}/(n)_{{{2*k}}} = {pred_simple}, ratio = {ratio}")

# Now specifically test the hypothesis for separated pairs
print("\n" + "="*70)
print("SEPARATED PAIRS HYPOTHESIS: E[prod] = 2^r / (n)_{2r}")
print("="*70)

for n in range(5, 9):
    print(f"\nn={n}:")
    results = compute_all_pair_expectations(n)

    for S, val in sorted(results.items(), key=lambda x: (len(x[0]), x[0])):
        if val == 0:
            continue
        cs, gaps, ppc = classify_by_pairs(S)

        # Only look at subsets where ALL clusters have size 2 (isolated pairs)
        if any(c != 2 for c in cs):
            continue

        r = len(cs)  # number of separated pairs
        pred = Fraction(2**r, falling_factorial(n, 2*r))
        match = "OK" if val == pred else f"MISMATCH (ratio={val/pred})"

        if r >= 2:  # Skip single pairs (trivially true)
            print(f"  S={sorted(S)}, {r} sep pairs, gaps={gaps}: E={val}, pred={pred} {match}")

# Connected cumulant analysis
print("\n" + "="*70)
print("CONNECTED CUMULANTS")
print("="*70)
print("If separated pairs factorize as 2^r/(n)_{2r},")
print("then connected cumulants come only from OVERLAPPING clusters.")

for n in range(5, 9):
    print(f"\nn={n}:")
    results = compute_all_pair_expectations(n)

    for S, val in sorted(results.items(), key=lambda x: (len(x[0]), x[0])):
        if val == 0:
            continue
        cs, gaps, ppc = classify_by_pairs(S)

        # Look at CONNECTED subsets (single cluster)
        if len(cs) != 1:
            continue

        L = cs[0]  # cluster size
        k = L // 2  # number of "pair slots" covered by this cluster... actually |S|=L, pairs=L-1
        # The cluster covers positions j, j+1, ..., j+L-1
        # This is L positions, containing L-1 adjacent pairs
        # But as a domino subset, it's a union of SOME pairs, not necessarily all

        # Actually, for a contiguous block of size L to be a domino subset,
        # every position must be in some pair. So L must be even.
        if L % 2 == 1:
            continue

        r = L // 2  # number of pairs
        pred_sep = Fraction(2**r, falling_factorial(n, 2*r))
        connected_excess = val - pred_sep
        ratio = val / pred_sep if pred_sep != 0 else "undef"

        if L >= 4:
            print(f"  contiguous L={L}: E={val}, sep_pred={pred_sep}, "
                  f"ratio={ratio}, excess={connected_excess}")

# The key question: for a contiguous block of 4 (2 overlapping pairs),
# what is E[Z_j Z_{j+1} Z_{j+2} Z_{j+3}]?
print("\n" + "="*70)
print("CONTIGUOUS BLOCK EXPECTATIONS")
print("="*70)

for n in range(4, 9):
    results = compute_all_pair_expectations(n)

    print(f"\nn={n}:")

    for L in range(2, 8, 2):  # even block sizes
        # Find a contiguous block of size L
        S = frozenset(range(L))
        if S in results and results[S] != 0:
            val = results[S]
            r = L // 2
            sep_pred = Fraction(2**r, falling_factorial(n, 2*r))
            ratio_to_sep = val / sep_pred if sep_pred != 0 else "?"

            # Also compare to 2/(n)_L
            cluster_pred = Fraction(2, falling_factorial(n, L))
            ratio_to_cluster = val / cluster_pred if cluster_pred != 0 else "?"

            print(f"  block L={L}: E={val}")
            print(f"    vs 2^{r}/(n)_{{{2*r}}} = {sep_pred} (ratio {ratio_to_sep})")
            print(f"    vs 2/(n)_{{{L}}} = {cluster_pred} (ratio {ratio_to_cluster})")

print("\nDone!")
