#!/usr/bin/env python3
"""
gk_structure_89c.py — Investigate WHY g_k(m) is always degree 3 for k≥3
opus-2026-03-15-S89c

Hypothesis: g_k(m) counts weighted domino tilings. Each tiling has k dominos
on positions 0..m+2k-2. The weight depends on the cluster structure (how
dominos group into consecutive runs).

Let's decompose g_k(m) by cluster type and see if this explains degree 3.
"""

from fractions import Fraction
from itertools import combinations
from math import factorial, comb
from functools import reduce

def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)

def domino_tilings(m, k):
    """Generate all ways to place k non-overlapping dominos on positions 0..m+2k-2.
    Each domino covers {j, j+1} for some j.
    Returns list of tuples of domino start positions.
    """
    n_pos = m + 2*k - 1  # total positions
    tilings = []

    def backtrack(start, remaining, current):
        if remaining == 0:
            tilings.append(tuple(current))
            return
        for j in range(start, n_pos - 1):
            # Place domino at {j, j+1}
            # Make sure it doesn't overlap previous
            if current and j <= current[-1] + 1:
                continue
            backtrack(j + 2, remaining - 1, current + [j])

    backtrack(0, k, [])
    return tilings

def cluster_type(tiling):
    """Decompose a tiling into clusters of consecutive dominos.
    Returns sorted tuple of cluster sizes.
    """
    if not tiling:
        return ()
    clusters = []
    current_cluster = 1
    for i in range(1, len(tiling)):
        if tiling[i] == tiling[i-1] + 2:  # Adjacent dominos (consecutive positions)
            current_cluster += 1
        else:
            clusters.append(current_cluster)
            current_cluster = 1
    clusters.append(current_cluster)
    return tuple(sorted(clusters, reverse=True))

# Count tilings by cluster type for various k and m
print("="*70)
print("DOMINO TILINGS BY CLUSTER TYPE")
print("="*70)

for k in range(1, 7):
    print(f"\nk={k}:")
    for m in range(1, 8):
        tilings = domino_tilings(m, k)
        type_counts = {}
        for t in tilings:
            ct = cluster_type(t)
            type_counts[ct] = type_counts.get(ct, 0) + 1

        total = len(tilings)
        # Expected: C(m+k-1, k) = C(m-1+k, k) for non-overlapping placements on m+2k-1 positions
        expected = comb(m + k - 1, k)
        print(f"  m={m}: total tilings={total} (expected C({m+k-1},{k})={expected})", end="")
        if total != expected:
            print(" MISMATCH!", end="")
        print()

        for ct in sorted(type_counts):
            print(f"    cluster type {ct}: {type_counts[ct]}")

# The number of domino tilings with cluster type λ = (λ_1, ..., λ_r)
# depends on m and the λ_i. Each cluster of size s covers 2s consecutive positions.
# Between clusters, there must be a gap of at least 1.
# After accounting for the clusters (which take 2Σλ_i = 2k positions) and gaps (r-1),
# we have m-1 remaining positions to distribute as gaps.
# The number of ways = C(m - 1 + r - 1, r) = C(m + r - 2, r) for a specific ORDERED partition,
# but we need to account for permutations of identical clusters.

print("\n" + "="*70)
print("DECOMPOSITION: tilings by number of clusters r")
print("="*70)

for k in range(1, 7):
    print(f"\nk={k}:")
    for m in range(1, 8):
        tilings = domino_tilings(m, k)
        by_r = {}
        for t in tilings:
            ct = cluster_type(t)
            r = len(ct)
            by_r[r] = by_r.get(r, 0) + 1
        print(f"  m={m}: by #clusters: {dict(sorted(by_r.items()))}")

# For k dominos with r clusters (ordered partitions of k into r parts):
# - The number of compositions of k into r parts = C(k-1, r-1)
# - Each composition λ = (λ_1, ..., λ_r) has clusters using 2k positions
# - Gaps: need r-1 mandatory gaps (at least 1 each) plus distribute m-1 remaining
# - Wait: total positions = m + 2k - 1. Clusters use 2k. Gaps use at least r-1.
#   Free positions = (m + 2k - 1) - 2k - (r-1) = m - r.
#   Distribute m-r into r+1 bins (before first, between, after last): C(m-r+r, r) = C(m, r).
#   But wait, r+1 bins with m-r stars: stars-and-bars = C(m-r + r, r) = C(m, r).

# So the number of tilings with r clusters of specified sizes = C(m, r)
# But we need to count compositions, not partitions:
# - Number of compositions of k into r parts = C(k-1, r-1)
# - For EACH composition, the number of placements = C(m, r)
# - Total tilings with r clusters = C(k-1, r-1) · C(m, r)

# Verify:
print("\n" + "="*70)
print("VERIFY: tilings with r clusters = C(k-1,r-1)·C(m,r)")
print("="*70)

for k in range(1, 7):
    print(f"\nk={k}:")
    all_match = True
    for m in range(1, 8):
        tilings = domino_tilings(m, k)
        by_r = {}
        for t in tilings:
            ct = cluster_type(t)
            r = len(ct)
            by_r[r] = by_r.get(r, 0) + 1

        for r in range(1, k+1):
            actual = by_r.get(r, 0)
            predicted = comb(k-1, r-1) * comb(m, r)
            if actual != predicted:
                print(f"  m={m}, r={r}: actual={actual}, predicted={predicted} MISMATCH!")
                all_match = False
    if all_match:
        print(f"  ✓ All match!")

# Total tilings: Σ_r C(k-1,r-1)·C(m,r) = C(m+k-1,k) (Vandermonde convolution)
# Verification is already in the counts above.

# So the total tiling count is C(m+k-1, k), and the r-cluster count is C(k-1,r-1)·C(m,r).
# Now the WEIGHT of each tiling depends on its cluster structure.
# The weight w(tiling) = E[∏_{j∈tiling} Z_j] for a permutation of n = m+2k elements.
# This is a product that depends on the cluster structure.

# For EACH cluster of size s, the contribution is a function of s and n.
# For isolated dominos (s=1), the contribution is E[Z_j Z_{j+1}] = 2/(n(n-1)).
# For a cluster of size 2 (4 consecutive), the contribution is E[Z_j Z_{j+1} Z_{j+2} Z_{j+3}].
# Etc.

# Let h_s(n) = E[∏_{j=0}^{2s-1} Z_j] for a cluster of s adjacent dominos in n elements.
# Then for a tiling with clusters of sizes s_1, ..., s_r at positions that are
# far apart, the expectation should factorize:
# E[∏Z] ≈ ∏ h_{s_i}(n).
# But they're NOT independent because the positions share the same permutation.
# However, for dominos at positions that are far apart, the correlations are 0
# (tridiagonal structure), so the factorization IS exact for non-adjacent clusters.

# So: g_k(m) = Σ_tilings w(tiling) × (n)_{2k} / 2
#            = Σ_r Σ_{compositions (s_1,...,s_r) of k} C(m,r) × ∏ h_{s_i}(n) × (n)_{2k} / 2

# Wait, but n = m + 2k. So (n)_{2k} depends on m. And h_{s_i}(n) depends on n = m + 2k.
# Let me think about this more carefully.

# Actually, the weight of a domino tiling is E[∏_{j∈S} Z_j] where S is the set of
# 2k positions covered, and the expectation is over Perm(n) with n = m + 2k.
# For non-overlapping clusters, this factorizes.

# For a single cluster of s dominos covering positions p, p+1, ..., p+2s-1:
# E[Z_p · Z_{p+1} · ... · Z_{p+2s-1}] = h_s(n)
# This is independent of the position p (by translation invariance of uniform perms).
# Actually, IS it independent of p? For uniform random permutations of {0,...,n-1},
# the statistics X_j, Y_j depend on the permutation values, not positions.
# By exchangeability, E[Z_p·...·Z_{p+2s-1}] should be the same for all p (as long as
# all indices are valid). Let me verify this.

# Actually, E[Z_j Z_{j+1}] = 2/(n(n-1)) is the same for all j, which supports this.

# So the weight of a tiling with cluster sizes s_1, ..., s_r is ∏_i h_{s_i}(n).
# And g_k(m) · 2/(n)_{2k} = Σ_tilings ∏ h_{s_i}(n)
#                          = Σ_r C(m,r) · Σ_{compositions of k into r parts} ∏ h_{s_i}(n)
# But wait, C(m,r) counts the number of ways to place r clusters with gaps,
# and each composition (s_1,...,s_r) of k contributes ∏ h_{s_i}(n).

# g_k(m) = (n)_{2k}/2 · Σ_r C(m,r) · Σ_{compositions of k into r parts} ∏ h_{s_i}(n)

# Define p_r(k) = Σ_{compositions of k into r parts} ∏ h_{s_i}(n)
# Then g_k(m) = (n)_{2k}/2 · Σ_r C(m,r) · p_r(k)

# Since C(m,r) is a polynomial of degree r in m, and r ≤ k, the degree of
# g_k(m) is at most k. But we see degree 3 for all k ≥ 3. This means
# p_r(k) = 0 for all r ≥ 4!

# Let me check: for k ≥ 4, do terms with r ≥ 4 clusters contribute?
# p_4(k) = Σ_{s_1+...+s_4=k, all s_i≥1} h_{s_1}·h_{s_2}·h_{s_3}·h_{s_4}
# This would need 4 clusters, each of size ≥ 1. If all h_s are nonzero, this should be nonzero.

# Unless the normalization by (n)_{2k}/2 absorbs the high-degree terms.
# Wait, n = m + 2k, so (n)_{2k} = (m+2k)(m+2k-1)...(m+1). This is degree 2k in m!
# So g_k(m) = (degree 2k in m) × Σ_r C(m,r) × (function of n=m+2k)
# The h_{s_i}(n) are functions of n = m+2k, which is itself a function of m.

# This is why it's not obvious that the degree is 3 — there's massive cancellation.

# Let me just compute h_s(n) explicitly for small s.
print("\n" + "="*70)
print("CLUSTER WEIGHTS h_s(n) = E[Z_0·Z_1·...·Z_{2s-1}] (from brute force)")
print("="*70)

from itertools import permutations as perms

for n in range(3, 9):
    print(f"\nn={n}:")
    for s in range(1, n//2 + 1):
        size = 2 * s
        if size > n - 1:
            break
        # Compute E[Z_0·Z_1·...·Z_{2s-1}]
        total = 0
        N = factorial(n)
        for perm in perms(range(n)):
            prod = 1
            for j in range(size):
                xj = 1 if perm[j+1] == perm[j] + 1 else 0
                yj = 1 if perm[j+1] == perm[j] - 1 else 0
                prod *= (xj - yj)
            total += prod
        h = Fraction(total, N)
        print(f"  h_{s}({n}) = {h} = {float(h):.10f}")
        # Compare with 2/(n)_{2s}
        ff = falling_factorial(n, 2*s)
        ratio = h * ff / 2
        print(f"    h_{s}·(n)_{{{2*s}}}/2 = {ratio}")

print("\nDone!")
