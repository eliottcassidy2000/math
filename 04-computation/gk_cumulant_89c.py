#!/usr/bin/env python3
"""
gk_cumulant_89c.py — Analyze the cumulant/correlation structure behind g_k
opus-2026-03-15-S89c

KEY DISCOVERY: In the binomial basis C(m,r), the true g_k(m) has support
ONLY on r=0,1,2,3. The r≥4 coefficients are exactly zero.

This means: the weighted sum over all r-cluster tilings is zero for r≥4.
But each individual tiling has nonzero weight!
So there must be exact cancellation — which is a CUMULANT phenomenon.

The naive formula assumes clusters are independent:
  weight(tiling with clusters s_1,...,s_r) = ∏ h_{s_i}(n)
  where h_s(n) = 2/(n)_{2s}.

The true weight involves joint expectations that don't factorize.
The difference (true - independent) is captured by cumulants.

For r clusters, the joint expectation E[C_1·...·C_r] can be expanded via
the moment-cumulant formula:
  E[C_1·...·C_r] = Σ_partitions ∏_blocks κ(C_{i_1},...,C_{i_j})

The cumulant κ_r(C_1,...,C_r) involves inclusion-exclusion over partitions.
For r≥4, the claim is that the sum over all compositions of k into r parts
of these cumulant contributions is exactly -∏h_{s_i}, canceling the independent part.

Let's investigate this by computing ACTUAL joint expectations for multi-cluster
configurations and comparing with the independent model.
"""

from fractions import Fraction
from itertools import permutations
from math import factorial, comb

def compute_expectation(n, positions_list):
    """Compute E[∏_clusters Z_cluster] where each cluster is a list of positions.
    Z_cluster = ∏_{j in cluster} Z_j where Z_j = X_j - Y_j.
    Returns exact fraction.
    """
    N = factorial(n)
    total = 0
    all_positions = []
    for cluster in positions_list:
        all_positions.extend(cluster)

    for perm in permutations(range(n)):
        prod = 1
        for j in all_positions:
            if j + 1 >= n:
                prod = 0
                break
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            prod *= (xj - yj)
        total += prod
    return Fraction(total, N)

def h_s(n, s):
    """Single cluster weight: h_s(n) = 2/(n)_{2s}."""
    ff = 1
    for i in range(2*s):
        ff *= (n - i)
    return Fraction(2, ff)

# First: compute exact joint expectations for small cases
# and compare with independent model
print("="*70)
print("JOINT vs INDEPENDENT CLUSTER EXPECTATIONS")
print("="*70)

for n in range(6, 10):
    print(f"\nn={n}:")
    h1 = h_s(n, 1)  # = 2/(n(n-1))
    print(f"  h_1 = {h1}")

    # 2 clusters of size 1, gap g: positions {0,1} and {2+g, 3+g}
    for g in range(0, min(4, n-5)):
        pos1 = [0, 1]
        pos2 = [2+g, 3+g]
        joint = compute_expectation(n, [pos1, pos2])
        independent = h1 * h1
        ratio = joint / independent if independent != 0 else None
        cumulant = joint - independent  # κ_2 for two clusters

        print(f"  2 clusters, gap={g}: joint={joint}, indep={independent}")
        print(f"    ratio={ratio}, cumulant κ₂={cumulant}")

    # 3 clusters of size 1
    for g1 in range(0, min(3, n-7)):
        for g2 in range(0, min(3, n-7-g1)):
            pos1 = [0, 1]
            p2_start = 2 + g1
            pos2 = [p2_start, p2_start + 1]
            p3_start = p2_start + 2 + g2
            pos3 = [p3_start, p3_start + 1]
            if p3_start + 1 >= n:
                continue

            joint = compute_expectation(n, [pos1, pos2, pos3])
            independent = h1 ** 3

            # 3-point cumulant: κ₃ = joint - Σ pairs κ₂·h₁ - h₁³
            pair12 = compute_expectation(n, [pos1, pos2])
            pair13 = compute_expectation(n, [pos1, pos3])
            pair23 = compute_expectation(n, [pos2, pos3])

            # Moment-cumulant: E[123] = κ₃ + κ₂(12)·κ₁(3) + κ₂(13)·κ₁(2) + κ₂(23)·κ₁(1) + κ₁(1)·κ₁(2)·κ₁(3)
            # where κ₁(i) = h₁, κ₂(ij) = pair_ij - h₁²
            k2_12 = pair12 - h1**2
            k2_13 = pair13 - h1**2
            k2_23 = pair23 - h1**2

            kappa3 = joint - k2_12*h1 - k2_13*h1 - k2_23*h1 - h1**3
            print(f"  3 clusters, gaps=({g1},{g2}): joint={joint}")
            print(f"    κ₃={kappa3}")

# Now: what IS the effective weight for r-cluster tilings?
# The sum over all compositions of k into r parts of the true joint weight
# gives the C(m,r) coefficient in g_k(m).
# For the true g_k, this is 0 for r≥4.

print("\n" + "="*70)
print("EFFECTIVE r-CLUSTER WEIGHTS (summed over compositions)")
print("="*70)

# For small n, compute g_k(m) directly and extract the effective weights
# The effective weight W_r(k) satisfies:
#   g_k(m) = Σ_r W_r(k) · C(m,r)
# So W_r(k) = Δ^r g_k(0) (forward differences)

gk_coeffs = {
    1: (0, 0, 3, 0),
    2: (0, 3, 0, 0),
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}

def g_true(k, m):
    if k == 1: return m
    if k == 2: return m*m
    if k in gk_coeffs:
        a, b, c, d = gk_coeffs[k]
        return Fraction(a*m**3 + b*m**2 + c*m + d, 3)
    return None

for k in range(1, 10):
    vals = [g_true(k, m) for m in range(k+4)]
    # Forward differences
    curr = list(vals)
    deltas = [curr[0]]
    for r in range(1, len(curr)):
        curr = [curr[i+1] - curr[i] for i in range(len(curr)-1)]
        if curr:
            deltas.append(curr[0])

    # Trim trailing zeros
    while deltas and deltas[-1] == 0:
        deltas.pop()

    # Naive weights
    naive_w = [comb(k-1, r-1) * 2**(r-1) if r >= 1 else 0 for r in range(k+1)]

    print(f"\nk={k}:")
    print(f"  True Δ^r: {deltas}")
    print(f"  Naive:    {naive_w[:len(deltas)+2]}")
    if len(deltas) <= 4:
        print(f"  DEGREE ≤ 3 ✓")

# Now: the key question is what the r=3 coefficient represents.
# For r=3 clusters of sizes s_1, s_2, s_3 with s_1+s_2+s_3=k:
# The effective weight includes cumulant corrections.
# The independent contribution would be:
#   Σ_{compositions} ∏ h_{s_i} · (n)_{2k}/2 = C(k-1,2) · (some product average)
# But the naive formula uses 2^{r-1} = 4 per 3-cluster configuration.

# Let's compute the ACTUAL per-composition weights for k=4
print("\n" + "="*70)
print("PER-COMPOSITION WEIGHTS for k=4")
print("="*70)

k = 4
for n in range(9, 12):
    m = n - 2*k
    if m < 1:
        continue
    print(f"\nn={n}, m={m}:")
    ff_2k = 1
    for i in range(2*k):
        ff_2k *= (n - i)

    # Compositions of 4 into r parts
    # r=1: (4,) — single cluster of size 4
    # r=2: (1,3), (2,2), (3,1)
    # r=3: (1,1,2), (1,2,1), (2,1,1)
    # r=4: (1,1,1,1)

    # For each composition, the weight is E[cluster1·cluster2·...] · (n)_{2k}/2
    # where clusters are at specific positions in a tiling.

    # For r=4 (all singletons): positions {0,1}, {2,3}, {4,5}, {6,7}
    # These are at fixed positions — but we need to sum over ALL placements.
    # Actually, C(m,4) gives the number of ways to place 4 singleton clusters
    # with gaps ≥ 1 between them.

    # The weight per placement should be h_1^4 · (n)_{2k}/2 times correlation correction.
    # But since all gaps ≥ 1 (actually ≥ 2 in position space):
    # Let me just compute for specific positions.

    # Place 4 dominos at positions 0,1 | 3,4 | 6,7 | 9,10 (gap=1 each)
    if n > 10:
        pos_list = [[0,1], [3,4], [6,7], [9,10]]
        joint = compute_expectation(n, pos_list)
        indep = h_s(n,1)**4
        print(f"  4 singletons at 01|34|67|9,10: joint={joint}")
        print(f"    independent = {indep}")
        print(f"    ratio = {joint/indep}")
        print(f"    (n)_8/2 = {ff_2k//2}")
        print(f"    weighted = {joint * ff_2k // 2}")

# For r=3 with k=4: compositions are (1,1,2), (1,2,1), (2,1,1)
# Total C(3,2) = 3 compositions, each giving C(m,3) placements
print("\n" + "="*70)
print("UNDERSTANDING THE r=3 WEIGHT FOR k=4")
print("="*70)

for n in [9, 10, 11]:
    m = n - 8
    if m < 3:
        continue
    print(f"\nn={n}, m={m}:")
    ff_8 = 1
    for i in range(8):
        ff_8 *= (n - i)

    # Composition (1,1,2): clusters of size 1,1,2
    # Place at positions 0,1 | 3,4 | 6,7,8,9 (size-2 cluster at end)
    if n >= 10:
        pos = [[0,1], [3,4], [6,7,8,9]]
        joint = compute_expectation(n, pos)
        indep = h_s(n,1) * h_s(n,1) * h_s(n,2)
        print(f"  (1,1,2) at 01|34|6789: joint={joint}")
        print(f"    independent = {indep}")
        print(f"    ratio = {joint/indep}")
        print(f"    weighted joint = {joint * ff_8 / 2}")

    # Composition (1,2,1)
    if n >= 10:
        pos = [[0,1], [3,4,5,6], [8,9]]
        joint = compute_expectation(n, pos)
        indep = h_s(n,1) * h_s(n,2) * h_s(n,1)
        print(f"  (1,2,1) at 01|3456|89: joint={joint}")
        print(f"    independent = {indep}")
        print(f"    ratio = {joint/indep}")

    # Composition (2,1,1)
    if n >= 10:
        pos = [[0,1,2,3], [5,6], [8,9]]
        joint = compute_expectation(n, pos)
        indep = h_s(n,2) * h_s(n,1) * h_s(n,1)
        print(f"  (2,1,1) at 0123|56|89: joint={joint}")
        print(f"    independent = {indep}")
        print(f"    ratio = {joint/indep}")

# KEY TEST: is the ratio independent of position/gaps?
print("\n" + "="*70)
print("GAP INDEPENDENCE: does ratio depend on gap sizes?")
print("="*70)

n = 10
h1 = h_s(n, 1)
h2 = h_s(n, 2)
print(f"n={n}, h_1={h1}, h_2={h2}")

# Two size-1 clusters with varying gaps
print("\nTwo s=1 clusters:")
for g in range(0, 5):
    p1 = [0, 1]
    p2 = [2+g, 3+g]
    if p2[-1] >= n:
        break
    joint = compute_expectation(n, [p1, p2])
    print(f"  gap={g}: ratio = {joint / (h1*h1)}")

# s=1 and s=2 clusters with varying gaps
print("\ns=1 and s=2 clusters:")
for g in range(0, 4):
    p1 = [0, 1]
    p2 = [2+g, 3+g, 4+g, 5+g]
    if p2[-1] >= n:
        break
    joint = compute_expectation(n, [p1, p2])
    print(f"  gap={g}: ratio = {joint / (h1*h2)}")

# Two s=2 clusters with varying gaps
print("\nTwo s=2 clusters:")
for g in range(0, 3):
    p1 = [0, 1, 2, 3]
    p2 = [4+g, 5+g, 6+g, 7+g]
    if p2[-1] >= n:
        break
    joint = compute_expectation(n, [p1, p2])
    print(f"  gap={g}: ratio = {joint / (h2*h2)}")

# Three s=1 clusters with varying gaps
print("\nThree s=1 clusters (equal gaps):")
for g in range(0, 3):
    p1 = [0, 1]
    p2 = [2+g, 3+g]
    p3 = [4+2*g, 5+2*g]
    if p3[-1] >= n:
        break
    joint = compute_expectation(n, [p1, p2, p3])
    print(f"  gaps=({g},{g}): ratio = {joint / (h1**3)}")

print("\nThree s=1 clusters (different gaps):")
for g1 in range(0, 3):
    for g2 in range(0, 3):
        p1 = [0, 1]
        p2 = [2+g1, 3+g1]
        p3 = [4+g1+g2, 5+g1+g2]
        if p3[-1] >= n:
            continue
        joint = compute_expectation(n, [p1, p2, p3])
        print(f"  gaps=({g1},{g2}): ratio = {joint / (h1**3)}")

print("\nDone!")
