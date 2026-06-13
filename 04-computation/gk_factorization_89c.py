#!/usr/bin/env python3
"""
gk_factorization_89c.py — Test whether non-adjacent domino cluster expectations factorize
opus-2026-03-15-S89c

Key question: does E[Z_0Z_1 · Z_lZ_{l+1}] = E[Z_0Z_1] · E[Z_lZ_{l+1}] for l≥3?

If yes, the g_k structure has a simple product decomposition.
If no, the corrections determine the degree-3 polynomial structure.
"""

from fractions import Fraction
from itertools import permutations
from math import factorial

def compute_joint(n, positions):
    """Compute E[∏_{j∈positions} Z_j] exactly."""
    N = factorial(n)
    total = 0
    for perm in permutations(range(n)):
        prod = 1
        for j in positions:
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            prod *= (xj - yj)
        total += prod
    return Fraction(total, N)

print("="*70)
print("TEST: FACTORIZATION OF NON-ADJACENT CLUSTER EXPECTATIONS")
print("="*70)

for n in range(6, 10):
    print(f"\nn={n}:")
    h1 = compute_joint(n, [0, 1])
    print(f"  h_1(n) = E[Z_0Z_1] = {h1} = {float(h1):.10f}")

    for gap in range(0, n-4):
        l = 2 + gap  # second domino starts at position l
        if l + 1 >= n - 1:
            break
        joint = compute_joint(n, [0, 1, l, l+1])
        product = h1 * compute_joint(n, [l, l+1])  # Should be h1²
        ratio = joint / (h1 * h1) if h1 != 0 else None
        match = "✓ FACTORIZES" if joint == h1 * h1 else "✗ DOES NOT FACTORIZE"
        print(f"  gap={gap}: E[Z_0Z_1·Z_{l}Z_{l+1}] = {joint} = {float(joint):.10f}")
        print(f"           h_1² = {h1*h1} = {float(h1*h1):.10f}")
        print(f"           ratio = {ratio}  {match}")

# Investigate the correction factor
print("\n" + "="*70)
print("CORRECTION FACTORS: E[Z_0Z_1Z_lZ_{l+1}] / h_1²")
print("="*70)

for n in range(6, 10):
    h1 = compute_joint(n, [0, 1])
    h1_sq = h1 * h1
    print(f"\nn={n}: h_1 = {h1}")
    for gap in range(0, n-4):
        l = 2 + gap
        if l + 1 >= n - 1:
            break
        joint = compute_joint(n, [0, 1, l, l+1])
        ratio = joint / h1_sq
        # Express as function of n and gap
        print(f"  gap={gap}: ratio = {ratio} = {float(ratio):.10f}")

# Now test factorization for larger clusters
print("\n" + "="*70)
print("THREE-CLUSTER FACTORIZATION")
print("="*70)

for n in range(8, 10):
    h1 = compute_joint(n, [0, 1])
    print(f"\nn={n}:")
    # Three isolated dominos at 0,1 | 3,4 | 6,7
    joint_3 = compute_joint(n, [0, 1, 3, 4, 6, 7])
    pair_03 = compute_joint(n, [0, 1, 3, 4])
    pair_06 = compute_joint(n, [0, 1, 6, 7])
    pair_36 = compute_joint(n, [3, 4, 6, 7])
    h1_cubed = h1**3

    print(f"  E[Z_0Z_1·Z_3Z_4·Z_6Z_7] = {joint_3} = {float(joint_3):.12f}")
    print(f"  h_1³ = {h1_cubed} = {float(h1_cubed):.12f}")
    print(f"  ratio = {joint_3/h1_cubed}")
    print(f"  E[01·34]·h_1 = {pair_03*h1} = {float(pair_03*h1):.12f}")
    print(f"  h_1·E[34·67] = {h1*pair_36} = {float(h1*pair_36):.12f}")

# Focus on the correction: for two dominos separated by gap g in perm of size n,
# compute the ratio E[Z_0Z_1Z_{2+g}Z_{3+g}] / (E[Z_0Z_1])²
# Conjecture: this ratio depends on n and g, and determines the polynomial structure.
print("\n" + "="*70)
print("RATIO TABLE: E[pair_gap_g] / h_1²")
print("="*70)

print(f"\n{'n':>3} {'gap=0':>15} {'gap=1':>15} {'gap=2':>15} {'gap=3':>15} {'gap=4':>15}")
for n in range(5, 10):
    h1 = compute_joint(n, [0, 1])
    h1_sq = h1 * h1
    row = f"{n:>3}"
    for gap in range(5):
        l = 2 + gap
        if l + 1 >= n - 1:
            row += f" {'---':>15}"
        else:
            joint = compute_joint(n, [0, 1, l, l+1])
            ratio = joint / h1_sq
            row += f" {float(ratio):>15.8f}"
    print(row)

# The gap=0 case: two consecutive dominos = one cluster of size 2
# E[Z_0Z_1Z_2Z_3] = h_2(n) = 2/(n)_4
# h_1² = (2/(n)_2)² = 4/((n)_2)²
# ratio = h_2/h_1² = 2·((n)_2)² / (4·(n)_4) = ((n)_2)² / (2·(n)_4)
# (n)_2 = n(n-1), (n)_4 = n(n-1)(n-2)(n-3)
# ratio = n²(n-1)² / (2n(n-1)(n-2)(n-3)) = n(n-1) / (2(n-2)(n-3))

print("\nGap=0 predicted: n(n-1)/(2(n-2)(n-3)):")
for n in range(5, 10):
    pred = Fraction(n*(n-1), 2*(n-2)*(n-3))
    print(f"  n={n}: {pred} = {float(pred):.8f}")

# For gap=1: E[Z_0Z_1Z_3Z_4] — positions 0,1 and 3,4 with gap 1
# This involves 5 consecutive positions. Let's see the ratio.
# The ratio should be a rational function of n.

print("\nGap≥1 ratios as fractions:")
for n in range(6, 10):
    h1 = compute_joint(n, [0, 1])
    h1_sq = h1 * h1
    for gap in range(1, n-4):
        l = 2 + gap
        if l + 1 >= n - 1:
            break
        joint = compute_joint(n, [0, 1, l, l+1])
        ratio = joint / h1_sq
        print(f"  n={n}, gap={gap}: ratio = {ratio}")

print("\nDone!")
