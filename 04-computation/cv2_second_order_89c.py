#!/usr/bin/env python3
"""
cv2_second_order_89c.py — Pin down the exact second-order term in CV²
opus-2026-03-14-S89c

We proved CV² = 2/n - 2/(n(n-1)) + R(n) where R(n) = O(1/n²).
Now compute R(n) exactly and find the limit of n²·R(n).

The |S|=4 contribution should give the second-order term.
Specifically, from two disjoint adjacent pairs {j,j+1,k,k+1}:
  E[Z_jZ_{j+1}Z_kZ_{k+1}] = E[Z_jZ_{j+1}]·E[Z_kZ_{k+1}] + correction
"""

from fractions import Fraction
from itertools import permutations, combinations
from math import factorial

def compute_full_expansion(n):
    """Compute E[∏(1+Z_j)] by exact subset expansion."""
    N = factorial(n)
    m = n - 1  # number of Z variables

    # Precompute Z values for all permutations
    all_Z = []
    for perm in permutations(range(n)):
        Z = []
        for j in range(m):
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            Z.append(xj - yj)
        all_Z.append(Z)

    # Compute contribution by subset size
    contributions = {}
    for size in range(0, m+1):
        total = Fraction(0)
        for S in combinations(range(m), size):
            # E[∏_{j∈S} Z_j]
            moment_sum = 0
            for Z in all_Z:
                prod = 1
                for j in S:
                    prod *= Z[j]
                moment_sum += prod
            total += Fraction(moment_sum, N)
        contributions[size] = total

    return contributions

print("="*70)
print("EXACT SUBSET-SIZE CONTRIBUTIONS TO E[∏(1+Z_j)]")
print("="*70)

for n in range(3, 10):
    print(f"\nn={n}:")
    contributions = compute_full_expansion(n)

    running_sum = Fraction(0)
    for size in sorted(contributions):
        c = contributions[size]
        running_sum += c
        if c != 0:
            print(f"  |S|={size}: {c} = {float(c):.10f}")

    # The result should be W(n)/n!
    print(f"  Total: {running_sum} = {float(running_sum):.10f}")

    # Extract residual after |S|=0 and |S|=2
    order_2_term = contributions.get(2, Fraction(0))
    residual = running_sum - 1 - order_2_term
    print(f"  |S|=0 contribution: 1")
    print(f"  |S|=2 contribution: {order_2_term} = {float(order_2_term):.10f}")
    print(f"  Residual (|S|≥4): {residual} = {float(residual):.10f}")
    print(f"  n²×residual: {float(n*n*residual):.10f}")

    cv2 = running_sum - 1
    leading = Fraction(2, n) - Fraction(2, n*(n-1))
    R = cv2 - leading
    print(f"  R(n) = CV² - 2/n + 2/(n(n-1)) = {R} = {float(R):.10f}")
    print(f"  n²·R(n) = {float(n*n*R):.10f}")
    print(f"  n³·R(n) = {float(n*n*n*R):.10f}")

# Now analyze the |S|=4 contribution more carefully
print("\n" + "="*70)
print("|S|=4 STRUCTURAL ANALYSIS")
print("="*70)

for n in range(5, 9):
    N = factorial(n)
    m = n - 1

    # Precompute Z values
    all_Z = []
    for perm in permutations(range(n)):
        Z = []
        for j in range(m):
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            Z.append(xj - yj)
        all_Z.append(Z)

    print(f"\nn={n}:")

    # Classify |S|=4 subsets by structure
    two_adj_pairs = Fraction(0)  # {j,j+1,k,k+1} with k>=j+2
    four_consec = Fraction(0)    # {j,j+1,j+2,j+3}
    three_plus_one = Fraction(0) # {j,j+1,j+2,k} etc
    other = Fraction(0)

    for S in combinations(range(m), 4):
        moment_sum = sum(Z[S[0]]*Z[S[1]]*Z[S[2]]*Z[S[3]] for Z in all_Z)
        moment = Fraction(moment_sum, N)

        # Classify
        s = sorted(S)
        gaps = [s[i+1]-s[i] for i in range(3)]

        if gaps == [1, 1, 1]:
            four_consec += moment
        elif gaps[0] == 1 and gaps[2] == 1 and gaps[1] >= 2:
            # Two adjacent pairs
            two_adj_pairs += moment
        elif sum(1 for g in gaps if g == 1) == 2:
            # Three consecutive + 1 separate
            three_plus_one += moment
        elif sum(1 for g in gaps if g == 1) == 1:
            # One adjacent pair + 2 singletons
            other += moment
        else:
            other += moment

    total_4 = four_consec + two_adj_pairs + three_plus_one + other
    print(f"  Four consecutive: {four_consec} = {float(four_consec):.10f}")
    print(f"  Two adj pairs:    {two_adj_pairs} = {float(two_adj_pairs):.10f}")
    print(f"  Three+one:        {three_plus_one} = {float(three_plus_one):.10f}")
    print(f"  Other:            {other} = {float(other):.10f}")
    print(f"  Total |S|=4:      {total_4} = {float(total_4):.10f}")

    # Expected from independent adjacent pairs:
    # Each pair contributes 2/(n(n-1)), so two pairs: 4/(n(n-1))²
    # Number of disjoint pairs {j,j+1},{k,k+1} with k>=j+2: C(n-2,2) - (n-3) = ... hmm
    # Actually: choose two from n-2 adjacent pairs, minus overlapping = C(n-2,2) - (n-3)
    n_adj_pairs = n - 2
    n_disjoint = n_adj_pairs * (n_adj_pairs - 1) // 2 - (n_adj_pairs - 1)
    # Wait: the adjacent pairs are {0,1}, {1,2}, ..., {n-3,n-2}.
    # Two pairs are "overlapping" if they share an element, i.e., {j,j+1} and {j+1,j+2}.
    # Number of such overlapping pairs: n-3.
    # Total C(n-2, 2) pairs of adjacent pairs; non-overlapping: C(n-2,2) - (n-3).
    n_pairs_of_adj = (n-2)*(n-3)//2
    n_overlapping = n - 3
    n_nonoverlapping = n_pairs_of_adj - n_overlapping

    independent_approx = n_nonoverlapping * Fraction(4, n*(n-1)*n*(n-1))
    print(f"  # non-overlapping adj pair-pairs: {n_nonoverlapping}")
    print(f"  Independent approx for two-pair: {float(independent_approx):.10f}")
    print(f"  Actual two-pair: {float(two_adj_pairs):.10f}")

print("\n" + "="*70)
print("SECOND ORDER COEFFICIENT ANALYSIS")
print("="*70)

# Compute R(n) = CV² - 2/n + 2/(n(n-1)) for larger n using W values from C code
# W values from nud_weight.c:
W_known = {
    3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752,
    9: 439670, 10: 4327904,
}

for n in sorted(W_known):
    W = W_known[n]
    nf = factorial(n)
    cv2 = Fraction(W, nf) - 1
    leading = Fraction(2, n) - Fraction(2, n*(n-1))
    R = cv2 - leading
    print(f"n={n}: R = {R} = {float(R):.12f}")
    print(f"   n²·R = {float(n*n*R):.10f}")
    print(f"   n³·R = {float(n*n*n*R):.10f}")

# The |S|=4 "two disjoint adjacent pairs" should give the dominant second-order term
# E[Z_jZ_{j+1}·Z_kZ_{k+1}] for non-overlapping pairs
# If approximately independent: ≈ (2/(n(n-1)))² = 4/(n²(n-1)²)
# Sum over ≈ C(n-2,2) such pairs minus overlaps ≈ n²/2 pairs
# Total ≈ n²/2 × 4/n⁴ = 2/n²
# But also need 4-consecutive terms: E[Z_jZ_{j+1}Z_{j+2}Z_{j+3}]
# = P(all 4 same sign) × 2 = 2/(n(n-1)(n-2)(n-3))
# Sum over n-4 such: 2(n-4)/(n(n-1)(n-2)(n-3)) ≈ 2/n³

# More precise: the |S|=4 contribution includes:
# 1. C(n-2,2)-(n-3) non-overlapping adjacent pairs × exact expectation
# 2. (n-3) overlapping adjacent pairs (these form 4-consecutive)
#    But 4-consecutive = overlapping adjacent pairs = {j,j+1} and {j+1,j+2} gives {j,j+1,j+2}
#    which is |S|=3, not |S|=4. Wait, overlapping pairs {j,j+1,j+2,j+3}={j,j+1}∪{j+2,j+3}
#    is non-overlapping! {j,j+1} and {j+1,j+2} gives {j,j+1,j+2} which is |S|=3.

print("\nDone!")
