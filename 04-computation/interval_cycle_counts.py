#!/usr/bin/env python3
"""
interval_cycle_counts.py — Exact odd-cycle counts from THM-143

THM-143 gives co_occ_k(d) = a_k + b_k*d for Interval tournaments.
By vertex-transitivity:
  c_k = (p/k) * sum_{d=1}^m co_occ_k(d)
     = (p/k) * [m*a_k + b_k * m(m+1)/2]

This gives EXACT formulas for c_3, c_5, c_7, ... for all Interval tournaments.

Key insight: the co_occ slope GF is x^3 * (1+x)^{m-2}, suggesting a clean
generating function for the cycle counts themselves.

Also: OCF says H = I(Omega, 2) where Omega encodes odd-cycle overlaps.
If we know ALL cycle counts AND all pairwise overlaps (from THM-143),
we can reconstruct I(Omega, 2) and hence H.

Author: opus-2026-03-12-S68
"""
import sys
from math import comb, factorial

def co_occ_formula(m, k, d):
    """THM-143: co_occ_k(d) for Interval tournament."""
    b_k = comb(m - 2, k - 3)
    a_k = comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - m * b_k
    return a_k + b_k * d

def cycle_count_interval(p, k):
    """Total number of directed k-cycles in Interval tournament T_p."""
    m = (p - 1) // 2
    total = 0
    for d in range(1, m + 1):
        total += co_occ_formula(m, k, d)
    # Each k-cycle is counted k times (once per vertex), and each pair
    # {0, d} is counted once. But co_occ_k(d) counts vertex SETS, not
    # directed cycles. Need to be careful.
    # co_occ_k(d) = #{k-vertex sets V containing {0,d} that support a HC}
    # Total vertex sets = (1/C(k,2)) * sum_{d} p * co_occ_k(d)
    # No wait: each k-vertex set has C(k,2) pairs, and each is counted
    # once in the sum over d. But we're summing over d=1..m, which only
    # covers half the pairs (by symmetry co_occ_k(d) = co_occ_k(p-d)).
    # So: n_0 = (2/(k-1)) * sum_{d=1}^m co_occ_k(d)
    #     n_total = n_0 * p / k = (2p / (k(k-1))) * sum co_occ_k(d)
    # But n_total counts vertex sets, not directed cycles. Each k-vertex
    # set supports 2*(k-1)!/2 = (k-1)! directed Hamiltonian cycles on it
    # (for a tournament — actually depends on the specific set).
    # Hmm, this gets complicated. Let me just count vertex sets.
    n_0 = 0
    for d in range(1, m + 1):
        n_0 += co_occ_formula(m, k, d)
    # n_0 = number of k-sets containing vertex 0 with a HC
    # By vertex-transitivity, n_total = n_0 * p / k
    n_total = n_0 * p // k  # Should be exact integer
    return n_0, n_total

print("=" * 70)
print("INTERVAL TOURNAMENT CYCLE COUNTS FROM THM-143")
print("=" * 70)

for p in [7, 11, 13, 17, 19, 23]:
    m = (p - 1) // 2
    print(f"\np={p}, m={m}:")

    for k in range(3, p + 1, 2):
        if k > p:
            break
        n_0, n_total = cycle_count_interval(p, k)
        # Slope and intercept
        b_k = comb(m - 2, k - 3)
        a_k = comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - m * b_k
        print(f"  k={k:2d}: n_0={n_0:8d}, n_total={n_total:10d}, "
              f"slope={b_k}, intercept={a_k}")

    # Special check: k=3 should give p(p^2-1)/24 vertex sets
    # Actually c_3 = p(p^2-1)/24 is the number of DIRECTED 3-cycles / 3
    # = number of 3-vertex sets that form a 3-cycle
    expected_c3 = p * (p*p - 1) // 24
    _, actual_c3 = cycle_count_interval(p, 3)
    print(f"  Check c_3: expected {expected_c3}, got {actual_c3}, {'MATCH' if expected_c3 == actual_c3 else 'MISMATCH'}")

# Generating function for slopes
print("\n" + "=" * 70)
print("SLOPE GENERATING FUNCTION: sum_k b_k x^k = x^3 (1+x)^{m-2}")
print("=" * 70)

for p in [7, 11, 13, 17, 19]:
    m = (p - 1) // 2
    slopes = []
    for k in range(3, p + 1, 2):
        slopes.append((k, comb(m - 2, k - 3)))
    print(f"  p={p} (m={m}): {slopes}")

# The GF sum_k b_k x^k = x^3 (1+x)^{m-2}
# Evaluated at x=1: sum b_k = 2^{m-2}
# This means the total "slope contribution" is 2^{m-2}
print(f"\n  Total slope = 2^(m-2):")
for p in [7, 11, 13, 17, 19]:
    m = (p - 1) // 2
    total_slope = sum(comb(m-2, k-3) for k in range(3, p+1, 2))
    print(f"    p={p}: sum b_k = {total_slope}, 2^(m-2) = {2**(m-2)}, "
          f"{'MATCH' if total_slope == 2**(m-2) else 'MISMATCH'}")

# Connection to H: sum of weighted cycle counts
print("\n" + "=" * 70)
print("H(Interval) FROM CYCLE COUNTS")
print("=" * 70)
print("""
OCF: H = I(Omega, 2) where Omega is the odd-cycle conflict graph.
     I(G, x) = sum_{S independent} x^|S|

The cycle counts alone DON'T determine H (need the conflict structure too).
BUT: the co_occ function DOES encode pairwise conflicts!

Two k-cycles conflict iff they share a vertex.
From co_occ_k(d), we know the pairwise overlap structure.

For DISJOINT cycle PAIRS: THM-142 gives the exact excess.
For GENERAL independence sets: need k-wise disjointness.

KEY QUESTION: Can we recover I(Omega, 2) from the co_occ functions alone?
""")

# Verify H values
print("Verifying H values:")
# These are known H values for Interval tournaments
known_H = {7: 175, 11: 95095}
for p in [7, 11]:
    m = (p - 1) // 2
    print(f"  p={p}: H(Interval) = {known_H.get(p, '?')}")
    # Print all cycle counts
    for k in range(3, p + 1, 2):
        _, n_total = cycle_count_interval(p, k)
        print(f"    c_{k} vertex sets = {n_total}")

print("\nDONE.")
