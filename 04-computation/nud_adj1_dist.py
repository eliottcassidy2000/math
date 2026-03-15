#!/usr/bin/env python3
"""Compute distribution of unit ascents (adj1) in NUD permutations.
N(n,j) = #{σ ∈ NUD(n) : adj1(σ) = j}
W(n) = Σ_j 2^j N(n,j)

Use bitmask DP: dp[mask][last][adj1_count] = number of permutations.
For n up to ~14.
"""
from collections import defaultdict
import sys

def compute_nud_adj1_dist(n):
    """Returns dict: j -> N(n,j) and W(n) check."""
    # dp[mask][last] = dict of adj1_count -> count
    # Too memory-heavy with dict. Use array approach.
    # Since max adj1 ≤ n-1, use array of size n.

    # dp[mask][last][j] = count
    # mask: bitmask of used elements (0..n-1)
    # last: last element placed (0..n-1)
    # j: number of unit ascents so far

    # Memory: 2^n * n * n * 8 bytes. For n=14: 16384 * 14 * 14 * 8 = 25 MB. OK.
    # For n=16: 65536 * 16 * 16 * 8 = 131 MB. Tight but OK.

    max_adj1 = n  # conservative upper bound
    dp = {}  # (mask, last) -> list of length max_adj1

    # Initialize: place each element as first
    for v in range(n):
        mask = 1 << v
        counts = [0] * max_adj1
        counts[0] = 1
        dp[(mask, v)] = counts

    full_mask = (1 << n) - 1

    for step in range(1, n):
        new_dp = {}
        for (mask, last), counts in dp.items():
            if bin(mask).count('1') != step:
                continue
            for v in range(n):
                if mask & (1 << v):
                    continue
                # Check NUD constraint: last - v != 1 (no unit descent)
                if last - v == 1:
                    continue
                new_mask = mask | (1 << v)
                # Check unit ascent: v - last == 1
                is_asc = 1 if (v - last == 1) else 0

                new_counts = [0] * max_adj1
                for j in range(max_adj1 - is_asc):
                    if counts[j] > 0:
                        new_counts[j + is_asc] += counts[j]

                key = (new_mask, v)
                if key in new_dp:
                    for j in range(max_adj1):
                        new_dp[key][j] += new_counts[j]
                else:
                    new_dp[key] = new_counts

        # Merge into dp
        dp.update(new_dp)
        # Remove old entries with wrong popcount
        dp = {k: v for k, v in dp.items() if bin(k[0]).count('1') == step + 1}

    # Sum over all complete permutations
    dist = [0] * max_adj1
    for (mask, last), counts in dp.items():
        if mask == full_mask:
            for j in range(max_adj1):
                dist[j] += counts[j]

    return dist

print("n | NUD(n) | dist(adj1) | W(n) = Σ 2^j N(n,j)")
print("-" * 80)

from fractions import Fraction

for n in range(1, 16):
    dist = compute_nud_adj1_dist(n)
    nud = sum(dist)
    wn = sum(2**j * dist[j] for j in range(len(dist)))

    # Trim trailing zeros
    while dist and dist[-1] == 0:
        dist.pop()

    print(f"n={n:2d}: NUD={nud:>12d}  W={wn:>15d}  dist={dist}")

    # Also print mean and variance of adj1
    if nud > 0:
        mean = Fraction(sum(j * dist[j] for j in range(len(dist))), nud)
        var = Fraction(sum(j**2 * dist[j] for j in range(len(dist))), nud) - mean**2
        print(f"       E[adj1]={float(mean):.6f}  Var[adj1]={float(var):.6f}")

# Print the bivariate table N(n,j) nicely
print("\n=== N(n,j) triangle ===")
all_dist = {}
for n in range(1, 16):
    all_dist[n] = compute_nud_adj1_dist(n)
    while all_dist[n] and all_dist[n][-1] == 0:
        all_dist[n].pop()

for n in range(1, 16):
    d = all_dist[n]
    print(f"n={n:2d}: {d}")

# Check: does N(n,j) satisfy any recurrence?
print("\n=== Checking if rows of N(n,j) relate to previous rows ===")
for n in range(3, 16):
    d = all_dist[n]
    dp = all_dist[n-1]
    dpp = all_dist[n-2] if n >= 3 else [0]
    # Try: N(n,j) = (n-1)*N(n-1,j) + (n-2)*N(n-2,j) + correction_j
    # (NUD recurrence is NUD(n) = (n-1)*NUD(n-1) + (n-2)*NUD(n-2))
    max_j = max(len(d), len(dp) + 1, len(dpp) + 1)
    corr = []
    for j in range(max_j):
        nj = d[j] if j < len(d) else 0
        np1j = dp[j] if j < len(dp) else 0
        np2j = dpp[j] if j < len(dpp) else 0
        c = nj - (n-1)*np1j - (n-2)*np2j
        corr.append(c)
    while corr and corr[-1] == 0:
        corr.pop()
    print(f"n={n:2d}: correction = {corr}")

# Try: N(n,j) = (n-1)*N(n-1,j) + (n-2)*N(n-2,j) + A*N(n-1,j-1) + B*N(n-2,j-1)
# The corrections above should reveal this structure
print("\n=== Checking correction at j-1 level ===")
for n in range(3, 16):
    d = all_dist[n]
    dp = all_dist[n-1]
    dpp = all_dist[n-2] if n >= 3 else [0]
    max_j = max(len(d), len(dp) + 1, len(dpp) + 1)
    for j in range(max_j):
        nj = d[j] if j < len(d) else 0
        np1j = dp[j] if j < len(dp) else 0
        np2j = dpp[j] if j < len(dpp) else 0
        corr = nj - (n-1)*np1j - (n-2)*np2j
        # Check if corr = A*N(n-1,j-1) + B*N(n-2,j-1) for some A, B
        if j >= 1:
            np1jm1 = dp[j-1] if j-1 < len(dp) else 0
            np2jm1 = dpp[j-1] if j-1 < len(dpp) else 0
            if np1jm1 > 0 or np2jm1 > 0:
                pass  # will check below
    break  # Just print for one n to see

# Actually, let me look at the correction vectors more carefully
print("\n=== Correction N(n,j) - (n-1)*N(n-1,j) - (n-2)*N(n-2,j) ===")
print("     j=0, j=1, j=2, j=3, j=4, ...")
for n in range(3, 14):
    d = all_dist[n]
    dp = all_dist[n-1]
    dpp = all_dist[n-2]
    corr = []
    for j in range(max(len(d), len(dp), len(dpp))):
        nj = d[j] if j < len(d) else 0
        np1j = dp[j] if j < len(dp) else 0
        np2j = dpp[j] if j < len(dpp) else 0
        corr.append(nj - (n-1)*np1j - (n-2)*np2j)
    while corr and corr[-1] == 0:
        corr.pop()
    print(f"  n={n:2d}: {corr}")
