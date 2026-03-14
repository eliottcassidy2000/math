#!/usr/bin/env python3
"""
h7_check_v2.py — opus-2026-03-14-S71e

EXHAUSTIVE check: Is H=7 achievable at n=7?

Strategy: enumerate ALL 2^21 tournaments at n=7 and find those with H=7.
This is ~2M tournaments. HP counting takes O(2^n * n) = O(896) per tournament.
Total: ~1.8 billion operations. Might take a few minutes.

Faster approach: since H is always odd, H=7 ↔ (H-1)/2=3 ↔ α₁+2α₂+4α₃=3.
We can first check dc3 quickly to filter.
"""

import sys
import time
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

n = 7
num_edges = n*(n-1)//2  # = 21

def count_hp_fast(adj_bits, n):
    """Count Hamiltonian paths using bitmask DP, with bit-packed adjacency."""
    dp = [0] * (n * (1 << n))
    for v in range(n):
        dp[v * (1 << n) + (1 << v)] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[v * (1 << n) + mask]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj_bits[v] & (1 << u):
                    dp[u * (1 << n) + (mask | (1 << u))] += val

    full = (1 << n) - 1
    return sum(dp[v * (1 << n) + full] for v in range(n))

print("=" * 70)
print("EXHAUSTIVE SEARCH FOR H=7 AT n=7")
print("=" * 70)

# Precompute edge pairs
edges = []
for i in range(n):
    for j in range(i+1, n):
        edges.append((i, j))

# We'll sample strategically: only check tournaments where dc3 is small
# dc3 is determined by score sequence: dc3 = C(n,3) - sum C(s_i, 2)
# C(7,3) = 35
# For dc3=3: sum C(s_i,2) = 32

# H=7 needs α₁+2α₂+4α₃ = 3
# Since α₁ ≥ 0, α₂ ≥ 0, α₃ ≥ 0:
# (α₁,α₂,α₃) ∈ {(3,0,0), (1,1,0)}
# (1,1,0) impossible since α₂ ≤ C(α₁,2) = 0
# So need: α₁=3, α₂=0, α₃=0

# α₁ = dc3 + dc5 + dc7 = 3, no disjoint pairs
# The simplest: dc3=3, dc5=0, dc7=0

# Alternative: dc3=2, dc5=1, dc7=0 → α₁=3
# Or: dc3=1, dc5=1, dc7=1 → α₁=3
# Or: dc3=0, dc5=0, dc7=3 → α₁=3 (impossible, only 1 7-cycle set)
# etc.

# Let's do a FULL exhaustive search using optimized code
# 2^21 = 2097152. With fast HP counting, should take ~5 min.

print(f"\n  Checking all {2**num_edges} tournaments at n={n}...")
print(f"  Looking for H=7 (needs α₁=3, α₂=0, α₃=0)")
print()

t0 = time.time()
h7_count = 0
h_low = Counter()
total_checked = 0

# Process in chunks
for bits in range(2**num_edges):
    # Build adjacency as bitmasks
    adj_bits = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj_bits[i] |= (1 << j)
            else:
                adj_bits[j] |= (1 << i)
            idx += 1

    # Quick filter: count dc3
    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations of 3-cycle
                if (adj_bits[i] & (1<<j)) and (adj_bits[j] & (1<<k)) and (adj_bits[k] & (1<<i)):
                    dc3 += 1
                elif (adj_bits[i] & (1<<k)) and (adj_bits[k] & (1<<j)) and (adj_bits[j] & (1<<i)):
                    dc3 += 1

    # For H=7, need α₁ ≤ 3, so dc3 ≤ 3
    if dc3 > 5:  # Be generous to catch dc3=3 with some dc5
        total_checked += 1
        if total_checked % 500000 == 0:
            elapsed = time.time() - t0
            print(f"  Checked {total_checked}/{2**num_edges} ({100*total_checked/2**num_edges:.1f}%) in {elapsed:.1f}s")
        continue

    H = count_hp_fast(adj_bits, n)

    if H <= 15:
        h_low[H] += 1

    if H == 7:
        h7_count += 1
        # Get details
        print(f"  *** H=7 FOUND at bits={bits}! ***")
        scores = sorted([bin(adj_bits[i]).count('1') for i in range(n)])
        print(f"  Scores: {scores}, dc3={dc3}")

    total_checked += 1
    if total_checked % 500000 == 0:
        elapsed = time.time() - t0
        print(f"  Checked {total_checked}/{2**num_edges} ({100*total_checked/2**num_edges:.1f}%) in {elapsed:.1f}s")
        print(f"  H≤15 so far: {sorted(h_low.items())}")

elapsed = time.time() - t0
print(f"\n  COMPLETE: checked all {total_checked} tournaments in {elapsed:.1f}s")
print(f"\n  H=7 found: {h7_count} tournaments")
print(f"\n  Low H distribution:")
for h in sorted(h_low.keys()):
    print(f"    H={h:3d}: {h_low[h]:8d}")

if h7_count == 0:
    print(f"\n  *** CONFIRMED: H=7 is IMPOSSIBLE at n=7 ***")
    print(f"  HYP-992 is correct (at n=7).")
else:
    print(f"\n  *** H=7 IS POSSIBLE at n=7 ***")
    print(f"  HYP-992 is WRONG!")

print("\nDone.")
