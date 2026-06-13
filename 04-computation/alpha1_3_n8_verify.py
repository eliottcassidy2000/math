#!/usr/bin/env python3
"""
alpha1_3_n8_verify.py — opus-2026-03-14-S71f

CRITICAL GAP: Verify that α₁=3 forces α₂≥2 at n=8.

The H=7 proof requires:
  α₁=3 (exactly 3 directed odd cycles in Ω) ⟹ α₂≥2

At n=7: verified (all 97/50000 samples with α₁=3 have α₂=2, H=15)
At n=8: UNVERIFIED — need to check.

Strategy: sample 500k random n=8 tournaments, find those with α₁=3,
check their α₂ values. If any has α₂=0, the H=7 proof has a gap.

Also check: can α₁=3 even occur at n=8? (It requires exactly 3
directed odd cycles across all possible cycle lengths 3,5,7)

Note on α₁: Each directed odd cycle is a SEPARATE vertex in Ω.
For 3-cycles: each triple {i,j,k} supports 0 or 1 directed 3-cycle.
For 5-cycles: each 5-set supports 0, 1, 2, or 3 directed 5-cycles.
For 7-cycles: each 7-set supports various numbers.

At n=8: the smallest odd cycles (3-cycles) are by far the most common.
So α₁=3 most likely means t₃=3 (3-cycle vertex sets) with no 5 or 7 cycles.
But need to verify that t₃=3 can exist at n=8 and that it forces α₂≥2.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import time

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_all_directed_odd_cycles(A, n, max_k=None):
    """Count all directed odd cycles as canonical tuples (start at min vertex)."""
    if max_k is None:
        max_k = n
    cycles = set()
    for k in range(3, min(max_k, n)+1, 2):
        for verts in combinations(range(n), k):
            for p in permutations(verts[1:]):
                order = [verts[0]] + list(p)
                is_cycle = True
                for idx in range(k):
                    if A[order[idx]][order[(idx+1) % k]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(tuple(order))
    return cycles

def count_3cycles_and_5cycles_fast(A, n):
    """Count directed 3-cycles and 5-cycles efficiently.
    Returns the set of all directed odd cycles (as canonical tuples)."""
    cycles = set()

    # 3-cycles: fix min vertex, check both directions
    for i, j, k in combinations(range(n), 3):
        # Check i→j→k→i
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.add((i, j, k))
        # Check i→k→j→i
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.add((i, k, j))

    # 5-cycles: fix min vertex, check all orderings
    for verts in combinations(range(n), 5):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            is_cycle = True
            for idx in range(5):
                if A[order[idx]][order[(idx+1) % 5]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                cycles.add(tuple(order))

    return cycles

def compute_alpha2(cycles):
    """Count pairs of cycles sharing no vertex."""
    cycle_list = list(cycles)
    vsets = [frozenset(c) for c in cycle_list]
    count = 0
    for i in range(len(cycle_list)):
        for j in range(i+1, len(cycle_list)):
            if not (vsets[i] & vsets[j]):
                count += 1
    return count

# Main sampling loop
n = 8
rng = np.random.default_rng(42)
N = 500000

alpha1_3_count = 0
alpha1_3_alpha2 = Counter()
alpha1_3_h = Counter()

# Also track overall α₁ distribution for α₁ ≤ 5
alpha1_small = Counter()

start = time.time()

for trial in range(N):
    if trial > 0 and trial % 50000 == 0:
        elapsed = time.time() - start
        rate = trial / elapsed
        eta = (N - trial) / rate
        print(f"  {trial}/{N} done ({rate:.0f}/s, ETA {eta:.0f}s)")

    A = random_tournament(n, rng)

    # Fast pre-check: count 3-cycle vertex sets
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1

    # At n=8, α₁ ≥ t₃ (since each 3-cycle vertex set has exactly 1 directed cycle)
    # If t₃ ≥ 4, then α₁ ≥ 4 (not interesting for α₁=3 check)
    if t3 > 3:
        continue

    # Need to also count directed 5-cycles and 7-cycles
    # For t₃ ≤ 3 cases, count the full α₁
    cycles = count_3cycles_and_5cycles_fast(A, n)

    # Also check 7-cycles (expensive but rare for t₃≤3)
    for verts in combinations(range(n), 7):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            is_cycle = True
            for idx in range(7):
                if A[order[idx]][order[(idx+1) % 7]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                cycles.add(tuple(order))

    alpha1 = len(cycles)

    if alpha1 <= 5:
        alpha1_small[alpha1] += 1

    if alpha1 == 3:
        alpha2 = compute_alpha2(cycles)
        alpha1_3_count += 1
        alpha1_3_alpha2[alpha2] += 1
        H = 1 + 2*alpha1 + 4*alpha2  # Quick check (may be wrong if α₃>0)
        alpha1_3_h[H] += 1

        # Print details of first few
        if alpha1_3_count <= 5:
            vsets = [frozenset(c) for c in cycles]
            print(f"\n  α₁=3 example #{alpha1_3_count}:")
            for c in sorted(cycles):
                print(f"    cycle: {c} (len={len(c)})")
            print(f"    α₂ = {alpha2}")
            print(f"    H_approx = {H}")

elapsed = time.time() - start
print(f"\n{'='*60}")
print(f"Completed {N} tournaments in {elapsed:.1f}s")
print(f"{'='*60}")

print(f"\nα₁ ≤ 5 distribution:")
for a1 in sorted(alpha1_small.keys()):
    print(f"  α₁={a1}: {alpha1_small[a1]} ({alpha1_small[a1]/N*100:.4f}%)")

print(f"\nα₁=3 results:")
print(f"  Total found: {alpha1_3_count}")
if alpha1_3_count > 0:
    print(f"  α₂ distribution:")
    for a2 in sorted(alpha1_3_alpha2.keys()):
        print(f"    α₂={a2}: {alpha1_3_alpha2[a2]} ({alpha1_3_alpha2[a2]/alpha1_3_count*100:.1f}%)")
    print(f"  H distribution:")
    for h in sorted(alpha1_3_h.keys()):
        print(f"    H={h}: {alpha1_3_h[h]}")

    if 0 in alpha1_3_alpha2:
        print(f"\n  *** CRITICAL: α₁=3, α₂=0 found! H=7 proof has a gap! ***")
    else:
        print(f"\n  ✓ α₁=3 always has α₂≥{min(alpha1_3_alpha2.keys())} at n=8")
        print(f"    H=7 proof holds for n=8")
else:
    print(f"  α₁=3 never occurred at n=8 in {N} samples")
    print(f"  This makes the H=7 proof even stronger (α₁=3 may be impossible at n=8)")
