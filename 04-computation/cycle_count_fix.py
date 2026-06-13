"""
cycle_count_fix.py -- kind-pasteur-2026-03-13-S62

CRITICAL FINDING: tr(A^k)/k does NOT count simple directed k-cycles
for k >= 7 in tournaments. Non-simple closed walks of length 7 exist
in 7-vertex tournaments.

Example: 0->1->2->3->0->4->5->0 is a valid closed walk of length 7
visiting vertex 0 twice. This requires 0 beating {1,5} and losing to {3,5}.

For k=3 and k=5: ALL closed walks in tournaments are simple (proved by
exhaustive case analysis on repeat positions — every repeat forces a
bidirectional arc, which is impossible in a tournament).

For k=7: gaps of 4 between repeat positions allow enough room for the
walk to return without forcing bidirectional arcs.

This means:
- c3 = tr(A^3)/3 is CORRECT (always)
- c5 = tr(A^5)/5 is CORRECT (always)
- c7 = tr(A^7)/7 OVERCOUNTS at n >= 7
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_ham_cycles_exact(A, n):
    """Count directed Hamiltonian cycles using Held-Karp DP.
    Fix starting vertex 0, count paths 0->...->v->0."""
    if n < 3:
        return 0

    dp = {}
    # Start at vertex 0
    dp[(1, 0)] = 1

    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):  # Must include vertex 0
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if v == 0 and mask_size < n:
                    continue  # Don't return to 0 early
                prev_mask = mask ^ (1 << v)
                if v == 0:
                    prev_mask = mask  # For the final step, we just check who can reach 0
                    continue
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v] and dp.get((prev_mask, u), 0):
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total

    # Count cycles: full mask, end at some v with v->0 arc
    full = (1 << n) - 1
    total_cycles = 0
    for v in range(1, n):
        if A[v][0] and dp.get((full, v), 0):
            total_cycles += dp[(full, v)]

    return total_cycles

def count_simple_k_cycles(A_full, n, k, verts):
    """Count directed simple k-cycles on vertex set verts."""
    size = len(verts)
    assert size == k

    # Build sub-adjacency
    sub_A = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(k):
            sub_A[i][j] = A_full[verts[i]][verts[j]]

    return count_ham_cycles_exact(sub_A, k)

# ============================================================
# TEST 1: Verify c3 = tr(A^3)/3 is correct
# ============================================================
print("=" * 70)
print("TEST 1: c3 = tr(A^3)/3 vs exact count (n=5, n=7)")
print("=" * 70)

np.random.seed(42)

for n in [5, 7]:
    total_bits = n*(n-1)//2
    mismatches = 0
    n_samples = 200

    for trial in range(n_samples):
        bits = np.random.randint(0, 1 << total_bits)
        A = bits_to_adj(bits, n)

        # tr(A^3)/3
        c3_trace = int(np.trace(np.linalg.matrix_power(A, 3))) // 3

        # Exact count
        c3_exact = 0
        for combo in combinations(range(n), 3):
            c3_exact += count_simple_k_cycles(A, n, 3, list(combo))

        if c3_trace != c3_exact:
            mismatches += 1

    print(f"  n={n}: c3 mismatches = {mismatches}/{n_samples}")

# ============================================================
# TEST 2: Verify c5 = tr(A^5)/5 is correct
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: c5 vs exact count")
print("=" * 70)

np.random.seed(42)

for n in [5, 7]:
    total_bits = n*(n-1)//2
    mismatches = 0
    n_samples = 100

    for trial in range(n_samples):
        bits = np.random.randint(0, 1 << total_bits)
        A = bits_to_adj(bits, n)

        # tr(A^5)/5 for full tournament
        c5_trace = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

        # Exact count: sum over all 5-vertex subsets
        c5_exact = 0
        for combo in combinations(range(n), 5):
            c5_exact += count_simple_k_cycles(A, n, 5, list(combo))

        if c5_trace != c5_exact:
            mismatches += 1
            if mismatches <= 3:
                print(f"    MISMATCH at n={n}: trace={c5_trace}, exact={c5_exact}")

    print(f"  n={n}: c5 mismatches = {mismatches}/{n_samples}")

# ============================================================
# TEST 3: c7 = tr(A^7)/7 vs exact count at n=7
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: c7 = tr(A^7)/7 vs exact count at n=7")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

overcount_dist = Counter()

for trial in range(200):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    # tr(A^7)/7
    c7_trace = int(np.trace(np.linalg.matrix_power(A, 7))) // 7

    # Exact: Hamiltonian cycles in full tournament
    c7_exact = count_ham_cycles_exact(A, n)

    overcount = c7_trace - c7_exact
    overcount_dist[overcount] += 1

    if trial < 10:
        print(f"  Trial {trial}: tr/7={c7_trace}, exact={c7_exact}, "
              f"overcount={overcount}")

print(f"\n  Overcount distribution: {dict(sorted(overcount_dist.items()))}")
print(f"  ALWAYS overcount = 0: {set(overcount_dist.keys()) == {0}}")

# ============================================================
# TEST 4: Fix alpha_1 and verify (H-1)/2 mod 2
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: (H-1)/2 mod 2 = alpha_1 mod 2 with EXACT cycle counts")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

matches = 0
n_samples = 100

for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Exact alpha_1: count ALL simple directed odd cycles
    alpha_1 = 0

    # 3-cycles (always correct via trace)
    for combo in combinations(range(n), 3):
        sub = np.zeros((3, 3), dtype=int)
        verts = list(combo)
        for a in range(3):
            for b in range(3):
                sub[a][b] = A[verts[a]][verts[b]]
        c = int(np.trace(sub @ sub @ sub)) // 3
        alpha_1 += c

    # 5-cycles (always correct via trace)
    for combo in combinations(range(n), 5):
        sub = np.zeros((5, 5), dtype=int)
        verts = list(combo)
        for a in range(5):
            for b in range(5):
                sub[a][b] = A[verts[a]][verts[b]]
        c = int(np.trace(np.linalg.matrix_power(sub, 5))) // 5
        alpha_1 += c

    # 7-cycles: use EXACT counting
    c7_exact = count_ham_cycles_exact(A, n)
    alpha_1 += c7_exact

    expected = alpha_1 % 2
    actual = ((H - 1) // 2) % 2

    if expected == actual:
        matches += 1

print(f"  Matches: {matches}/{n_samples}")

print("\n\nDone.")
