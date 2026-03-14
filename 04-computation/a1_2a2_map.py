"""
a1_2a2_map.py — Map achievable values of a1 + 2*a2 at n=7
kind-pasteur-2026-03-14-S64

Since H = 1 + 2*(a1 + 2*a2), and H is always odd,
the achievable H values are exactly {1 + 2*k : k in achievable a1+2*a2}.

H=21 impossible iff a1+2*a2=10 impossible.

Strategy: Large random sample with proper RNG, then
check which values of a1+2*a2 = (H-1)/2 appear.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
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
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    n = 7
    N = 50000

    print(f"Sampling {N} random tournaments at n={n}...")
    print(f"Computing H via Hamiltonian path count (no cycle decomposition needed)")

    rng = np.random.default_rng(2026_0313)

    k_values = defaultdict(int)  # k = (H-1)/2 -> count
    h_values = defaultdict(int)

    for trial in range(N):
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        H = count_ham_paths(A)
        assert H % 2 == 1, f"H={H} is even!"
        k = (H - 1) // 2
        k_values[k] += 1
        h_values[H] += 1

        if (trial + 1) % 10000 == 0:
            print(f"  {trial+1}/{N}")

    print(f"\nAchievable k = (H-1)/2 values:")
    max_k = max(k_values.keys())
    min_k = min(k_values.keys())

    gaps = []
    for k in range(min_k, max_k + 1):
        count = k_values.get(k, 0)
        h = 2 * k + 1
        marker = ""
        if count == 0:
            gaps.append(k)
            if k <= 50:  # Only mark small gaps
                marker = f"  <-- GAP (H={h})"
        if k <= 60 or count == 0:
            print(f"  k={k:3d} (H={h:3d}): {count:5d}{marker}")

    print(f"\nGap values of k (a1+2*a2) up to 60: {[g for g in gaps if g <= 60]}")
    print(f"Corresponding forbidden H: {[2*g+1 for g in gaps if g <= 60]}")

    # Density analysis
    print(f"\nDensity by range:")
    for lo, hi in [(0, 10), (10, 20), (20, 30), (30, 50), (50, 70), (70, 100)]:
        total_in_range = sum(k_values.get(k, 0) for k in range(lo, hi+1))
        achievable = sum(1 for k in range(lo, hi+1) if k_values.get(k, 0) > 0)
        possible = hi - lo + 1
        print(f"  k in [{lo:3d}, {hi:3d}]: {achievable}/{possible} values achieved, {total_in_range} tournaments")

    print("\nDone.")

if __name__ == "__main__":
    main()
