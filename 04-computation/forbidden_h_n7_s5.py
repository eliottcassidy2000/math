#!/usr/bin/env python3
"""
Sample n=7 tournaments to see which H values appear.
Goal: check if H=21, 35, 39 (forbidden at n=6) remain forbidden at n=7.

opus-2026-05-28-S5
"""
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import random
from collections import Counter

random.seed(42)


def count_HP(T, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1<<n):
        for v in range(n):
            if dp[S][v] == 0: continue
            if not ((S >> v) & 1): continue
            for u in range(n):
                if (S >> u) & 1: continue
                if T[v][u]:
                    dp[S | (1<<u)][u] += dp[S][v]
    full = (1<<n) - 1
    return sum(dp[full][v] for v in range(n))


def main():
    n = 7
    SAMPLES = 200000
    h_counts = Counter()
    for trial in range(SAMPLES):
        T = [[False]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    T[i][j] = True
                else:
                    T[j][i] = True
        H = count_HP(T, n)
        h_counts[H] += 1
    print(f"n={n}, {SAMPLES} samples:")
    # show all encountered H, sorted
    seen = sorted(h_counts.keys())
    print(f"  Distinct H values encountered: {len(seen)}")
    print(f"  Min H = {min(seen)}, Max H = {max(seen)}")
    # in range [1, 189], check missing odd values
    target_missing = [7, 21, 35, 39]
    print(f"  Counts for target potentially-forbidden H values:")
    for h in target_missing:
        print(f"    H={h}: {h_counts[h]}")
    # check all odd H in range [1, max_seen]
    max_seen = max(seen)
    odd_missing = [h for h in range(1, max_seen+1, 2) if h not in h_counts]
    print(f"  Odd H values in [1,{max_seen}] NOT seen in sample (could be sampling artifacts for rare ones):")
    print(f"  Count of missing: {len(odd_missing)}")
    print(f"  First 30 missing: {odd_missing[:30]}")
    # show H=21 appearance count and other small H
    print(f"  Count for small H values (sample frequencies):")
    for h in sorted([k for k in h_counts.keys() if k <= 50]):
        print(f"    H={h}: {h_counts[h]}")


if __name__ == "__main__":
    main()
