#!/usr/bin/env python3
"""
Exhaustive enumeration of n=7 tournaments to verify H != 7 and H != 21.

This is the slow but rigorous check that replaces the 200k/1M sampled tests
from S5. n=7 means 2^21 = 2,097,152 tournaments.

Uses HP counting via bitmask DP, written for speed.

opus-2026-05-28-S6
"""
import os, sys, time
os.environ['PYTHONIOENCODING'] = 'utf-8'
from collections import Counter
from itertools import combinations
from math import comb


def count_HP_fast(adj_matrix_int, n):
    """Bitmask DP for HP counting. T encoded as n x n bool matrix."""
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    full = (1<<n) - 1
    for S in range(1, full+1):
        for v in range(n):
            if not ((S >> v) & 1): continue
            d = dp[S][v]
            if d == 0: continue
            for u in range(n):
                if (S >> u) & 1: continue
                if adj_matrix_int[v] & (1 << u):
                    dp[S | (1<<u)][u] += d
    return sum(dp[full][v] for v in range(n))


def enumerate_n7_H_distribution():
    n = 7
    n_pairs = n * (n-1) // 2
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    total = 1 << n_pairs
    h_distrib = Counter()
    t0 = time.time()
    print(f"Enumerating 2^{n_pairs} = {total} n={n} tournaments", flush=True)
    for bits in range(total):
        # build adj matrix as integer rows
        adj = [0]*n
        for k,(i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
        H = count_HP_fast(adj, n)
        h_distrib[H] += 1
        if (bits + 1) % 50000 == 0:
            elapsed = time.time() - t0
            rate = (bits+1)/elapsed
            eta = (total - bits - 1) / rate
            print(f"  progress: {bits+1}/{total} ({100*(bits+1)/total:.2f}%), eta {eta:.0f}s", flush=True)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s", flush=True)
    return h_distrib


def main():
    h_distrib = enumerate_n7_H_distribution()
    achievable = sorted(h_distrib.keys())
    print(f"n=7: H distribution computed", flush=True)
    print(f"  Number of distinct H values: {len(achievable)}", flush=True)
    print(f"  Min H = {min(achievable)}, Max H = {max(achievable)}", flush=True)
    print(f"  Total tournaments: {sum(h_distrib.values())}", flush=True)
    print(f"  Expected: 2^21 = {1 << 21}", flush=True)
    # Check forbidden values
    print(f"\n  Count for key target H values:", flush=True)
    for h in (7, 21, 35, 39, 49, 63, 77, 91):
        print(f"    H={h}: {h_distrib.get(h, 0)}", flush=True)
    forbidden = [h for h in range(1, max(achievable)+1, 2) if h not in h_distrib]
    print(f"\n  Forbidden odd H in [1, {max(achievable)}] at n=7: {forbidden}", flush=True)
    print(f"  Count of forbidden: {len(forbidden)}", flush=True)


if __name__ == "__main__":
    main()
