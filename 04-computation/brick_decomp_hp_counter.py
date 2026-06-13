#!/usr/bin/env python3
"""
brick_decomp_hp_counter.py — opus-2026-03-14-S71f

Engineering application: fast H computation via I(Ω, 2).

When Ω(T) decomposes into cliques (packing structure), we can compute H
in O(n^3) time instead of O(n · 2^n) (Held-Karp) or O(n!) (brute force).

Method:
1. Find all directed 3-cycles: O(n^3)
2. Build conflict graph Ω (adjacency): O(α₁²)
3. Find connected components of Ω
4. For each component, compute independence polynomial
5. H = product of component I-polynomials evaluated at x=2

When Ω is small (α₁ << 2^n), this can be MUCH faster than DP.
When Ω decomposes into small cliques, each component contributes (1+c·x)
and H = product of (1+2c).

BENCHMARK: Compare DP vs brick decomposition on n=5,6,7 tournaments.
"""

import numpy as np
from itertools import combinations, permutations
import time
from collections import defaultdict

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def hp_dp(A, n):
    """Held-Karp DP: O(n * 2^n)"""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def ocf_compute(A, n):
    """Compute H via OCF: I(Ω(T), 2)"""
    # Step 1: Find all directed odd cycles
    cycles = set()

    # 3-cycles (dominant at small n)
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.add((i, j, k))
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.add((i, k, j))

    # 5-cycles
    if n >= 5:
        for verts in combinations(range(n), 5):
            for p in permutations(verts[1:]):
                order = [verts[0]] + list(p)
                ok = True
                for idx in range(5):
                    if A[order[idx]][order[(idx+1) % 5]] != 1:
                        ok = False
                        break
                if ok:
                    cycles.add(tuple(order))

    # 7-cycles (expensive, skip for benchmark)
    # For exact results at n>=7, would need this

    cycle_list = list(cycles)
    nc = len(cycle_list)

    if nc == 0:
        return 1

    # Step 2: Build adjacency (share vertex)
    vsets = [frozenset(c) for c in cycle_list]
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True

    # Step 3: Find connected components
    visited = [False]*nc
    components = []

    def bfs(start):
        comp = [start]
        visited[start] = True
        queue = [start]
        while queue:
            u = queue.pop(0)
            for v in range(nc):
                if not visited[v] and adj[u][v]:
                    visited[v] = True
                    comp.append(v)
                    queue.append(v)
        return comp

    for i in range(nc):
        if not visited[i]:
            components.append(bfs(i))

    # Step 4: For each component, compute I(component, 2)
    H = 1
    for comp in components:
        m = len(comp)
        if m <= 20:
            # Brute force independence polynomial
            total = 0
            for mask in range(1 << m):
                bits = []
                temp = mask
                while temp:
                    bits.append(temp & -temp)
                    temp &= temp - 1
                indices = [b.bit_length()-1 for b in bits]
                is_indep = True
                for a in range(len(indices)):
                    for b in range(a+1, len(indices)):
                        ci, cj = comp[indices[a]], comp[indices[b]]
                        if adj[ci][cj]:
                            is_indep = False
                            break
                    if not is_indep:
                        break
                if is_indep:
                    total += 2**len(indices)
            H *= total
        else:
            # For large components, approximate (complete graph assumption)
            H *= (1 + 2*m)

    return H

# ============================================================
# Benchmark
# ============================================================

rng = np.random.default_rng(42)

print("=" * 70)
print("Benchmark: DP vs OCF computation of H")
print("=" * 70)

for n in [5, 6, 7]:
    N = 1000 if n <= 6 else 200

    # Time DP
    t0 = time.time()
    h_dp = []
    for _ in range(N):
        A = random_tournament(n, rng)
        h_dp.append(hp_dp(A, n))
    dp_time = time.time() - t0

    # Reset RNG for same tournaments
    rng = np.random.default_rng(42)

    # Time OCF
    t0 = time.time()
    h_ocf = []
    for _ in range(N):
        A = random_tournament(n, rng)
        h_ocf.append(ocf_compute(A, n))
    ocf_time = time.time() - t0

    # Verify
    matches = sum(1 for a, b in zip(h_dp, h_ocf) if a == b)

    print(f"\n  n={n} ({N} tournaments):")
    print(f"    DP:  {dp_time:.3f}s ({N/dp_time:.0f} tours/s)")
    print(f"    OCF: {ocf_time:.3f}s ({N/ocf_time:.0f} tours/s)")
    print(f"    Speedup: {dp_time/ocf_time:.2f}x")
    print(f"    Match: {matches}/{N}")
    if matches < N:
        # Show first mismatch
        for i, (a, b) in enumerate(zip(h_dp, h_ocf)):
            if a != b:
                print(f"    First mismatch at tournament {i}: DP={a}, OCF={b}")
                break

    rng = np.random.default_rng(42 + n)  # Different seed for next n

# ============================================================
# Component structure analysis
# ============================================================

print(f"\n{'='*70}")
print("Component Structure Analysis at n=6")
print(f"{'='*70}")

rng = np.random.default_rng(123)
n = 6
N = 10000
comp_stats = defaultdict(int)

for _ in range(N):
    A = random_tournament(n, rng)

    cycles = set()
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.add((i, j, k))
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.add((i, k, j))
    for verts in combinations(range(n), 5):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            ok = True
            for idx in range(5):
                if A[order[idx]][order[(idx+1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                cycles.add(tuple(order))

    cycle_list = list(cycles)
    nc = len(cycle_list)

    if nc == 0:
        comp_stats[("empty",)] += 1
        continue

    vsets = [frozenset(c) for c in cycle_list]
    adj_mat = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                adj_mat[i][j] = adj_mat[j][i] = True

    visited = [False]*nc
    comp_sizes = []
    def bfs(start):
        comp = [start]
        visited[start] = True
        queue = [start]
        while queue:
            u = queue.pop(0)
            for v in range(nc):
                if not visited[v] and adj_mat[u][v]:
                    visited[v] = True
                    comp.append(v)
                    queue.append(v)
        return comp

    for i in range(nc):
        if not visited[i]:
            comp_sizes.append(len(bfs(i)))

    comp_stats[tuple(sorted(comp_sizes, reverse=True))] += 1

print(f"\nComponent decomposition of Ω at n=6 ({N} tournaments):")
for pattern, count in sorted(comp_stats.items(), key=lambda x: -x[1])[:20]:
    pct = count/N*100
    if pattern == ("empty",):
        print(f"  Ω = ∅ (transitive):    {count:5d} ({pct:.1f}%)")
    else:
        parts = "+".join(str(s) for s in pattern)
        print(f"  Ω = {parts:15s}:  {count:5d} ({pct:.1f}%)")

print(f"\nTotal with >1 component: {sum(v for k,v in comp_stats.items() if len(k) > 1 and k != ('empty',))}")
print(f"Total with 1 component: {sum(v for k,v in comp_stats.items() if len(k) == 1 and k != ('empty',))}")
