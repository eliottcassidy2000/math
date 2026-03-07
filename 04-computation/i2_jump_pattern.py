#!/usr/bin/env python3
"""
i_2 jump pattern analysis at n=5,6,7 (exhaustive).

Key finding at n=8 (sampling): for each alpha_1 in {4,6,8,10},
the achievable i_2 values SKIP the value needed for H=21.

This script checks if the same pattern holds at smaller n.

Instance: opus-2026-03-07-S40
"""

from itertools import combinations
from collections import Counter, defaultdict
import time


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)
        yield adj


def find_all_odd_cycles(adj, n):
    """Find all directed odd cycles up to length n (if n is odd) or n-1."""
    cycles = []

    # 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i] >> j) & 1 and (adj[j] >> k) & 1 and (adj[k] >> i) & 1:
                    cycles.append((frozenset([i, j, k]), 3))
                elif (adj[i] >> k) & 1 and (adj[k] >> j) & 1 and (adj[j] >> i) & 1:
                    cycles.append((frozenset([i, j, k]), 3))

    # 5-cycles (if n >= 5)
    if n >= 5:
        for verts in combinations(range(n), 5):
            v = list(verts)
            dp = {}
            dp[(1, 0)] = 1
            for S in range(1, 32):
                for i in range(5):
                    if not (S & (1 << i)):
                        continue
                    if (S, i) not in dp:
                        continue
                    c = dp[(S, i)]
                    for j in range(5):
                        if S & (1 << j):
                            continue
                        if (adj[v[i]] >> v[j]) & 1:
                            key = (S | (1 << j), j)
                            dp[key] = dp.get(key, 0) + c
            for j in range(1, 5):
                if (31, j) in dp and (adj[v[j]] >> v[0]) & 1:
                    for _ in range(dp[(31, j)]):
                        cycles.append((frozenset(verts), 5))

    # 7-cycles (if n >= 7)
    if n >= 7:
        for verts in combinations(range(n), 7):
            v = list(verts)
            dp = {}
            dp[(1, 0)] = 1
            for S in range(1, 128):
                for i in range(7):
                    if not (S & (1 << i)):
                        continue
                    if (S, i) not in dp:
                        continue
                    c = dp[(S, i)]
                    for j in range(7):
                        if S & (1 << j):
                            continue
                        if (adj[v[i]] >> v[j]) & 1:
                            key = (S | (1 << j), j)
                            dp[key] = dp.get(key, 0) + c
            for j in range(1, 7):
                if (127, j) in dp and (adj[v[j]] >> v[0]) & 1:
                    for _ in range(dp[(127, j)]):
                        cycles.append((frozenset(verts), 7))

    return cycles


def compute_i2(cycles):
    """Count vertex-disjoint pairs."""
    m = len(cycles)
    i2 = 0
    for a in range(m):
        for b in range(a+1, m):
            if not (cycles[a][0] & cycles[b][0]):
                i2 += 1
    return i2


def held_karp(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
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
                if (adj[v] >> u) & 1:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))


def analyze_n(n):
    print(f"\n{'='*60}")
    print(f"n={n}: exhaustive analysis of (alpha_1, i_2) pairs")
    print(f"{'='*60}")

    alpha_i2_data = defaultdict(list)

    start = time.time()
    for adj in all_tournaments(n):
        cycles = find_all_odd_cycles(adj, n)
        alpha1 = len(cycles)
        i2 = compute_i2(cycles)
        H = held_karp(adj, n)
        alpha_i2_data[alpha1].append((i2, H))

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s")

    for a1 in sorted(alpha_i2_data.keys()):
        data = alpha_i2_data[a1]
        i2_vals = sorted(set(i2 for i2, _ in data))
        h_vals = sorted(set(h for _, h in data))
        i2_dist = Counter(i2 for i2, _ in data)

        # Check if H=21 is achievable
        needed_i2 = None
        if 2*a1 <= 20:  # alpha_1 + 2*i_2 = 10
            rem = 10 - a1
            if rem >= 0 and rem % 2 == 0:
                needed_i2 = rem // 2

        tag = ""
        if needed_i2 is not None:
            if needed_i2 in i2_vals:
                # Check if H=21 actually occurs
                h_at_target = [h for i2, h in data if i2 == needed_i2]
                if 21 in h_at_target:
                    tag = " *** H=21 ACHIEVED! ***"
                else:
                    tag = f" [i_2={needed_i2} exists but H≠21: H={Counter(h_at_target)}]"
            else:
                tag = f" [BLOCKED: need i_2={needed_i2} but i_2 ∈ {{{','.join(map(str,i2_vals))}}}]"

        if a1 <= 12:
            print(f"  alpha_1={a1}: i_2 ∈ {{{','.join(map(str,i2_vals))}}}, "
                  f"H ∈ {{{','.join(map(str,h_vals))}}}, {len(data)} tours{tag}")


def main():
    for n in [5, 6, 7]:
        analyze_n(n)


if __name__ == "__main__":
    main()
