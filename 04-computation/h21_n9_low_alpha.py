#!/usr/bin/env python3
"""
At n=9, check if H=21 is achievable for cycle-rich tournaments.

Key findings so far:
- t3=3 with every vertex in 3-cycle => alpha_3 >= 1 => Part C blocks
- t3=5,6,7 can have alpha_1 as low as 6-9

Need to check: what is the FULL H value for these low-alpha_1 tournaments?
H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
from math import comb
import random

def count_3cycles(adj, n):
    scores = [sum(adj[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def find_directed_cycles_dp(adj, n, k):
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        full = (1 << k) - 1
        for S in range(1, full + 1):
            for i in range(k):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(k):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def held_karp(adj, n):
    """Compute H(T) via Held-Karp."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])


def main():
    n = 9
    random.seed(42)

    print(f"=== n={n}: Full H analysis for cycle-rich, low-alpha_1 tournaments ===")

    h_vals = []
    found = 0

    for trial in range(2000000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        scores = [sum(adj[i]) for i in range(n)]
        if 0 in scores or (n-1) in scores:
            continue

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 10 or t3 < 3:
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        found += 1

        # Compute H via Held-Karp
        H = held_karp(adj, n)
        h_vals.append(H)

        if H == 21:
            print(f"  !!! H=21 FOUND! t3={t3}, scores={tuple(sorted(scores))}")
            # Print adjacency matrix
            for i in range(n):
                print(f"    {adj[i]}")
            break

        if found <= 20:
            print(f"  #{found}: H={H}, t3={t3}, scores={tuple(sorted(scores))}")

    print(f"\nTotal trials: {trial+1}, cycle-rich found: {found}")
    if h_vals:
        print(f"H range: [{min(h_vals)}, {max(h_vals)}]")
        print(f"Min H values: {sorted(set(h_vals))[:20]}")
        h_counter = Counter(h_vals)
        print(f"H=21 count: {h_counter.get(21, 0)}")
        print(f"H < 21 count: {sum(1 for h in h_vals if h < 21)}")
        print(f"H=21 in spectrum: {21 in h_counter}")

        # Distribution of small H
        small_h = {h: c for h, c in h_counter.items() if h <= 30}
        print(f"\nH <= 30 distribution: {dict(sorted(small_h.items()))}")


if __name__ == "__main__":
    main()
