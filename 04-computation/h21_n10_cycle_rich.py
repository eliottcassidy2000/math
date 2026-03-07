#!/usr/bin/env python3
"""
At n=10, check cycle-rich tournaments: min H and min decomposition.

If min decomp > 10 at n=10: cycle-rich at n=10 also blocks H=21.
Pattern: n=8 min H=25, n=9 min H=45. Should grow further.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
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

def held_karp(adj, n):
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
    n = 10
    random.seed(42)

    print(f"=== n={n}: cycle-rich min-H check ===")

    h_vals = []
    found = 0
    min_H = float('inf')

    for trial in range(1000000):
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
        if t3 > 15:  # Relaxed threshold
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        found += 1

        H = held_karp(adj, n)
        h_vals.append(H)

        if H < min_H:
            min_H = H
            print(f"  New min H={H}, t3={t3}, scores={tuple(sorted(scores))}")

        if H == 21:
            print(f"  !!! H=21 FOUND !!!")
            break

        if found % 50 == 0:
            print(f"  ... {found} found, min H so far = {min_H}")

    print(f"\nTotal trials: {trial+1}, cycle-rich found: {found}")
    if h_vals:
        print(f"H range: [{min(h_vals)}, {max(h_vals)}]")
        print(f"H=21 found: {21 in h_vals}")
        sorted_h = sorted(set(h_vals))
        print(f"Smallest H values: {sorted_h[:15]}")

if __name__ == "__main__":
    main()
