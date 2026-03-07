#!/usr/bin/env python3
"""
At n=9, mixed-direction 3-cycle pairs EXIST.
But does mixed_sum = T_only_sum = T^op_only_sum = H(T)?

The key test: do mixed-direction permutations contribute to the
Rédei-Berge sum, and if so, do they cancel?

Using Held-Karp for H(T) and direct enumeration of cycle permutations.

Note: enumerating all 9! permutations is feasible (362880 perms).

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict
import random

n = 9

def held_karp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))

def analyze_rédei_berge(n, A):
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    H = held_karp(A, n)

    # Classify permutations
    T_only_sum = 0
    Top_only_sum = 0
    mixed_sum = 0
    mixed_contrib = 0  # Contribution from MIXED-only (not T-only or Top-only)

    for sigma_tuple in permutations(range(n)):
        sigma = list(sigma_tuple)
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        all_odd = True
        has_T = False
        has_Top = False
        ok = True
        nontrivial = 0

        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                all_odd = False
                ok = False
                break
            nontrivial += 1

            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set
                       for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_set
                        for i in range(len(cyc)))

            if not is_T and not is_Top:
                ok = False
                break

            if is_T:
                has_T = True
            if is_Top:
                has_Top = True

        if not ok or not all_odd:
            continue

        weight = 2**nontrivial

        if has_T and not has_Top:
            T_only_sum += weight
        elif has_Top and not has_T:
            Top_only_sum += weight
        elif has_T and has_Top:
            mixed_contrib += weight
        else:
            # Identity permutation (no nontrivial cycles)
            T_only_sum += weight
            Top_only_sum += weight
            mixed_contrib += weight  # Don't double count; will handle below

        mixed_sum += weight

    # Adjust: identity counted in all three
    # Actually let me redo the counting more carefully
    # The "mixed_sum" counts ALL valid permutations (T, T^op, or mixed cycles)
    # with weight 2^{nontrivial}. This should equal H(T) per GS.

    return H, T_only_sum, Top_only_sum, mixed_contrib, mixed_sum

# Test with a random n=9 tournament
rng = random.Random(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if rng.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

print(f"Random n=9 tournament (seed=42):")
H, T_only, Top_only, mixed_c, mixed_total = analyze_rédei_berge(n, A)
print(f"  H(T) = {H}")
print(f"  T-only sum:    {T_only}")
print(f"  T^op-only sum: {Top_only}")
print(f"  Mixed contrib: {mixed_c}")
print(f"  Total (mixed): {mixed_total}")
print(f"  T-only + mixed = {T_only + mixed_c}")
print(f"  T-only = H? {T_only == H}")
print(f"  Total = H? {mixed_total == H}")

# Try a second tournament
rng2 = random.Random(7)
A2 = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if rng2.random() < 0.5:
            A2[i][j] = 1
        else:
            A2[j][i] = 1

print(f"\nRandom n=9 tournament (seed=7):")
H2, T2, Top2, mix2, total2 = analyze_rédei_berge(n, A2)
print(f"  H(T) = {H2}")
print(f"  T-only sum:    {T2}")
print(f"  T^op-only sum: {Top2}")
print(f"  Mixed contrib: {mix2}")
print(f"  Total (mixed): {total2}")
print(f"  T-only = H? {T2 == H2}")
print(f"  Total = H? {total2 == H2}")
