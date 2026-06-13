#!/usr/bin/env python3
"""
MAJOR FINDING: Quiver mutation preserves H(T) when result is a tournament!

Verified exhaustive at n=4,5 and sampled at n=6: whenever quiver mutation
at vertex k produces a valid tournament, H is preserved. 0 exceptions.

This script:
1. Verifies exhaustively at n=4,5
2. Characterizes WHEN quiver mutation produces a tournament
3. Checks if mutation-preserving = isomorphism (same tournament up to relabeling?)
4. Tests at n=7 with sampling

The question: is "quiver mutation producing a tournament" equivalent to
some known operation (e.g., relabeling, or rotation of a regular tournament)?

Instance: kind-pasteur-2026-03-07-S34
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations
from math import comb
import random


def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def adj_to_bits(adj, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                bits |= (1 << idx)
            idx += 1
    return bits


def quiver_mutate(adj, n, k):
    new = [row[:] for row in adj]
    for i in range(n):
        if i == k:
            continue
        for j in range(n):
            if j == k or j == i:
                continue
            if adj[i][k] and adj[k][j]:
                new[i][j] += 1
    for i in range(n):
        if i == k:
            continue
        new[i][k], new[k][i] = new[k][i], new[i][k]
    for i in range(n):
        for j in range(i+1, n):
            cancel = min(new[i][j], new[j][i])
            new[i][j] -= cancel
            new[j][i] -= cancel
    return new


def is_tournament(adj, n):
    for i in range(n):
        for j in range(n):
            if i == j:
                if adj[i][j] != 0:
                    return False
            else:
                if adj[i][j] not in (0, 1) or adj[i][j] + adj[j][i] != 1:
                    return False
    return True


def count_ham_paths(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])


def score_seq(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))


def are_isomorphic(adj1, adj2, n):
    """Check if two tournaments are isomorphic (try all permutations for small n)."""
    s1 = score_seq(adj1, n)
    s2 = score_seq(adj2, n)
    if s1 != s2:
        return False
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(i+1, n):
                if adj1[i][j] != adj2[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False


def main():
    print("=== Quiver Mutation H-Invariance: Deep Analysis ===\n")

    for n in [4, 5]:
        print(f"\n--- n={n} (exhaustive) ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        total_mutations = 0
        tournament_results = 0
        h_preserved = 0
        isomorphic_count = 0
        self_mutation_count = 0  # mutation gives same tournament

        # Characterize when mutation produces tournament
        produces_tournament = {}  # (score_of_T, out_degree_of_k) -> count

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            scores = [sum(adj[i]) for i in range(n)]

            for k in range(n):
                total_mutations += 1
                mutated = quiver_mutate(adj, n, k)

                if is_tournament(mutated, n):
                    tournament_results += 1
                    H_mut = count_ham_paths(mutated, n)
                    if H_mut == H:
                        h_preserved += 1

                    # Check isomorphism
                    if are_isomorphic(adj, mutated, n):
                        isomorphic_count += 1

                    # Check if same tournament
                    if adj_to_bits(mutated, n) == bits:
                        self_mutation_count += 1

                    key = (tuple(sorted(scores)), scores[k])
                    produces_tournament[key] = produces_tournament.get(key, 0) + 1

        print(f"Total mutations: {total_mutations}")
        print(f"Produces tournament: {tournament_results}/{total_mutations} ({100*tournament_results/total_mutations:.1f}%)")
        print(f"H preserved: {h_preserved}/{tournament_results} ({100*h_preserved/max(1,tournament_results):.1f}%)")
        print(f"Isomorphic to original: {isomorphic_count}/{tournament_results} ({100*isomorphic_count/max(1,tournament_results):.1f}%)")
        print(f"Same tournament (identical): {self_mutation_count}/{tournament_results}")

        print(f"\nWhen does mutation produce a tournament?")
        for (sc, deg), cnt in sorted(produces_tournament.items()):
            print(f"  Score seq {sc}, mut at deg-{deg} vertex: {cnt} times")

    # n=7 sampling
    print(f"\n\n--- n=7 (sampling) ---")
    n = 7
    random.seed(42)
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    total = 0
    tournament_ok = 0
    h_preserved = 0
    h_violations = 0

    for trial in range(5000):
        bits = random.randint(0, max_bits - 1)
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)

        for k in range(n):
            total += 1
            mutated = quiver_mutate(adj, n, k)
            if is_tournament(mutated, n):
                tournament_ok += 1
                H_mut = count_ham_paths(mutated, n)
                if H_mut == H:
                    h_preserved += 1
                else:
                    h_violations += 1
                    print(f"  !!! VIOLATION at trial={trial}, k={k}: H={H} -> H_mut={H_mut}")

    print(f"Total: {total} mutations")
    print(f"Tournament: {tournament_ok}/{total} ({100*tournament_ok/total:.1f}%)")
    print(f"H preserved: {h_preserved}/{tournament_ok}")
    print(f"H violations: {h_violations}/{tournament_ok}")
    if h_violations == 0:
        print(f"*** CONFIRMED: Quiver mutation preserves H whenever result is a tournament ***")


if __name__ == "__main__":
    main()
