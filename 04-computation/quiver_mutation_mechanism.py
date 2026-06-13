#!/usr/bin/env python3
"""
Deep analysis of WHY quiver mutation preserves H(T).

Key question: When quiver mutation at vertex k produces a tournament,
what is the structural relationship between T and mu_k(T)?

Hypotheses to test:
1. Is mu_k(T) always T^op restricted somehow?
2. Does mu_k preserve the score sequence?
3. Does mu_k preserve the conflict graph Omega(T)?
4. Is mu_k related to the swap involution?
5. What is the graph-theoretic effect of mu_k on the arcs?

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
    return tuple(sum(adj[i]) for i in range(n))


def count_3cycles(adj, n):
    count = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            count += 1
    return count


def arc_diff(adj1, adj2, n):
    """Count arcs that differ between two tournaments."""
    diff = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj1[i][j] != adj2[i][j]:
                diff += 1
    return diff


def main():
    print("=== Quiver Mutation Mechanism Analysis ===\n")

    for n in [4, 5]:
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            scores_orig = score_seq(adj, n)
            c3_orig = count_3cycles(adj, n)

            for k in range(n):
                mutated = quiver_mutate(adj, n, k)
                if not is_tournament(mutated, n):
                    continue

                H_mut = count_ham_paths(mutated, n)
                scores_mut = score_seq(mutated, n)
                c3_mut = count_3cycles(mutated, n)
                n_diff = arc_diff(adj, mutated, n)

                # Analyze what changed
                # Which arcs flipped?
                flipped = []
                for i in range(n):
                    for j in range(i+1, n):
                        if adj[i][j] != mutated[i][j]:
                            if adj[i][j]:
                                flipped.append(f"{i}->{j} becomes {j}->{i}")
                            else:
                                flipped.append(f"{j}->{i} becomes {i}->{j}")

                # Score of vertex k
                deg_k_before = scores_orig[k]
                deg_k_after = scores_mut[k]

                if n == 4 or (n == 5 and bits < 32):
                    print(f"  bits={bits:0{num_edges}b} k={k}: "
                          f"H={H}->{H_mut} "
                          f"c3={c3_orig}->{c3_mut} "
                          f"scores={scores_orig}->{scores_mut} "
                          f"deg_k={deg_k_before}->{deg_k_after} "
                          f"#flipped={n_diff}")
                    if n_diff <= 4:
                        print(f"    Flipped arcs: {flipped}")

        # Summary statistics at n=5
        if n == 5:
            print(f"\n  Summary for n={n}:")
            score_changes = {}
            c3_changes = {}
            flip_counts = {}
            for bits in range(max_bits):
                adj = tournament_adj(n, bits)
                for k in range(n):
                    mutated = quiver_mutate(adj, n, k)
                    if not is_tournament(mutated, n):
                        continue
                    nd = arc_diff(adj, mutated, n)
                    sc_orig = sorted(score_seq(adj, n))
                    sc_mut = sorted(score_seq(mutated, n))
                    c3o = count_3cycles(adj, n)
                    c3m = count_3cycles(mutated, n)
                    flip_counts[nd] = flip_counts.get(nd, 0) + 1
                    score_changes[(tuple(sc_orig), tuple(sc_mut))] = \
                        score_changes.get((tuple(sc_orig), tuple(sc_mut)), 0) + 1
                    c3_changes[(c3o, c3m)] = c3_changes.get((c3o, c3m), 0) + 1

            print(f"  Arc flip counts: {dict(sorted(flip_counts.items()))}")
            print(f"  Score sequence changes:")
            for (so, sm), cnt in sorted(score_changes.items()):
                print(f"    {so} -> {sm}: {cnt}")
            print(f"  c3 changes: {dict(sorted(c3_changes.items()))}")


if __name__ == "__main__":
    main()
