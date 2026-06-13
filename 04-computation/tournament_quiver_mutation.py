#!/usr/bin/env python3
"""
Tournament-Quiver Mutation Analysis (kind-pasteur-S34).

Novel connection: A tournament is a quiver (complete, no 2-cycles).
Cluster algebra mutations on quivers transform the adjacency matrix.
We investigate:

1. How does quiver mutation relate to arc-flip?
2. Is there a cluster variable that equals H(T) or I(Omega(T), x)?
3. What is the exchange graph for tournaments under mutation?
4. Does the Laurent phenomenon constrain H(T) values?

Quiver mutation at vertex k:
- For each pair i,j with arrows i->k and k->j: add arrow i->j
- Reverse all arrows incident to k
- Remove any 2-cycles created

For TOURNAMENTS (complete, no 2-cycles), mutation at k:
- Reverses all arcs incident to k (like an arc-flip of the vertex-star)
- The "add i->j" step creates 2-cycles which are then removed
- Net effect: complex interaction between vertex-k reversal and transitive closure
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from math import comb
import random


def tournament_adj(n, bits):
    """Generate tournament adjacency from bits."""
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


def quiver_mutate(adj, n, k):
    """Perform quiver mutation at vertex k on a tournament.

    Standard quiver mutation:
    1. For each i->k->j path (2-path through k): add arrow i->j
    2. Reverse all arrows incident to k
    3. Remove 2-cycles

    Returns the resulting digraph (may NOT be a tournament anymore).
    """
    # Copy
    new_adj = [row[:] for row in adj]

    # Step 1: Add arrows for 2-paths through k
    for i in range(n):
        if i == k:
            continue
        for j in range(n):
            if j == k or j == i:
                continue
            if adj[i][k] and adj[k][j]:
                new_adj[i][j] += 1  # may create 2-cycle

    # Step 2: Reverse all arrows incident to k
    for i in range(n):
        if i == k:
            continue
        old_ik = new_adj[i][k]
        old_ki = new_adj[k][i]
        new_adj[i][k] = old_ki
        new_adj[k][i] = old_ik

    # Step 3: Remove 2-cycles (cancel pairs)
    for i in range(n):
        for j in range(i+1, n):
            # Remove min(a[i][j], a[j][i]) from both
            cancel = min(new_adj[i][j], new_adj[j][i])
            new_adj[i][j] -= cancel
            new_adj[j][i] -= cancel

    return new_adj


def is_tournament(adj, n):
    """Check if adj is a valid tournament (complete, no 2-cycles, all 0/1)."""
    for i in range(n):
        for j in range(n):
            if i == j:
                if adj[i][j] != 0:
                    return False
            else:
                if adj[i][j] not in (0, 1):
                    return False
                if adj[i][j] + adj[j][i] != 1:
                    return False
    return True


def tournament_to_bits(adj, n):
    """Convert adjacency matrix to bit encoding."""
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                bits |= (1 << idx)
            idx += 1
    return bits


def count_ham_paths(adj, n):
    """Count Hamiltonian paths via DP (Held-Karp)."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def vertex_star_flip(adj, n, k):
    """Flip all arcs incident to vertex k (simpler operation)."""
    new_adj = [row[:] for row in adj]
    for i in range(n):
        if i == k:
            continue
        new_adj[i][k], new_adj[k][i] = new_adj[k][i], new_adj[i][k]
    return new_adj


def main():
    print("=== Tournament Quiver Mutation Analysis ===")
    print()

    for n in [4, 5, 6]:
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        if n <= 5:
            trials = range(max_bits)
            label = "exhaustive"
        else:
            random.seed(42)
            trials = [random.randint(0, max_bits - 1) for _ in range(2000)]
            label = "2000 random"

        still_tournament = 0
        not_tournament = 0
        h_preserved = 0
        h_changed = 0
        h_diffs = {}

        star_flip_h_preserved = 0
        star_flip_h_changed = 0

        total = 0
        for bits in trials:
            adj = tournament_adj(n, bits if isinstance(bits, int) else bits)
            H_orig = count_ham_paths(adj, n)

            for k in range(n):
                # Quiver mutation at k
                mutated = quiver_mutate(adj, n, k)
                if is_tournament(mutated, n):
                    still_tournament += 1
                    H_mut = count_ham_paths(mutated, n)
                    if H_mut == H_orig:
                        h_preserved += 1
                    else:
                        h_changed += 1
                        diff = H_mut - H_orig
                        h_diffs[diff] = h_diffs.get(diff, 0) + 1
                else:
                    not_tournament += 1

                # Vertex star flip at k (simpler: just reverse all k's arcs)
                star_flipped = vertex_star_flip(adj, n, k)
                H_star = count_ham_paths(star_flipped, n)
                if H_star == H_orig:
                    star_flip_h_preserved += 1
                else:
                    star_flip_h_changed += 1

            total += 1

        print(f"  {total} tournaments ({label})")
        print(f"  Quiver mutation results:")
        print(f"    Still tournament: {still_tournament}/{still_tournament + not_tournament}")
        print(f"    NOT tournament: {not_tournament}/{still_tournament + not_tournament}")
        if still_tournament > 0:
            print(f"    H preserved: {h_preserved}/{still_tournament}")
            print(f"    H changed: {h_changed}/{still_tournament}")
            if h_diffs:
                print(f"    H diffs (first 10): {dict(list(sorted(h_diffs.items()))[:10])}")

        print(f"  Star-flip results (reverse all arcs at vertex k):")
        total_flips = star_flip_h_preserved + star_flip_h_changed
        print(f"    H preserved: {star_flip_h_preserved}/{total_flips} ({100*star_flip_h_preserved/total_flips:.1f}%)")
        print(f"    H changed: {star_flip_h_changed}/{total_flips} ({100*star_flip_h_changed/total_flips:.1f}%)")


if __name__ == "__main__":
    main()
