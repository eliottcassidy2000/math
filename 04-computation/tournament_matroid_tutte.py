#!/usr/bin/env python3
"""
Tournament Matroid and Tutte Polynomial Analysis (kind-pasteur-S34).

Novel connection: Define a matroid on the arcs of a tournament.

Matroid 1 (Acyclic matroid):
  Ground set E = arcs of T
  Independent sets = subsets of arcs containing no directed cycle
  This is NOT a matroid in general (not closed under subset restriction
  in the right way), but the GRAPHIC matroid of K_n ignoring orientation
  gives a valid matroid. The Tutte polynomial T(x,y) at (2,0) counts
  acyclic orientations.

Matroid 2 (Transversal matroid):
  The bipartite graph with rows = positions (1..n) and columns = vertices
  has edge (i,v) iff vertex v can be placed at position i in some
  Hamiltonian path of T. This gives a transversal matroid.

We test:
1. Does the Tutte polynomial of the cycle matroid of K_n relate to H(T)?
2. What invariant of the tournament does T(K_n; 2,0) give?
3. Can we define a tournament-specific matroid whose Tutte polynomial gives H?

Instance: kind-pasteur-2026-03-07-S34
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
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


def has_directed_cycle(adj, n, arc_set):
    """Check if subset of arcs contains a directed cycle."""
    # Build subgraph with only the given arcs
    sub = [[0]*n for _ in range(n)]
    for (i, j) in arc_set:
        sub[i][j] = 1
    # Check for cycles via DFS
    WHITE, GRAY, BLACK = 0, 1, 2
    color = [WHITE] * n
    def dfs(v):
        color[v] = GRAY
        for u in range(n):
            if sub[v][u]:
                if color[u] == GRAY:
                    return True
                if color[u] == WHITE and dfs(u):
                    return True
        color[v] = BLACK
        return False
    for v in range(n):
        if color[v] == WHITE:
            if dfs(v):
                return True
    return False


def count_acyclic_subsets(adj, n):
    """Count subsets of arcs that are acyclic (contain no directed cycle)."""
    arcs = []
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                arcs.append((i, j))
    m = len(arcs)
    count = 0
    for mask in range(1 << m):
        arc_set = [arcs[k] for k in range(m) if mask & (1 << k)]
        if not has_directed_cycle(adj, n, arc_set):
            count += 1
    return count


def count_acyclic_subsets_by_size(adj, n):
    """Count acyclic subsets by size."""
    arcs = []
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                arcs.append((i, j))
    m = len(arcs)
    counts = [0] * (m + 1)
    for mask in range(1 << m):
        size = bin(mask).count('1')
        arc_set = [arcs[k] for k in range(m) if mask & (1 << k)]
        if not has_directed_cycle(adj, n, arc_set):
            counts[size] += 1
    return counts


def main():
    print("=== Tournament Matroid / Tutte Polynomial Analysis ===\n")

    # Small n only (exponential in arcs)
    for n in [3, 4]:
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            scores = tuple(sorted(sum(adj[i]) for i in range(n)))

            # Count acyclic subsets (no directed cycle)
            acyclic_count = count_acyclic_subsets(adj, n)
            acyclic_by_size = count_acyclic_subsets_by_size(adj, n)

            # Count c3 (directed 3-cycles)
            c3 = 0
            for vs in combinations(range(n), 3):
                a, b, c = vs
                if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                   (adj[a][c] and adj[c][b] and adj[b][a]):
                    c3 += 1

            print(f"  bits={bits:0{num_edges}b} scores={scores} H={H} c3={c3} "
                  f"#acyclic_subsets={acyclic_count} "
                  f"by_size={acyclic_by_size[:6]}")

    # Check if acyclic_count correlates with H
    print(f"\n\n--- n=4: H vs acyclic subset count ---")
    n = 4
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_acyclic = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        ac = count_acyclic_subsets(adj, n)
        if H not in h_to_acyclic:
            h_to_acyclic[H] = set()
        h_to_acyclic[H].add(ac)

    for H in sorted(h_to_acyclic):
        print(f"  H={H}: acyclic counts = {sorted(h_to_acyclic[H])}")


if __name__ == "__main__":
    main()
