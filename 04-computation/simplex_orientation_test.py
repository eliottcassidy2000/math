#!/usr/bin/env python3
"""
Tournament as Simplex Orientation (kind-pasteur-S34).

A tournament on n vertices is an orientation of K_n, the 1-skeleton
of the (n-1)-simplex Delta^{n-1}.

Novel ideas to test:
1. The "acyclicity number" = max acyclic subgraph = min feedback arc set
   This is the NUMBER OF FORWARD ARCS in some optimal linear ordering.
   The "min feedback arc set" fas(T) relates to H via:
   H = F_{n-1} (all forward paths) but also counts paths not respecting
   any particular ordering.

2. The "Euler characteristic" of the acyclic complex:
   The acyclic complex Acyc(T) = {subsets of arcs containing no directed cycle}
   is a simplicial complex. Its f-vector and Euler characteristic relate to
   the cycle structure.

3. Tournament flow polytope:
   Fix source s and sink t. The flow polytope P(T, s, t) consists of
   all unit flows from s to t in the tournament (nonneg weights on arcs,
   conservation at internal vertices). The volume of this polytope counts
   something related to paths from s to t.

4. Connection to Morse theory:
   A discrete Morse function on Delta^{n-1} assigns values to vertices
   compatible with the orientation. The number of critical cells relates
   to the Betti numbers of the complex.
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations, permutations
from math import comb


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


def min_feedback_arc_set(adj, n):
    """Compute the minimum feedback arc set size.
    = n*(n-1)/2 - max number of forward arcs over all linear orderings.
    Brute force over permutations for small n."""
    m = n * (n - 1) // 2
    max_forward = 0
    for perm in permutations(range(n)):
        forward = 0
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]]:
                    forward += 1
        max_forward = max(max_forward, forward)
    return m - max_forward


def count_acyclic_subsets(adj, n):
    """Count acyclic subsets of arcs (no directed cycle)."""
    arcs = [(i, j) for i in range(n) for j in range(n) if adj[i][j]]
    m = len(arcs)

    def has_cycle(arc_set):
        sub = [[0]*n for _ in range(n)]
        for (i, j) in arc_set:
            sub[i][j] = 1
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

    count = 0
    for mask in range(1 << m):
        arc_set = [arcs[k] for k in range(m) if mask & (1 << k)]
        if not has_cycle(arc_set):
            count += 1
    return count


def acyclic_euler_char(adj, n):
    """Compute Euler characteristic of acyclic complex.
    chi = sum_{k >= 0} (-1)^k * f_k where f_k = # acyclic subsets of size k."""
    arcs = [(i, j) for i in range(n) for j in range(n) if adj[i][j]]
    m = len(arcs)

    def has_cycle(arc_set):
        sub = [[0]*n for _ in range(n)]
        for (i, j) in arc_set:
            sub[i][j] = 1
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

    euler = 0
    for mask in range(1 << m):
        size = bin(mask).count('1')
        arc_set = [arcs[k] for k in range(m) if mask & (1 << k)]
        if not has_cycle(arc_set):
            euler += (-1)**size
    return euler


def main():
    print("=== Tournament as Simplex Orientation ===\n")

    # Part 1: Min feedback arc set vs H
    print("--- Part 1: Min Feedback Arc Set vs H ---")
    for n in [4, 5]:
        print(f"\n  n={n}:")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        fas_to_h = {}
        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            fas = min_feedback_arc_set(adj, n)

            if fas not in fas_to_h:
                fas_to_h[fas] = set()
            fas_to_h[fas].add(H)

        for fas in sorted(fas_to_h):
            h_vals = sorted(fas_to_h[fas])
            print(f"    fas={fas}: H values = {h_vals}")

    # Part 2: Acyclic complex at n=3,4
    print("\n\n--- Part 2: Acyclic Complex ---")
    for n in [3, 4]:
        print(f"\n  n={n}:")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            ac_count = count_acyclic_subsets(adj, n)
            ec = acyclic_euler_char(adj, n)
            scores = tuple(sorted(sum(adj[i]) for i in range(n)))
            print(f"    bits={bits:0{num_edges}b} scores={scores} H={H} "
                  f"#acyclic={ac_count} euler_chi={ec}")

    # Part 3: Is Euler char of acyclic complex related to H?
    print("\n\n--- Part 3: Euler char of acyclic complex vs H ---")
    n = 4
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_ec = {}
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        ec = acyclic_euler_char(adj, n)
        if H not in h_to_ec:
            h_to_ec[H] = set()
        h_to_ec[H].add(ec)

    for H in sorted(h_to_ec):
        vals = sorted(h_to_ec[H])
        print(f"  H={H}: euler_chi(Acyc) = {vals}")


if __name__ == "__main__":
    main()
