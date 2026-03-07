#!/usr/bin/env python3
"""
Fast typed independence polynomial at n=7.

Uses DFS-based cycle finding instead of brute force permutations.

opus-2026-03-07-S36
"""
from itertools import combinations
import random
from collections import defaultdict

def tournament_from_seed(n, seed):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_directed_cycles_dfs(A, n, vertices=None):
    """Find all directed cycles using DFS. Returns list of canonical tuples."""
    if vertices is None:
        vertices = list(range(n))
    vset = set(vertices)
    cycles = set()

    for length in range(3, len(vertices)+1, 2):  # odd only
        for subset in combinations(vertices, length):
            sset = set(subset)
            # Try to find directed Hamiltonian cycles on this subset
            # DFS from smallest vertex
            start = subset[0]

            def dfs(path, visited):
                if len(path) == length:
                    # Check if last -> first
                    if A[path[-1]][start]:
                        canon = tuple(path)
                        # Canonical: rotate so min is first
                        min_idx = canon.index(min(canon))
                        canon = canon[min_idx:] + canon[:min_idx]
                        cycles.add(canon)
                    return
                last = path[-1]
                for nxt in subset:
                    if nxt not in visited and A[last][nxt]:
                        visited.add(nxt)
                        path.append(nxt)
                        dfs(path, visited)
                        path.pop()
                        visited.remove(nxt)

            dfs([start], {start})

    return list(cycles)

def I_at_2_fast(cycles_list):
    """Independence polynomial at x=2 via divide-and-conquer."""
    nc = len(cycles_list)
    if nc == 0:
        return 1
    cvsets = [frozenset(c) for c in cycles_list]
    adj = {}
    for i in range(nc):
        adj[i] = frozenset(j for j in range(nc) if j != i and cvsets[i] & cvsets[j])

    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return 1
        v = max(verts, key=lambda u: len(adj[u] & verts))
        without = solve(verts - {v})
        with_v = solve(verts - (adj[v] & verts) - {v})
        result = without + 2 * with_v
        memo[verts] = result
        return result

    return solve(frozenset(range(nc)))

def main():
    print("=== Typed Independence Polynomial at n=7 ===\n")

    n = 7
    for seed in range(10):
        A = tournament_from_seed(n, seed)
        cycles = find_directed_cycles_dfs(A, n)
        nc = len(cycles)
        cvsets = [frozenset(c) for c in cycles]

        # Classify by length
        by_len = defaultdict(list)
        for i, c in enumerate(cycles):
            by_len[len(c)].append(i)

        # Build adjacency
        adj_set = {}
        for i in range(nc):
            adj_set[i] = set()
            for j in range(nc):
                if i != j and cvsets[i] & cvsets[j]:
                    adj_set[i].add(j)

        # Count typed independent sets up to size 3
        typed_counts = defaultdict(int)

        # Size 1
        for i in range(nc):
            typed_counts[(len(cycles[i]),)] += 1

        # Size 2
        for i in range(nc):
            for j in range(i+1, nc):
                if j not in adj_set[i]:
                    key = tuple(sorted([len(cycles[i]), len(cycles[j])]))
                    typed_counts[key] += 1

        # Size 3
        for i in range(nc):
            for j in range(i+1, nc):
                if j in adj_set[i]:
                    continue
                for k in range(j+1, nc):
                    if k in adj_set[i] or k in adj_set[j]:
                        continue
                    key = tuple(sorted([len(cycles[i]), len(cycles[j]), len(cycles[k])]))
                    typed_counts[key] += 1

        H = I_at_2_fast(cycles)

        print(f"seed={seed}: {len(by_len[3])}x3 + {len(by_len[5])}x5 + {len(by_len[7])}x7 = {nc} cycles, H={H}")

        H_from_typed = 1 + sum(2**len(key) * typed_counts[key] for key in typed_counts)
        print(f"  H from typed: {H_from_typed} (match: {H_from_typed == H})")

        for key in sorted(typed_counts.keys()):
            S = sum(l - 1 for l in key)
            print(f"  type {key}: count={typed_counts[key]}, S={S}")

        # Check which types have same S
        s_groups = defaultdict(list)
        for key in typed_counts:
            S = sum(l - 1 for l in key)
            s_groups[S].append((key, typed_counts[key]))
        for S in sorted(s_groups.keys()):
            if len(s_groups[S]) > 1:
                print(f"  ** S={S} shared by: {[(k,c) for k,c in s_groups[S]]}")
        print()

if __name__ == "__main__":
    main()
