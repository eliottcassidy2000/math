#!/usr/bin/env python3
"""
opus-2026-05-29-S8

Audit a concrete n=8 counterexample to the false universal claim H(T) != 63.

The adjacency matrix is Example 1 from 05-knowledge/results/h63_verify.out.
This script independently checks:
  1. H(T) by Held-Karp DP
  2. H(T) by direct permutation enumeration
  3. I(Omega(T), 2) using directed odd cycles as Omega vertices
  4. basic structure of Omega(T), including connected components
"""

from __future__ import annotations

import itertools
from collections import Counter, deque


ADJ = [
    [0, 1, 0, 0, 0, 1, 0, 1],
    [0, 0, 1, 1, 1, 1, 1, 1],
    [1, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 1],
    [1, 0, 0, 1, 1, 1, 0, 1],
    [0, 0, 0, 1, 1, 0, 0, 0],
]


def count_hamiltonian_paths_dp(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full])


def count_hamiltonian_paths_direct(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for perm in itertools.permutations(range(n)):
        if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
            total += 1
    return total


def directed_odd_cycles(adj: list[list[int]]) -> list[tuple[int, ...]]:
    n = len(adj)
    cycles: list[tuple[int, ...]] = []
    for length in range(3, n + 1, 2):
        for perm in itertools.permutations(range(n), length):
            if perm[0] != min(perm):
                continue
            if all(adj[perm[i]][perm[(i + 1) % length]] for i in range(length)):
                cycles.append(perm)
    return cycles


def conflict_graph(cycles: list[tuple[int, ...]]) -> list[set[int]]:
    vertices = [set(c) for c in cycles]
    graph = [set() for _ in cycles]
    for i in range(len(cycles)):
        for j in range(i + 1, len(cycles)):
            if vertices[i] & vertices[j]:
                graph[i].add(j)
                graph[j].add(i)
    return graph


def independence_counts(graph: list[set[int]]) -> list[int]:
    # At n=8, at most two vertex-disjoint odd cycles can coexist, but keep this generic.
    counts = [1]
    for k in range(1, len(graph) + 1):
        count = 0
        for combo in itertools.combinations(range(len(graph)), k):
            if all(j not in graph[i] for i, j in itertools.combinations(combo, 2)):
                count += 1
        if count == 0:
            break
        counts.append(count)
    return counts


def connected_components(graph: list[set[int]]) -> list[list[int]]:
    unseen = set(range(len(graph)))
    comps: list[list[int]] = []
    while unseen:
        start = unseen.pop()
        comp = [start]
        q = deque([start])
        while q:
            v = q.popleft()
            for w in graph[v]:
                if w in unseen:
                    unseen.remove(w)
                    comp.append(w)
                    q.append(w)
        comps.append(sorted(comp))
    return sorted(comps, key=lambda c: (-len(c), c[0]))


def component_is_complete(graph: list[set[int]], comp: list[int]) -> bool:
    return all(j in graph[i] for i, j in itertools.combinations(comp, 2))


def strongly_connected_component_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reachable_from(source: int, reverse: bool = False) -> set[int]:
        seen = {source}
        q = deque([source])
        while q:
            v = q.popleft()
            for w in range(n):
                edge = adj[w][v] if reverse else adj[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        v = next(iter(remaining))
        comp = reachable_from(v) & reachable_from(v, reverse=True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def main() -> None:
    n = len(ADJ)

    h_dp = count_hamiltonian_paths_dp(ADJ)
    h_direct = count_hamiltonian_paths_direct(ADJ)
    cycles = directed_odd_cycles(ADJ)
    graph = conflict_graph(cycles)
    alpha = independence_counts(graph)
    ip_at_2 = sum((2**k) * a for k, a in enumerate(alpha))
    comps = connected_components(graph)
    score_sequence = sorted(sum(row) for row in ADJ)
    cycle_lengths = Counter(map(len, cycles))

    print("H=63 counterexample audit (opus-2026-05-29-S8)")
    print("=" * 60)
    print(f"n = {n}")
    print(f"score sequence = {score_sequence}")
    print(f"SCC sizes = {strongly_connected_component_sizes(ADJ)}")
    print()
    print("Adjacency matrix:")
    for row in ADJ:
        print("  " + " ".join(map(str, row)))
    print()
    print(f"H(T) by DP = {h_dp}")
    print(f"H(T) by direct permutation enumeration = {h_direct}")
    print(f"Counts agree = {h_dp == h_direct}")
    print()
    print(f"directed odd cycles |Omega vertices| = {len(cycles)}")
    print(f"cycle length distribution = {dict(sorted(cycle_lengths.items()))}")
    print(f"independence counts alpha_k = {alpha}")
    print(f"I(Omega(T), 2) = {ip_at_2}")
    print(f"OCF agrees with H = {ip_at_2 == h_dp}")
    print()
    print(f"Omega connected components = {[len(c) for c in comps]}")
    for idx, comp in enumerate(comps, 1):
        edges = sum(1 for i, j in itertools.combinations(comp, 2) if j in graph[i])
        complete = component_is_complete(graph, comp)
        lengths = Counter(len(cycles[i]) for i in comp)
        print(
            f"  component {idx}: size={len(comp)}, edges={edges}, "
            f"complete={complete}, cycle_lengths={dict(sorted(lengths.items()))}"
        )
    print()
    print("Conclusion:")
    print("  This concrete tournament has H(T)=63 at n=8.")
    print("  Therefore any universal theorem H(T) != 63 is false.")
    print("  The true statement is only that H=63 is absent in the n<=7 data.")


if __name__ == "__main__":
    main()
