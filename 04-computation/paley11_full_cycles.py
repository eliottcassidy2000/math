#!/usr/bin/env python3
"""paley11_full_cycles.py -- Complete directed cycle count for Paley T_11.

Session: kind-pasteur-2026-03-20-S1

Uses Held-Karp DP on each vertex subset for efficiency.
For a subset of size k, counts directed Hamiltonian cycles by:
1. Fix smallest vertex as start
2. Run Held-Karp DP to find paths from start visiting all vertices in subset
3. Count paths ending at vertices adjacent back to start

Also computes the full alpha vector (independence polynomial of Omega).
"""

from itertools import combinations
from collections import Counter
from math import comb

def paley_adj(p):
    """Paley tournament adjacency matrix for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in qr:
            j = (i + s) % p
            adj[i][j] = 1
    return adj

def interval_adj(p):
    """Interval tournament adjacency matrix."""
    S = set(range(1, (p+1)//2))
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S:
            j = (i + s) % p
            adj[i][j] = 1
    return adj

def count_ham_cycles_on_subset(adj_global, verts):
    """Count directed Hamiltonian cycles on vertex subset using Held-Karp.
    Fix verts[0] as start, count cycles via DP.
    Returns number of distinct directed cycles.
    """
    k = len(verts)
    if k < 3:
        return 0

    # Build sub-adjacency
    sub = [[adj_global[verts[i]][verts[j]] for j in range(k)] for i in range(k)]

    # DP: dp[mask][v] = number of paths from vertex 0 to vertex v
    # visiting exactly the vertices in mask (0-indexed within subset)
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1  # Start at vertex 0, mask = {0}

    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]

    # Count cycles: paths visiting all vertices that return to vertex 0
    full_mask = (1 << k) - 1
    cycles = 0
    for v in range(1, k):  # v != 0
        if dp[full_mask][v] > 0 and sub[v][0]:
            cycles += dp[full_mask][v]

    return cycles


def compute_full_alpha(adj, n, directed_cycles_by_vset):
    """Compute alpha vector from directed cycle data.
    directed_cycles_by_vset: dict mapping frozenset -> count of directed cycles
    """
    # Create one Omega vertex per directed cycle
    # Two cycles conflict iff their vertex sets overlap
    # All directed cycles on same vertex set conflict with each other
    # and with any cycle whose vertex set overlaps

    # For alpha_k, we need independent sets of size k in Omega
    # An independent set is a collection of cycles with pairwise disjoint vertex sets

    # This is equivalent to counting collections of VERTEX-DISJOINT vertex sets,
    # weighted by the product of cycle counts on each set.

    # Formally: alpha_k = sum over k-tuples (S_1,...,S_k) of pairwise disjoint vertex sets
    #           product_{i=1}^k count[S_i]

    # At k=1: alpha_1 = sum of all cycle counts = total directed cycles
    # At k=2: alpha_2 = sum over disjoint pairs (S,T) of count[S]*count[T]

    vsets = list(directed_cycles_by_vset.keys())
    counts = [directed_cycles_by_vset[s] for s in vsets]
    num_vsets = len(vsets)

    # alpha_1
    alpha_1 = sum(counts)

    # alpha_2: pairs of disjoint vertex sets
    alpha_2 = 0
    for i in range(num_vsets):
        for j in range(i+1, num_vsets):
            if not (vsets[i] & vsets[j]):
                alpha_2 += counts[i] * counts[j]

    # alpha_3: triples of pairwise disjoint vertex sets
    alpha_3 = 0
    for i in range(num_vsets):
        for j in range(i+1, num_vsets):
            if vsets[i] & vsets[j]:
                continue
            for k in range(j+1, num_vsets):
                if (vsets[i] & vsets[k]) or (vsets[j] & vsets[k]):
                    continue
                alpha_3 += counts[i] * counts[j] * counts[k]

    return [1, alpha_1, alpha_2, alpha_3]


def main():
    p = 11
    print(f"Computing full directed cycle structure for Paley T_{p}...")
    adj = paley_adj(p)

    # Count directed Hamiltonian cycles for each odd-size vertex subset
    directed_by_vset = {}
    total_by_length = Counter()

    for size in range(3, p+1, 2):
        vset_count = 0
        dir_count = 0
        for verts in combinations(range(p), size):
            c = count_ham_cycles_on_subset(adj, verts)
            if c > 0:
                directed_by_vset[frozenset(verts)] = c
                vset_count += 1
                dir_count += c
                total_by_length[size] += c

        print(f"  {size}-cycles: {dir_count} directed cycles on {vset_count} vertex sets "
              f"(avg {dir_count/max(vset_count,1):.2f} per set)")

    total_directed = sum(total_by_length.values())
    print(f"\n  Total directed odd cycles: {total_directed}")
    print(f"  By length: {dict(sorted(total_by_length.items()))}")

    # Compute alpha vector
    print(f"\n  Computing alpha vector...")
    alpha = compute_full_alpha(adj, p, directed_by_vset)
    H = sum(a * 2**k for k, a in enumerate(alpha))
    print(f"  alpha = {alpha}")
    print(f"  I(Omega, 2) = {H}")
    print(f"  Expected H = 95095")
    print(f"  Match: {H == 95095}")

    # Also compute for interval tournament
    print(f"\n\nComputing for Interval T_{p}...")
    adj_int = interval_adj(p)

    int_by_vset = {}
    int_by_length = Counter()

    for size in range(3, p+1, 2):
        vset_count = 0
        dir_count = 0
        for verts in combinations(range(p), size):
            c = count_ham_cycles_on_subset(adj_int, verts)
            if c > 0:
                int_by_vset[frozenset(verts)] = c
                vset_count += 1
                dir_count += c
                int_by_length[size] += c

        print(f"  {size}-cycles: {dir_count} directed cycles on {vset_count} vertex sets "
              f"(avg {dir_count/max(vset_count,1):.2f} per set)")

    total_int = sum(int_by_length.values())
    print(f"\n  Total directed odd cycles: {total_int}")
    print(f"  By length: {dict(sorted(int_by_length.items()))}")

    alpha_int = compute_full_alpha(adj_int, p, int_by_vset)
    H_int = sum(a * 2**k for k, a in enumerate(alpha_int))
    print(f"  alpha = {alpha_int}")
    print(f"  I(Omega, 2) = {H_int}")
    print(f"  Expected H = 93027")
    print(f"  Match: {H_int == 93027}")

    # Comparison
    print(f"\n\n{'='*60}")
    print(f"COMPARISON: Paley vs Interval at p={p}")
    print(f"{'='*60}")
    print(f"  {'':>20} {'Paley':>12} {'Interval':>12} {'Ratio':>10}")
    for size in sorted(set(total_by_length.keys()) | set(int_by_length.keys())):
        pc = total_by_length.get(size, 0)
        ic = int_by_length.get(size, 0)
        ratio = pc/max(ic,1)
        print(f"  {size}-cycles:       {pc:>12} {ic:>12} {ratio:>10.3f}")
    print(f"  {'Total':>20} {total_directed:>12} {total_int:>12} {total_directed/max(total_int,1):>10.3f}")
    print(f"  {'alpha_1':>20} {alpha[1]:>12} {alpha_int[1]:>12}")
    print(f"  {'alpha_2':>20} {alpha[2]:>12} {alpha_int[2]:>12}")
    print(f"  {'alpha_3':>20} {alpha[3]:>12} {alpha_int[3]:>12}")
    print(f"  {'H':>20} {H:>12} {H_int:>12} {H/max(H_int,1):>10.6f}")


if __name__ == "__main__":
    main()
