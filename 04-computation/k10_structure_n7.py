#!/usr/bin/env python3
"""
Analyze K_10 structure at n=7: why do 10 pairwise-sharing odd cycles
always have exactly 2 disjoint pairs (i_2=2)?

At n=7, alpha_1=10 always gives H=29 = 1 + 20 + 8 = 1 + 2*10 + 4*2.
So i_2=2 always. This means among the 10 cycles, there are always
exactly 2 vertex-disjoint pairs.

This script examines the structure of these 10-cycle configurations.

Instance: opus-2026-03-07-S40
"""

from itertools import combinations
from collections import Counter
import time


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)
        yield adj


def find_directed_3_cycles(adj, n):
    """Find all directed 3-cycles as (vertex_set, direction) pairs."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i] >> j) & 1 and (adj[j] >> k) & 1 and (adj[k] >> i) & 1:
                    cycles.append(frozenset([i, j, k]))
                elif (adj[i] >> k) & 1 and (adj[k] >> j) & 1 and (adj[j] >> i) & 1:
                    cycles.append(frozenset([i, j, k]))
    return cycles


def find_directed_5_cycles(adj, n):
    """Find all directed 5-cycles (Hamiltonian cycles on 5-subsets)."""
    cycles = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S_mask in range(1, 32):
            for i in range(5):
                if not (S_mask & (1 << i)):
                    continue
                if (S_mask, i) not in dp:
                    continue
                c = dp[(S_mask, i)]
                for j in range(5):
                    if S_mask & (1 << j):
                        continue
                    if (adj[v[i]] >> v[j]) & 1:
                        key = (S_mask | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        # Count Hamiltonian cycles starting at v[0]
        for j in range(1, 5):
            if (31, j) in dp and (adj[v[j]] >> v[0]) & 1:
                count = dp[(31, j)]
                for _ in range(count):
                    cycles.append(frozenset(verts))
    return cycles


def find_directed_7_cycles(adj, n):
    """Find directed 7-cycles (Hamiltonian cycles on all 7 vertices)."""
    if n != 7:
        return []
    cycles = []
    dp = {}
    dp[(1, 0)] = 1
    for S in range(1, 128):
        for v in range(7):
            if not (S & (1 << v)):
                continue
            if (S, v) not in dp:
                continue
            c = dp[(S, v)]
            for u in range(7):
                if S & (1 << u):
                    continue
                if (adj[v] >> u) & 1:
                    key = (S | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    for j in range(1, 7):
        if (127, j) in dp and (adj[j] >> 0) & 1:
            count = dp[(127, j)]
            for _ in range(count):
                cycles.append(frozenset(range(7)))
    return cycles


def compute_i2(cycles):
    """Count number of vertex-disjoint pairs (i_2 in independence polynomial)."""
    m = len(cycles)
    i2 = 0
    for a in range(m):
        for b in range(a+1, m):
            if not (cycles[a] & cycles[b]):
                i2 += 1
    return i2


def main():
    n = 7
    print(f"Analyzing alpha_1=10 tournaments at n={n}...")
    start = time.time()

    results = []
    for adj in all_tournaments(n):
        c3 = find_directed_3_cycles(adj, n)
        c5 = find_directed_5_cycles(adj, n)
        c7 = find_directed_7_cycles(adj, n)
        all_cycles = c3 + c5 + c7
        alpha1 = len(all_cycles)

        if alpha1 != 10:
            continue

        i2 = compute_i2(all_cycles)
        t3 = len(c3)
        t5 = len(c5)
        t7 = len(c7)

        # Analyze which pairs are disjoint
        disjoint_pairs = []
        for a in range(len(all_cycles)):
            for b in range(a+1, len(all_cycles)):
                if not (all_cycles[a] & all_cycles[b]):
                    ca = all_cycles[a]
                    cb = all_cycles[b]
                    la = len(ca)
                    lb = len(cb)
                    disjoint_pairs.append((la, lb))

        results.append((t3, t5, t7, i2, sorted(disjoint_pairs)))

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s. Found {len(results)} tournaments with alpha_1=10.\n")

    # Summarize
    combo_dist = Counter((t3, t5, t7) for t3, t5, t7, _, _ in results)
    i2_dist = Counter(i2 for _, _, _, i2, _ in results)
    pair_dist = Counter(tuple(dp) for _, _, _, _, dp in results)

    print("=" * 60)
    print("alpha_1=10: (t3, t5, t7) distribution")
    print("=" * 60)
    for combo, cnt in sorted(combo_dist.items()):
        print(f"  {combo}: {cnt} tournaments")

    print(f"\ni_2 distribution: {dict(sorted(i2_dist.items()))}")

    print(f"\nDisjoint pair types (cycle lengths):")
    for pairs, cnt in sorted(pair_dist.items()):
        print(f"  {pairs}: {cnt}")

    # Check: are the 2 disjoint pairs always between specific cycle types?
    print("\n" + "=" * 60)
    print("Detail: what are the 2 disjoint pairs?")
    print("=" * 60)
    for r in results[:20]:  # show first 20
        t3, t5, t7, i2, dp = r
        print(f"  ({t3},{t5},{t7}), i2={i2}, disjoint pairs: {dp}")


if __name__ == "__main__":
    main()
