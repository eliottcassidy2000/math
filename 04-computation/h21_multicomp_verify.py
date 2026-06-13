#!/usr/bin/env python3
"""
Verify multi-component factorizations at n=7.
Check that H=9 comes from (3,3) and H=15 from (3,5).

Also verify at n=5,6 for completeness.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict
import time


def held_karp(n, adj_bits):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj_bits[v] & (1 << u):
                    dp[S | (1 << u)][u] += c
    return sum(dp[(1 << n) - 1])


def get_three_cycles(n, adj_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] & (1 << j)):
                A[i][j] = 1
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append(frozenset({a, b, c}))
    return cycles


def get_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [set() for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        queue = [start]
        visited[start] = True
        while queue:
            v = queue.pop()
            comp.append(v)
            for u in adj[v]:
                if not visited[u]:
                    visited[u] = True
                    queue.append(u)
        components.append(comp)
    return components


def analyze_multicomp(n):
    """Find multi-component tournaments and analyze their structure."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    total = 1 << m

    multi_examples = defaultdict(list)  # H -> list of (bits, comp_cycle_sets)
    single_comp_H = set()

    t0 = time.time()
    for bits in range(total):
        adj_bits = [0]*n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj_bits[j] |= (1 << i)
            else:
                adj_bits[i] |= (1 << j)

        cycles = get_three_cycles(n, adj_bits)
        if not cycles:
            continue

        comps = get_components(cycles)
        if len(comps) > 1:
            H = held_karp(n, adj_bits)
            comp_info = []
            for comp in comps:
                cycle_vsets = [cycles[i] for i in comp]
                # All vertices used by this component
                all_verts = set()
                for c in cycle_vsets:
                    all_verts |= c
                comp_info.append((len(comp), all_verts))
            if len(multi_examples[H]) < 3:
                multi_examples[H].append((bits, comp_info))
        else:
            H = held_karp(n, adj_bits)
            single_comp_H.add(H)

    elapsed = time.time() - t0
    print(f"n={n} ({elapsed:.1f}s):")
    print(f"  Single-component H values: {sorted(single_comp_H)}")
    print(f"  Multi-component H values: {sorted(multi_examples.keys())}")

    for H in sorted(multi_examples.keys()):
        exs = multi_examples[H]
        print(f"\n  H={H}: {len(exs)} examples (showing first)")
        for bits, comp_info in exs[:2]:
            comp_desc = [(ncycles, sorted(verts)) for ncycles, verts in comp_info]
            print(f"    bits={bits}: components = {comp_desc}")

            # For each component, check: are the vertices disjoint from other components?
            all_comp_verts = [verts for _, verts in comp_info]
            for i in range(len(all_comp_verts)):
                for j in range(i+1, len(all_comp_verts)):
                    overlap = all_comp_verts[i] & all_comp_verts[j]
                    if overlap:
                        print(f"      WARNING: components {i},{j} share vertices {overlap}")
                    else:
                        print(f"      Components {i},{j} are vertex-disjoint")

    return multi_examples


if __name__ == '__main__':
    for n in [5, 6]:
        analyze_multicomp(n)
        print()
