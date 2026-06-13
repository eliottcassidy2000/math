#!/usr/bin/env python3
"""
Why 5-cycles never bridge 3-cycle components.

THEOREM: If two 3-cycle components are vertex-disjoint in the 3-cycle
conflict graph, then no odd cycle uses vertices from both components.

PROOF SKETCH:
Let C1 = {v1, v2, v3} and C2 = {v4, v5, v6} be 3-cycles in different
components (disjoint vertex sets).

A 5-cycle using vertices from both must use some from {v1,v2,v3} and
some from {v4,v5,v6}. Say it uses a from C1 and b from C2, with a+b <= 5.

If it uses vertices from C1: those vertices form a sub-path of the 5-cycle.
The edges between C1 vertices must be consistent with the tournament.

Key: if the 5-cycle uses 2 vertices from C1 and 3 from C2, or vice versa,
then the 3 vertices from one group form a 3-element subset. If that subset
is a directed 3-cycle, it shares vertices with both C1/C2, connecting the
components -- contradiction since components are disjoint.

Actually, let me think about this more carefully. The 5-cycle vertex set
using vertices from both C1 and C2 doesn't necessarily mean the 3-vertex
subset forms a 3-cycle.

Let me check the actual multi-component tournaments to understand the
vertex structure.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict


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
    return cycles, A


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


def analyze_disconnected_structure():
    """Analyze the vertex structure of disconnected tournaments at n=6,7."""
    for n in [6, 7]:
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        m = len(edges)
        total = 1 << m

        comp_patterns = defaultdict(int)

        for bits in range(total):
            adj_bits = [0]*n
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    adj_bits[j] |= (1 << i)
                else:
                    adj_bits[i] |= (1 << j)

            three_cycles, A = get_three_cycles(n, adj_bits)
            if not three_cycles:
                continue

            comps = get_components(three_cycles)
            if len(comps) <= 1:
                continue

            comp_verts = []
            for comp in comps:
                verts = frozenset().union(*(three_cycles[i] for i in comp))
                comp_verts.append(verts)

            # Record: sizes of vertex sets, and whether they partition [n]
            vert_sizes = tuple(sorted(len(cv) for cv in comp_verts))
            all_used = frozenset().union(*comp_verts)
            unused = set(range(n)) - all_used
            pattern = (vert_sizes, len(unused))
            comp_patterns[pattern] += 1

        print(f"n={n}: disconnected tournament patterns (vertex sizes, unused count):")
        for (sizes, unused), cnt in sorted(comp_patterns.items()):
            print(f"  vertex sizes {sizes}, {unused} unused vertices: {cnt} tournaments")

        # Check: are the component vertex sets always disjoint?
        print(f"\n  Checking vertex disjointness:")
        overlap_count = 0
        for bits in range(total):
            adj_bits = [0]*n
            for k, (i, j) in enumerate(edges):
                if bits & (1 << k):
                    adj_bits[j] |= (1 << i)
                else:
                    adj_bits[i] |= (1 << j)

            three_cycles, A = get_three_cycles(n, adj_bits)
            if not three_cycles:
                continue

            comps = get_components(three_cycles)
            if len(comps) <= 1:
                continue

            comp_verts = []
            for comp in comps:
                verts = frozenset().union(*(three_cycles[i] for i in comp))
                comp_verts.append(verts)

            for i in range(len(comp_verts)):
                for j in range(i+1, len(comp_verts)):
                    if comp_verts[i] & comp_verts[j]:
                        overlap_count += 1

        print(f"  Overlapping component pairs: {overlap_count}")
        print()


if __name__ == '__main__':
    analyze_disconnected_structure()
