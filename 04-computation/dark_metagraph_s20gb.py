#!/usr/bin/env python3
"""
dark_metagraph_s20gb.py — Enumerate odd graphs and compute their metagraph
kind-pasteur-2026-03-25-S20gb

Treat odd graphs as first-class objects. Build their metagraph
(edge flip connects two odd-graph classes). Compute all the same
metrics we compute for tournaments.

An ODD GRAPH = graph where at least one vertex has odd degree.
Equivalently: NOT an even graph (not all degrees even).

The "dark metagraph" connects two odd-graph classes if they differ
by toggling one edge. Compare with the tournament metagraph.
"""

import sys
import numpy as np
from math import factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DARK METAGRAPH: Odd Graphs as First-Class Objects")
print("  kind-pasteur-2026-03-25-S20gb")
print("=" * 80)

for n in [4, 5, 6]:
    t0 = time.time()
    N = n
    edges = [(i,j) for i in range(N) for j in range(i+1,N)]
    m = len(edges)  # C(n,2)
    all_perms = list(permutations(range(N)))

    def graph_canon(bits):
        """Canonical form of undirected graph."""
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = edges.index((pi, pj))
                if bits & (1 << k):
                    nb |= (1 << nk)
            if best is None or nb < best:
                best = nb
        return best

    def is_even(bits):
        """Check if all vertices have even degree."""
        deg = [0] * N
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                deg[i] += 1
                deg[j] += 1
        return all(d % 2 == 0 for d in deg)

    def degree_seq(bits):
        deg = [0] * N
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                deg[i] += 1
                deg[j] += 1
        return tuple(sorted(deg))

    def edge_count(bits):
        return bin(bits).count('1')

    # Enumerate ALL graphs
    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} edges")
    print(f"{'#'*60}")

    print(f"  Enumerating...", end=" ", flush=True)
    canon_map = {}
    even_map = {}
    for bits in range(1 << m):
        cn = graph_canon(bits)
        canon_map[bits] = cn
        even_map[bits] = is_even(bits)

    graph_classes = sorted(set(canon_map.values()))
    even_classes = sorted(set(cn for bits, cn in canon_map.items() if even_map[bits]))
    odd_classes = sorted(set(cn for bits, cn in canon_map.items() if not even_map[bits]))

    V_all = len(graph_classes)
    V_even = len(even_classes)
    V_odd = len(odd_classes)

    print(f"done. V_all={V_all}, V_even={V_even}, V_odd={V_odd}")

    # Build edge-flip metagraph for ODD graphs
    odd_edges = set()
    odd_sl = 0
    even_edges = set()
    cross_edges = set()  # even <-> odd

    for bits in range(1 << m):
        cn1 = canon_map[bits]
        is_odd1 = not even_map[bits]

        for k in range(m):
            flipped = bits ^ (1 << k)
            cn2 = canon_map[flipped]
            is_odd2 = not even_map[flipped]

            if is_odd1 and is_odd2 and cn1 != cn2:
                odd_edges.add((min(cn1,cn2), max(cn1,cn2)))
            elif is_odd1 and is_odd2 and cn1 == cn2:
                odd_sl += 1
            elif not is_odd1 and not is_odd2 and cn1 != cn2:
                even_edges.add((min(cn1,cn2), max(cn1,cn2)))
            elif is_odd1 != is_odd2:
                cross_edges.add((min(cn1,cn2), max(cn1,cn2)))

    E_odd = len(odd_edges)
    E_even = len(even_edges)
    E_cross = len(cross_edges)

    print(f"\n  METAGRAPH COMPARISON:")
    print(f"    {'':>15} {'V':>6} {'E':>6} {'E/V':>6}")
    print(f"    {'Even (tourn.)':>15} {V_even:6d} {E_even:6d} {E_even/max(V_even,1):6.1f}")
    print(f"    {'Odd (dark)':>15} {V_odd:6d} {E_odd:6d} {E_odd/max(V_odd,1):6.1f}")
    print(f"    {'Cross':>15} {'':>6} {E_cross:6d}")
    print(f"    {'Total graphs':>15} {V_all:6d}")

    # Properties of odd graph classes
    odd_class_masks = defaultdict(list)
    for bits, cn in canon_map.items():
        if not even_map[bits]:
            odd_class_masks[cn].append(bits)

    print(f"\n  ODD GRAPH CLASS PROPERTIES:")
    print(f"  {'Edges':>5} {'DegSeq':>20} {'Size':>6}")
    for cn in sorted(odd_classes, key=lambda c: edge_count(c)):
        bits = cn
        ec = edge_count(bits)
        ds = degree_seq(bits)
        sz = len(odd_class_masks[cn])
        if V_odd <= 30 or ec <= 3 or ec >= m-3:
            print(f"  {ec:5d} {str(ds):>20} {sz:6d}")

    # Adjacency spectrum for odd metagraph
    if V_odd <= 100 and V_odd > 0:
        odd_list = sorted(odd_classes)
        oidx = {cn: i for i, cn in enumerate(odd_list)}
        A_odd = np.zeros((V_odd, V_odd))
        for a, b in odd_edges:
            if a in oidx and b in oidx:
                A_odd[oidx[a], oidx[b]] = 1
                A_odd[oidx[b], oidx[a]] = 1

        if V_odd > 1:
            evals = sorted(np.linalg.eigvalsh(A_odd), reverse=True)
            deg = A_odd.sum(axis=1)
            print(f"\n  ODD METAGRAPH SPECTRUM:")
            print(f"    Eigenvalues (top 5): {[f'{x:.2f}' for x in evals[:5]]}")
            print(f"    Degree: min={int(deg.min())}, max={int(deg.max())}, avg={deg.mean():.1f}")

            # Connected?
            from collections import deque
            visited = {0}
            queue = deque([0])
            while queue:
                u = queue.popleft()
                for v in range(V_odd):
                    if A_odd[u][v] and v not in visited:
                        visited.add(v)
                        queue.append(v)
            print(f"    Connected: {len(visited) == V_odd}")

    # Compare even and odd metagraphs
    if V_even > 0 and V_odd > 0:
        even_list = sorted(even_classes)
        eidx = {cn: i for i, cn in enumerate(even_list)}
        A_even = np.zeros((V_even, V_even))
        for a, b in even_edges:
            if a in eidx and b in eidx:
                A_even[eidx[a], eidx[b]] = 1
                A_even[eidx[b], eidx[a]] = 1

        if V_even > 1:
            evals_e = sorted(np.linalg.eigvalsh(A_even), reverse=True)
            deg_e = A_even.sum(axis=1)
            print(f"\n  EVEN METAGRAPH SPECTRUM (for comparison):")
            print(f"    Eigenvalues (top 5): {[f'{x:.2f}' for x in evals_e[:5]]}")
            print(f"    Degree: min={int(deg_e.min())}, max={int(deg_e.max())}, avg={deg_e.mean():.1f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
