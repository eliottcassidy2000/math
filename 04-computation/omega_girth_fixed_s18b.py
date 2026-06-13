#!/usr/bin/env python3
"""
omega_girth_fixed_s18b.py -- kind-pasteur-2026-03-21-S18b

FIXED girth computation for Omega(T).
The previous version had a BFS bug. This version uses correct cycle detection.
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, deque

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def find_all_odd_cycle_vertex_sets(A, n):
    cycle_sets = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            sub = list(subset)
            has_cycle = False
            for perm in permutations(sub[1:]):
                ordering = [sub[0]] + list(perm)
                is_cycle = True
                for idx in range(length):
                    if not A[ordering[idx]][ordering[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    has_cycle = True
                    break
            if has_cycle:
                cycle_sets.append(frozenset(subset))
    return cycle_sets

def build_omega_matrix(cycle_sets):
    """Build conflict graph as adjacency matrix."""
    nc = len(cycle_sets)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = True
    return adj, nc

def girth_from_matrix(adj, nc):
    """Compute girth using BFS from each vertex. Correct implementation."""
    if nc <= 2:
        # Need at least 3 vertices for a cycle
        # Check if nc=2 and there's an edge (that's not a cycle, just an edge)
        return float('inf')

    best = float('inf')

    for s in range(nc):
        # BFS from s
        dist = [-1] * nc
        parent = [-1] * nc
        dist[s] = 0
        queue = deque([s])

        while queue:
            u = queue.popleft()
            if dist[u] >= best // 2:
                break  # Can't find shorter cycle from this source
            for v in range(nc):
                if not adj[u][v]:
                    continue
                if dist[v] == -1:
                    # Tree edge
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    queue.append(v)
                elif parent[u] != v:
                    # Back/cross edge that isn't just the parent
                    # This forms a cycle of length dist[u] + dist[v] + 1
                    cycle_len = dist[u] + dist[v] + 1
                    if cycle_len < best:
                        best = cycle_len

    return best

def complement_matrix(adj, nc):
    """Build complement adjacency matrix."""
    comp = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if not adj[i][j]:
                comp[i][j] = comp[j][i] = True
    return comp

print("=" * 72)
print("  THE GIRTH OF OMEGA(T) — FIXED COMPUTATION")
print("  kind-pasteur-2026-03-21-S18b")
print("=" * 72)

# First: sanity check on known graphs
print("\n  SANITY CHECKS:")

# K_3 (complete on 3): girth should be 3
adj_k3 = [[False,True,True],[True,False,True],[True,True,False]]
print(f"  K_3 girth: {girth_from_matrix(adj_k3, 3)} (expected 3)")

# K_5: girth 3
adj_k5 = [[i!=j for j in range(5)] for i in range(5)]
print(f"  K_5 girth: {girth_from_matrix(adj_k5, 5)} (expected 3)")

# C_5 (5-cycle): girth 5
adj_c5 = [[False]*5 for _ in range(5)]
for i in range(5):
    adj_c5[i][(i+1)%5] = adj_c5[(i+1)%5][i] = True
print(f"  C_5 girth: {girth_from_matrix(adj_c5, 5)} (expected 5)")

# Path P_3: girth inf (tree)
adj_p3 = [[False,True,False],[True,False,True],[False,True,False]]
print(f"  P_3 girth: {girth_from_matrix(adj_p3, 3)} (expected inf)")

# Petersen: girth 5
pet_edges = [(0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,8),(4,9),(5,7),(5,8),(6,8),(6,9),(7,9)]
adj_pet = [[False]*10 for _ in range(10)]
for i,j in pet_edges:
    adj_pet[i][j] = adj_pet[j][i] = True
print(f"  Petersen girth: {girth_from_matrix(adj_pet, 10)} (expected 5)")

# ========================================================================
# MAIN COMPUTATION
# ========================================================================

for n in [3, 4, 5, 6]:
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")

    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    girth_to_H = defaultdict(list)
    girth_to_alpha1 = defaultdict(set)
    H_to_girth = defaultdict(set)
    anti_girth_to_H = defaultdict(list)

    total = 0
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1

        H = count_hp(A, n)
        cycles = find_all_odd_cycle_vertex_sets(A, n)
        nc = len(cycles)

        if nc <= 1:
            g = float('inf')
            ag = float('inf')
        else:
            adj, _ = build_omega_matrix(cycles)
            g = girth_from_matrix(adj, nc)

            # Anti-conflict girth
            comp = complement_matrix(adj, nc)
            has_comp_edge = any(comp[i][j] for i in range(nc) for j in range(i+1, nc))
            if not has_comp_edge:
                ag = float('inf')
            else:
                ag = girth_from_matrix(comp, nc)

        g_key = g if g != float('inf') else 'inf'
        ag_key = ag if ag != float('inf') else 'inf'

        girth_to_H[g_key].append(H)
        girth_to_alpha1[g_key].add(nc)
        H_to_girth[H].add(g_key)
        anti_girth_to_H[ag_key].append(H)
        total += 1

    print(f"\n  GIRTH OF OMEGA DISTRIBUTION:")
    for g in sorted(girth_to_H.keys(), key=lambda x: (float('inf') if x=='inf' else x)):
        Hs = girth_to_H[g]
        H_set = sorted(set(Hs))
        a1s = sorted(girth_to_alpha1[g])
        H_str = str(H_set) if len(H_set) <= 12 else f"{H_set[:5]}...{H_set[-3:]}"
        print(f"    girth={str(g):<5s}: count={len(Hs):<6d} H in {H_str}, alpha_1 in {a1s}")

    print(f"\n  H -> POSSIBLE GIRTHS:")
    for H_val in sorted(H_to_girth.keys()):
        girths = sorted(H_to_girth[H_val], key=lambda x: (float('inf') if x=='inf' else x))
        print(f"    H={H_val:<3d}: girths = {girths}")

    print(f"\n  ANTI-CONFLICT GIRTH (complement of Omega):")
    for g in sorted(anti_girth_to_H.keys(), key=lambda x: (float('inf') if x=='inf' else x)):
        Hs = anti_girth_to_H[g]
        H_set = sorted(set(Hs))
        H_str = str(H_set) if len(H_set) <= 10 else f"{H_set[:4]}...{H_set[-2:]}"
        print(f"    anti-girth={str(g):<5s}: count={len(Hs):<6d} H in {H_str}")

# ========================================================================
# SUMMARY
# ========================================================================
print(f"\n{'='*72}")
print(f"  SYNTHESIS: GIRTH OF OMEGA AND THE SIX PATTERNS")
print(f"{'='*72}")
