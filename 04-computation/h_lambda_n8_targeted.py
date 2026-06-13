#!/usr/bin/env python3
"""
h_lambda_n8_targeted.py — opus-2026-03-13-S71c

At n=8, 20k samples gave 0 ambiguous lambda fibers but almost no collisions.
Need to find lambda COLLISIONS first, then check if H differs.

Strategy: generate pairs via Vitali atom (1,1,2,2) reversal,
which preserves lambda by construction. Compare H values.
"""

import sys, time
import numpy as np
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

from itertools import combinations

def vitali_reversals(A, n):
    """Find all (1,1,2,2) reversals that preserve lambda.
    A Vitali atom reverses edges (a→b, c→d → d→c, b→a) where
    {a,b,c,d} forms a specific pattern in the lambda graph.

    Actually: ANY (1,1,2,2) reversal on a 4-vertex set preserves lambda.
    The Vitali reversal picks {i,j,k,l} with specific score distribution
    within the subgraph (scores 1,1,2,2 in the sub-tournament).
    """
    results = []
    for combo in combinations(range(n), 4):
        # Extract 4-vertex subtournament
        verts = list(combo)
        sub_scores = []
        for v in verts:
            s = sum(A[v][w] for w in verts if w != v)
            sub_scores.append(s)

        sorted_scores = sorted(sub_scores)
        if sorted_scores == [0, 1, 2, 3]:
            # Transitive sub-tournament — reversing creates (1,1,2,2)
            # The (1,1,2,2) reversal: reverse the unique HP direction
            # Actually: find the vertex with score 0 and score 3
            src = verts[sub_scores.index(3)]  # beats all 3
            sink = verts[sub_scores.index(0)]  # loses to all 3
            # Reverse edges src→sink only? No, (1,1,2,2) means reverse
            # a specific subset. Let me just reverse ALL edges to get (3,2,1,0).
            # That changes too many things.
            pass
        elif sorted_scores == [1, 1, 2, 2]:
            # Already a (1,1,2,2) subtournament
            # Reverse it to get another (1,1,2,2) or (0,1,2,3)?
            # The Vitali reversal: reverse ALL edges within this 4-set
            A_new = A.copy()
            for i_idx in range(4):
                for j_idx in range(i_idx+1, 4):
                    u, v = verts[i_idx], verts[j_idx]
                    A_new[u][v], A_new[v][u] = A_new[v][u], A_new[u][v]

            # Check if this preserves lambda
            results.append((combo, A_new))

    return results

n = 8
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: Searching for H-varying Vitali pairs...")
h_diff_count = 0
h_same_count = 0
dc7_dist = defaultdict(int)

t0 = time.time()
for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H_orig = count_ham_paths(A, n)

    # Find (1,1,2,2) subtournaments and reverse them
    for combo in combinations(range(n), 4):
        verts = list(combo)
        sub_scores = []
        for v in verts:
            s = sum(A[v][w] for w in verts if w != v)
            sub_scores.append(s)

        if sorted(sub_scores) != [1, 1, 2, 2]:
            continue

        # Reverse all edges within this 4-set
        A_new = A.copy()
        for i_idx in range(4):
            for j_idx in range(i_idx+1, 4):
                u, v = verts[i_idx], verts[j_idx]
                A_new[u][v], A_new[v][u] = A_new[v][u], A_new[u][v]

        H_new = count_ham_paths(A_new, n)
        dH = H_new - H_orig

        if dH != 0:
            h_diff_count += 1
            dc7_dist[abs(dH)] += 1
            if h_diff_count <= 10:
                c3 = int(np.trace(A @ A @ A)) // 3
                c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
                scores = sorted([int(sum(A[i])) for i in range(n)])
                print(f"  trial {trial}, combo={combo}: H={H_orig}→{H_new}, ΔH={dH}, c3={c3}, c5={c5}, scores={scores}")
        else:
            h_same_count += 1

        break  # Only first (1,1,2,2) per tournament

    if trial % 500 == 0:
        dt = time.time() - t0
        print(f"  ... trial {trial}: {dt:.1f}s, dH≠0: {h_diff_count}, dH=0: {h_same_count}")

dt = time.time() - t0
print(f"\n  Total: dH≠0: {h_diff_count}, dH=0: {h_same_count}, {dt:.1f}s")
print(f"  |ΔH| distribution: {dict(sorted(dc7_dist.items()))}")
print(f"  Note: ΔH = 2·Δc7 + 4·Δα₂(3,5) since at n=8, α₂ can include (3,5) pairs")

print("\nDone.")
