#!/usr/bin/env python3
"""
Transfer matrix structure at n=7 — sampled, since exhaustive is too expensive.
Focus on: eigenvalue patterns, sparsity, relationship to cycle counts.
"""
from itertools import permutations
import numpy as np
import random
from math import factorial

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def compute_M_diagonal(A, n):
    """Just the diagonal: M[a,a] = sum_P (-1)^{pos(a,P)} over HPs"""
    diag = [0] * n
    # Use DP to enumerate HPs
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = {v: 1}  # track which vertices are at which position

    # Actually, let's just track the sign contribution per vertex
    # M[a,a] = sum over HPs P of (-1)^{position of a in P}

    # Method: for each HP, find position of each vertex
    full = (1 << n) - 1

    # Reconstruct all HPs using DP
    # dp[(mask, last)] = count of paths ending at 'last' visiting 'mask'
    dp_count = {}
    for v in range(n):
        dp_count[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp_count:
                continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp_count[key] = dp_count.get(key, 0) + dp_count[(mask, v)]

    # For diagonal, we need sum_P (-1)^{pos(a,P)} for each a
    # This requires tracking position parity per vertex through DP

    # dp_sign[(mask, last, a_parity)] where a_parity tracks (-1)^{pos(a)}
    # Too complex for general n=7. Let me use a simpler approach.

    # For each vertex a, use the formula:
    # M[a,a] = sum_{S subset of V\{a}} (-1)^|S| * E_a(S) * B_a(V\{a}\S)
    # where E_a(S) = #HPs in T[S∪{a}] ending at a
    # B_a(R) = #HPs in T[{a}∪R] starting from a

    V = set(range(n))
    for a in range(n):
        others = sorted(V - {a})
        total = 0
        for mask in range(1 << len(others)):
            S = set()
            for idx in range(len(others)):
                if (mask >> idx) & 1:
                    S.add(others[idx])
            R = set(others) - S

            # E_a(S): HPs in T[S∪{a}] ending at a
            verts = sorted(S | {a})
            k = len(verts)
            E_a = count_HPs_ending_at(A, verts, a)

            # B_a(R): HPs in T[{a}∪R] starting from a
            verts2 = sorted({a} | R)
            B_a = count_HPs_starting_at(A, verts2, a)

            total += (-1)**len(S) * E_a * B_a
        diag[a] = total
    return diag

def count_HPs_ending_at(A, verts, target):
    """Count HPs through verts ending at target, using DP"""
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == target else 0

    idx_map = {v: i for i, v in enumerate(verts)}
    dp = {}
    for v in verts:
        if v == target and n > 1:
            continue  # target must be last, so don't start there unless n=1
        dp[(1 << idx_map[v], idx_map[v])] = 1

    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for vi in range(n):
            if not (mask & (1 << vi)) or (mask, vi) not in dp:
                continue
            for ui in range(n):
                if mask & (1 << ui): continue
                v, u = verts[vi], verts[ui]
                if A[v][u]:
                    key = (mask | (1 << ui), ui)
                    dp[key] = dp.get(key, 0) + dp[(mask, vi)]

    return dp.get((full, idx_map[target]), 0)

def count_HPs_starting_at(A, verts, source):
    """Count HPs through verts starting at source, using DP"""
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == source else 0

    idx_map = {v: i for i, v in enumerate(verts)}
    dp = {(1 << idx_map[source], idx_map[source]): 1}

    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for vi in range(n):
            if not (mask & (1 << vi)) or (mask, vi) not in dp:
                continue
            for ui in range(n):
                if mask & (1 << ui): continue
                v, u = verts[vi], verts[ui]
                if A[v][u]:
                    key = (mask | (1 << ui), ui)
                    dp[key] = dp.get(key, 0) + dp[(mask, vi)]

    return sum(dp.get((full, ui), 0) for ui in range(n))

random.seed(42)

print("=== Transfer matrix DIAGONAL at n=7 ===")
n = 7
for trial in range(5):
    A = random_tournament(n)
    H = ham_path_count_dp(A, n)
    diag = compute_M_diagonal(A, n)
    scores = [sum(A[i]) for i in range(n)]

    print(f"\n  trial {trial}: H={H}, sum(diag)={sum(diag)}")
    print(f"    scores: {scores}")
    print(f"    diag: {diag}")
    print(f"    diag range: [{min(diag)}, {max(diag)}]")
    print(f"    |diag| range: [{min(abs(d) for d in diag)}, {max(abs(d) for d in diag)}]")

    # Is sum(diag) = H?
    print(f"    sum(diag) == H? {sum(diag) == H}")
