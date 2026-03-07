#!/usr/bin/env python3
"""
Search for non-circulant position-uniform tournaments at n=7.

Since full enumeration (2^21) is too expensive, we use two strategies:
1. Check specific construction families
2. Random sampling with efficient Ham path counting via DP

For position-uniformity check, we need N[v,j] for all v,j.
Ham path DP: dp[S][v] = number of Ham paths using vertex set S ending at v.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np
import random

def ham_path_dp(A):
    """Count Ham paths by DP on subsets. Returns N[v][j] matrix."""
    n = len(A)
    # dp[S][v] = number of paths through vertex set S ending at v
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[u][v] and (S_prev, u) in dp:
                        total += dp[(S_prev, u)]
                if total > 0:
                    dp[(S, v)] = total

    # Total paths and position distribution
    full = (1 << n) - 1
    H = sum(dp.get((full, v), 0) for v in range(n))
    return H

def position_matrix_dp(A):
    """Compute N[v][j] efficiently using layered DP."""
    n = len(A)
    # dp[S][v] = number of paths using set S, of length |S|, ending at v
    # position = |S| - 1 (0-indexed) when the path reaches v
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    # N[v][j] = number of Ham paths where v is at position j
    N = [[0]*n for _ in range(n)]

    # Position 0: each vertex starts some paths
    # (counted after building full paths)

    # Build DP layer by layer
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[u][v] and (S_prev, u) in dp:
                        total += dp[(S_prev, u)]
                if total > 0:
                    dp[(S, v)] = total

    # Now compute N[v][j] by tracking positions
    # We need a different DP: dp2[S][v] = number of paths of vertex set S ending at v
    # and for each, v is at position |S|-1.
    # To get N[v][j], we need: for each (S, v) with |S| = j+1 and v is the last vertex,
    # count paths through S ending at v that can be EXTENDED to full Hamiltonian paths.

    # Forward DP already done. Need backward DP too.
    # Backward: dp_back[S][v] = number of paths using set S starting at v
    dp_back = {}
    for v in range(n):
        dp_back[(1 << v, v)] = 1

    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[v][u] and (S_prev, u) in dp_back:
                        total += dp_back[(S_prev, u)]
                if total > 0:
                    dp_back[(S, v)] = total

    # N[v][j] = sum over all S containing v with |S|=j+1,
    #           and all T = complement(S) union {v}... no.
    # Actually: N[v][j] = sum_{S: v in S, |S|=j+1} dp[(S,v)] * dp_back[(complement(S)|v, v)]
    # Wait, need to split: left part is S ending at v (position j), right part is complement starting at v.
    # But the right part uses complement(S) ∪ {v}... no, the right part uses vertices NOT in S except v is shared.

    # Better: N[v][j] = Σ over (S_left, S_right) where:
    #   S_left has j+1 vertices including v, last vertex is v
    #   S_right has n-j vertices including v, first vertex is v
    #   S_left ∪ S_right = {0,...,n-1}, S_left ∩ S_right = {v}

    full = (1 << n) - 1
    for v in range(n):
        for j in range(n):
            left_size = j + 1
            right_size = n - j
            # S_left must have left_size elements including v, ending at v
            # S_right must be complement(S_left) | (1<<v), starting at v

            total = 0
            # Enumerate all S_left of size left_size containing v
            others = [u for u in range(n) if u != v]
            for chosen in combinations(others, left_size - 1):
                S_left = (1 << v)
                for u in chosen:
                    S_left |= (1 << u)
                S_right = (full ^ S_left) | (1 << v)

                fwd = dp.get((S_left, v), 0)
                bwd = dp_back.get((S_right, v), 0)
                total += fwd * bwd

            N[v][j] = total

    return N

def is_position_uniform(N, n, H):
    if H % n != 0:
        return False
    target = H // n
    return all(N[v][j] == target for v in range(n) for j in range(n))

# =====================================================================
print("=" * 70)
print("TESTING CONSTRUCTION FAMILIES AT n=7")
print("=" * 70)

n = 7

# Family 1: Circulant tournaments (known to have scalar M)
print("\n  Family 1: Circulant tournaments")
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

for gs in sorted(gen_sets):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in gs:
                A[i][j] = 1
    H = ham_path_dp(A)
    N = position_matrix_dp(A)
    pu = is_position_uniform(N, n, H)
    print(f"    gen={sorted(gs)}: H={H}, PU={pu}")

# Family 2: "Almost circulant" — circulant with one edge flipped
print("\n  Family 2: Circulant + 1 flip (sample)")
count_pu = 0
count_tested = 0
for gs in sorted(gen_sets)[:2]:  # Just first 2 circulants
    A_base = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in gs:
                A_base[i][j] = 1

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            A_flip = [row[:] for row in A_base]
            A_flip[i][j] = 1 - A_flip[i][j]
            A_flip[j][i] = 1 - A_flip[j][i]
            H = ham_path_dp(A_flip)
            count_tested += 1
            if H % n == 0:
                N = position_matrix_dp(A_flip)
                if is_position_uniform(N, n, H):
                    count_pu += 1
                    scores = tuple(sorted(sum(row) for row in A_flip))
                    print(f"    gen={sorted(gs)}, flip ({i},{j}): H={H}, PU=True, scores={scores}")

print(f"  Tested {count_tested}, found {count_pu} PU")

# Family 3: Random tournaments
print("\n  Family 3: Random sampling (1000 tournaments)")
random.seed(42)
count_pu_rand = 0
count_h_div = 0
for trial in range(1000):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    H = ham_path_dp(A)
    if H % n == 0:
        count_h_div += 1
        N = position_matrix_dp(A)
        if is_position_uniform(N, n, H):
            count_pu_rand += 1
            scores = tuple(sorted(sum(row) for row in A))
            # Check if circulant
            is_circ = False
            for gs_check in gen_sets:
                A_check = [[0]*n for _ in range(n)]
                for i in range(n):
                    for j in range(n):
                        if i != j and (j-i)%n in gs_check:
                            A_check[i][j] = 1
                for perm in permutations(range(n)):
                    if all(A[perm[i]][perm[j]] == A_check[i][j] for i in range(n) for j in range(n)):
                        is_circ = True
                        break
                if is_circ:
                    break
            print(f"    Trial {trial}: H={H}, scores={scores}, circulant={is_circ}")

print(f"  H divisible by 7: {count_h_div}/1000")
print(f"  Position-uniform: {count_pu_rand}/1000")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
