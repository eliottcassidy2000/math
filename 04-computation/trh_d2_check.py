#!/usr/bin/env python3
"""
trh_d2_check.py — opus-2026-03-13-S71

Check whether d²=0 holds for TRH (interior boundary on regular paths)
for circulant vs non-circulant tournaments.
"""

import numpy as np
from itertools import combinations

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def interior_boundary(paths_m, paths_m1):
    if not paths_m or not paths_m1:
        return np.zeros((len(paths_m1) if paths_m1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    idx = {p: i for i, p in enumerate(paths_m1)}
    m = len(paths_m[0]) - 1
    B = np.zeros((len(paths_m1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i
    return B

def check_d2(A, label):
    n = A.shape[0]
    all_paths = {m: get_regular_paths(A, m) for m in range(n)}
    ok = True
    for m in range(2, n):
        if all_paths[m] and all_paths[m-1] and all_paths[m-2]:
            B_m = interior_boundary(all_paths[m], all_paths[m-1])
            B_m1 = interior_boundary(all_paths[m-1], all_paths[m-2])
            d2 = B_m1 @ B_m
            if np.max(np.abs(d2)) > 0:
                ok = False
    return ok

# ============================================================
print("="*70)
print("TRH d²=0 CHECK FOR CIRCULANT TOURNAMENTS")
print("="*70)

# All circulants at n=5
print("\nAll circulant tournaments at n=5:")
for S_tuple in combinations(range(1,5), 2):
    S = set(S_tuple)
    comp = {5-s for s in S}
    if S & comp == set() and S | comp == set(range(1,5)):
        A = circulant_tournament(5, S)
        ok = check_d2(A, f'C5_{list(S_tuple)}')
        print(f"  C5_{list(S_tuple)}: d²=0? {ok}")

# All circulants at n=7
print("\nAll circulant tournaments at n=7:")
for S_tuple in combinations(range(1,7), 3):
    S = set(S_tuple)
    comp = {7-s for s in S}
    if S & comp == set() and S | comp == set(range(1,7)):
        A = circulant_tournament(7, S)
        ok = check_d2(A, f'C7_{list(S_tuple)}')
        print(f"  C7_{list(S_tuple)}: d²=0? {ok}")

# All circulants at n=9
print("\nSelected circulant tournaments at n=9:")
for S_tuple in [[1,2,3,4], [1,3,5,7]]:
    S = set(S_tuple)
    A = circulant_tournament(9, S)
    ok = check_d2(A, f'C9_{S_tuple}')
    print(f"  C9_{S_tuple}: d²=0? {ok}")

# Paley
print("\nPaley tournaments:")
for p in [3, 7, 11]:
    if p % 4 != 3: continue
    QR = {a % p for a in range(1, p) if pow(a, (p-1)//2, p) == 1}
    A = circulant_tournament(p, QR)
    ok = check_d2(A, f'Paley p={p}')
    print(f"  Paley p={p}: d²=0? {ok}")

# ============================================================
print(f"\n{'='*70}")
print("TRH d²=0 CHECK FOR ALL n=5 TOURNAMENTS")
print("="*70)

from itertools import permutations as perms

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
d2_pass = 0
d2_fail = 0
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    if check_d2(A, ""):
        d2_pass += 1
    else:
        d2_fail += 1

print(f"  n=5: {d2_pass} pass, {d2_fail} fail out of {2**len(pairs)} tournaments")

# ============================================================
print(f"\n{'='*70}")
print("TRH d²=0 CHECK FOR ALL n=6 TOURNAMENTS")
print("="*70)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
d2_pass = 0
d2_fail = 0
fail_scores = []
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    if check_d2(A, ""):
        d2_pass += 1
    else:
        d2_fail += 1
        if len(fail_scores) < 5:
            scores = tuple(sorted(int(sum(A[v][w] for w in range(n) if w != v))
                                  for v in range(n)))
            fail_scores.append(scores)

print(f"  n=6: {d2_pass} pass, {d2_fail} fail out of {2**len(pairs)} tournaments")
print(f"  Example failing score sequences: {fail_scores}")

print("\nDONE.")
