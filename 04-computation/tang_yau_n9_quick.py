#!/usr/bin/env python3
"""
tang_yau_n9_quick.py - Quick enumeration of circulant tournaments at n=9

Lists all valid S (|S|=4, S cap (-S) = empty) and computes beta_5 only
(focused computation, not full Betti vector).

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            if j != i:
                A[i][j] = 1
    return A

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def compute_beta5_focused(A, n):
    """Compute beta_5 using only Om_4, Om_5, Om_6."""
    allowed = {}
    for p in range(-1, 7):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
            if not allowed[p]:
                return 0

    omega_basis = {}
    for p in [4, 5, 6]:
        basis = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_basis[p] = basis

    dim5 = omega_basis[5].shape[1] if omega_basis[5].ndim == 2 else 0
    if dim5 == 0:
        return 0

    bd5 = build_full_boundary_matrix(allowed[5], allowed[4]) @ omega_basis[5]
    rk5 = np.linalg.matrix_rank(bd5, tol=1e-8)

    dim6 = omega_basis[6].shape[1] if omega_basis[6].ndim == 2 else 0
    if dim6 == 0:
        rk6 = 0
    else:
        bd6 = build_full_boundary_matrix(allowed[6], allowed[5]) @ omega_basis[6]
        rk6 = np.linalg.matrix_rank(bd6, tol=1e-8)

    return (dim5 - rk5) - rk6

print("=" * 70)
print("CIRCULANT TOURNAMENTS AT n=9")
print("=" * 70)

n = 9
valid_sets = []
for S_tuple in combinations(range(1, n), 4):
    S = frozenset(S_tuple)
    nS = frozenset((n - s) % n for s in S)
    if S & nS:
        continue
    # Canonical: pick lexicographically smaller of S and -S
    S_sorted = tuple(sorted(S))
    nS_sorted = tuple(sorted(nS))
    if S_sorted > nS_sorted:
        continue  # skip duplicate (we'll test -S as S)
    valid_sets.append(S_sorted)

print(f"Total valid circulant tournament connection sets: {len(valid_sets)}")
print()

for S_tuple in valid_sets:
    S = set(S_tuple)
    A = build_circulant(n, S)
    H = H_tournament(A, n)
    t0 = time.time()
    b5 = compute_beta5_focused(A, n)
    dt = time.time() - t0
    marker = f" *** beta5={b5}" if b5 > 0 else ""
    print(f"  S={sorted(S)}: H={H}, beta_5={b5} ({dt:.1f}s){marker}")

print("\nDone.")
