#!/usr/bin/env python3
"""Compute Ω dims for the TRANSITIVE tournament at various n.
For transitive T_n: 0→1→2→...→n-1 with i→j for all i<j.

Expected: Ω_k = C(n, k+1) (binomial coefficient).
Verify this and check the chain complex structure.
"""
import numpy as np
from itertools import combinations
import sys
from math import comb
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
)

print("=" * 70)
print("TRANSITIVE TOURNAMENT: Ω DIMENSIONS")
print("=" * 70)

for n in range(3, 10):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1  # i→j for all i<j

    dims = [n]
    for d in range(1, n):
        ap = enumerate_allowed_paths(A, n, d)
        if d == 1:
            dims.append(len(ap))
        else:
            ap_dm1 = enumerate_allowed_paths(A, n, d-1)
            om = compute_omega_basis(A, n, d, ap, ap_dm1)
            dims.append(om.shape[1] if om.ndim == 2 else 0)

    binomial = [comb(n, k+1) for k in range(n)]
    match = dims == binomial
    chi = sum((-1)**k * dims[k] for k in range(n))

    print(f"  n={n}: Ω = {dims}")
    print(f"       C(n,k+1) = {binomial}")
    print(f"       Match: {'✓' if match else '✗'}, χ = {chi}")

    if n <= 8:
        betti = path_betti_numbers(A, n, max_dim=n-1)
        print(f"       β = {betti}")

# ===== Also check the REGULAR tournament (all out-deg = (n-1)/2) =====
print(f"\n\n{'='*70}")
print("REGULAR TOURNAMENT (ROTATIONAL): Ω DIMENSIONS")
print("="*70)

for n in [5, 7]:
    # Regular tournament: circulant with S = {1, 2, ..., (n-1)/2}
    k = (n-1)//2
    S = list(range(1, k+1))
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1

    dims = [n]
    for d in range(1, n):
        ap = enumerate_allowed_paths(A, n, d)
        if d == 1:
            dims.append(len(ap))
        else:
            ap_dm1 = enumerate_allowed_paths(A, n, d-1)
            om = compute_omega_basis(A, n, d, ap, ap_dm1)
            dims.append(om.shape[1] if om.ndim == 2 else 0)

    chi = sum((-1)**k * dims[k] for k in range(n))
    betti = path_betti_numbers(A, n, max_dim=n-1)

    # Check palindromicity
    is_palin = all(dims[k] == dims[n-1-k] for k in range(n))

    print(f"  n={n}, S={S}: Ω = {dims}")
    print(f"    Palindromic: {'✓' if is_palin else '✗'}, χ = {chi}")
    print(f"    β = {betti}")

    t3 = sum(1 for a, b, c in combinations(range(n), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    print(f"    t3 = {t3} (max = {comb(n,3)})")

print("\nDone.")
