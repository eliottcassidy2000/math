#!/usr/bin/env python3
"""
Compute Ω dimensions for the three regular n=7 tournament classes.

Track |A_d|, dim(Ω_d), rank(∂_d), and β_d for each.

opus-2026-03-13-S71b
"""

import sys, time
sys.path.insert(0, '04-computation')
from path_homology_v2 import (enumerate_allowed_paths, compute_omega_basis,
                                build_full_boundary_matrix)
import numpy as np
import itertools

def count_hp(adj, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and adj[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def is_regular(adj, n):
    return all(sum(adj[i]) == (n-1)//2 for i in range(n))

def all_tournaments_n7():
    n = 7
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

n = 7
max_d = n - 1

# Find reps
reps = {}
for adj in all_tournaments_n7():
    if not is_regular(adj, n):
        continue
    H = count_hp(adj, n)
    if H not in reps:
        reps[H] = adj
    if len(reps) == 3:
        break

for H in sorted(reps.keys()):
    adj = reps[H]
    print(f"\n{'='*60}")
    print(f"H = {H}")
    print(f"{'='*60}")

    # Compute all allowed paths and omega bases
    allowed = {}
    omega_basis = {}
    for d in range(-1, max_d + 2):
        if d < 0:
            allowed[d] = []
        else:
            allowed[d] = enumerate_allowed_paths(adj, n, d)

    for d in range(max_d + 2):
        omega_basis[d] = compute_omega_basis(adj, n, d, allowed[d], allowed[d-1])

    # Compute boundary maps on Omega and ranks
    A_dims = []
    Omega_dims = []
    ranks = [0]  # R_0 = 0 (no boundary from degree 0)
    betti = []

    for d in range(max_d + 1):
        a_d = len(allowed[d])
        o_d = omega_basis[d].shape[1] if omega_basis[d].ndim == 2 else 0
        A_dims.append(a_d)
        Omega_dims.append(o_d)

    # Compute ranks of ∂_d: Ω_d → Ω_{d-1}
    for d in range(1, max_d + 1):
        if Omega_dims[d] == 0 or Omega_dims[d-1] == 0:
            ranks.append(0)
            continue
        # Boundary in full space
        bd = build_full_boundary_matrix(allowed[d], allowed[d-1])
        # Restrict to Omega
        bd_omega = bd @ omega_basis[d]  # image in A_{d-1} coordinates
        # Project onto Omega_{d-1} coordinates
        omega_prev = omega_basis[d-1]
        # Change to Omega_{d-1} basis: bd_restricted = omega_prev^+ @ bd_omega
        bd_restricted = np.linalg.lstsq(omega_prev, bd_omega, rcond=None)[0]
        rk = np.linalg.matrix_rank(bd_restricted, tol=1e-8)
        ranks.append(rk)

    # Betti = Ω_d - R_d - R_{d+1}
    ranks.append(0)  # R_{max_d+1} = 0
    for d in range(max_d + 1):
        b = Omega_dims[d] - ranks[d] - ranks[d+1]
        betti.append(b)

    print(f"  d    |A_d|   Ω_d    R_d    β_d")
    for d in range(max_d + 1):
        print(f"  {d}    {A_dims[d]:5d}   {Omega_dims[d]:4d}   {ranks[d]:4d}   {betti[d]:4d}")
    print(f"  χ = {sum((-1)**d * betti[d] for d in range(max_d+1))}")
    print(f"  β = {tuple(betti)}")
