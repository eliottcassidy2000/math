#!/usr/bin/env python3
"""β_2 = 0 proof approach via tournament completeness.

KEY INSIGHT from oriented graph analysis:
β_2 > 0 requires "twin vertices" — vertices with identical neighborhoods.
In tournaments, a→b or b→a breaks the twin symmetry, creating DT 3-paths
that fill any potential 2-hole.

PROOF STRATEGY:
For a 2-cycle z ∈ ker(∂_2) ∩ Ω_2 in a tournament T:
1. Every 2-path (a,b,c) in z has a→b→c.
2. For any two vertices u,v in T, one beats the other (say u→v).
3. This edge u→v enables 3-paths (u,v,*,*) or (*,*,u,v) that provide
   boundaries in Ω_3.

QUESTION: Can we prove that for any z ∈ ker(∂_2)∩Ω_2, there exists
w ∈ Ω_3 with ∂_3(w) = z, by constructing w explicitly from the
tournament edges?

TEST: For each 2-cycle in a tournament at n=5,6, decompose the filling
chain w into contributions from specific edges/DT-paths.
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def analyze_filling_chain(A, n):
    """For each 2-cycle, find the filling 3-chain and analyze its structure."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        return None

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        return None

    # ker(∂_2|Ω_2) basis
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    # Full Ω_3
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if dim_om3 == 0:
        return 'no_om3'

    # ∂_3(Ω_3)
    bd3 = build_full_boundary_matrix(a3_list, a2_list)
    im3 = bd3 @ om3

    results = []
    for ci in range(ker_dim):
        z = ker_basis[ci]
        # Find w in Ω_3 with ∂_3(w) = z
        w_om3, res, _, _ = np.linalg.lstsq(im3, z, rcond=None)
        err = np.max(np.abs(im3 @ w_om3 - z))
        if err > 1e-6:
            results.append('unfilled')
            continue

        # w in A_3 coordinates
        w = om3 @ w_om3
        # How many nonzero terms?
        nz = np.sum(np.abs(w) > 1e-8)

        # Classify terms
        a2_set = set(a2_list)
        dt_terms = 0
        cancel_terms = 0
        for j in range(len(a3_list)):
            if abs(w[j]) > 1e-8:
                p = a3_list[j]
                a, b, c, d = p
                faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
                if all(f in a2_set for f in faces):
                    dt_terms += 1
                else:
                    cancel_terms += 1

        results.append({
            'nz': nz, 'dt': dt_terms, 'cancel': cancel_terms,
            'err': err
        })

    return results

# n=5 exhaustive
print("=" * 70)
print("FILLING CHAIN ANALYSIS AT n=5")
print("=" * 70)

fill_stats = Counter()
total_cycles = 0
total_with_cancel = 0

for A in all_tournaments_gen(5):
    info = analyze_filling_chain(A, 5)
    if info is None or info == 'no_om3':
        continue

    for r in info:
        if r == 'unfilled':
            fill_stats['unfilled'] += 1
        else:
            total_cycles += 1
            fill_stats[(r['nz'], r['dt'], r['cancel'])] += 1
            if r['cancel'] > 0:
                total_with_cancel += 1

print(f"Total filling chains: {total_cycles}")
print(f"With cancellation 3-chains: {total_with_cancel}")
print(f"\n(#terms, #DT, #cancel): count")
for k in sorted(fill_stats):
    print(f"  {k}: {fill_stats[k]}")

# n=6 — sample
print(f"\n\n{'='*70}")
print("FILLING CHAIN ANALYSIS AT n=6 (sample)")
print("="*70)

import random
random.seed(42)
fill_stats6 = Counter()
total6 = 0
cancel6 = 0
count6 = 0

for A in all_tournaments_gen(6):
    count6 += 1
    if random.random() > 0.02:  # ~2% sample = ~650 tournaments
        continue
    info = analyze_filling_chain(A, 6)
    if info is None or info == 'no_om3':
        continue
    for r in info:
        if r == 'unfilled':
            fill_stats6['unfilled'] += 1
        else:
            total6 += 1
            key = (r['nz'], r['dt'], r['cancel'])
            fill_stats6[key] += 1
            if r['cancel'] > 0:
                cancel6 += 1

print(f"Total filling chains (sample): {total6}")
print(f"With cancellation 3-chains: {cancel6}")
print(f"\n(#terms, #DT, #cancel): count")
for k in sorted(fill_stats6)[:20]:
    print(f"  {k}: {fill_stats6[k]}")

# KEY QUESTION: At n=6, is there a UNIQUE minimal filling chain,
# or are there many choices? If unique, the filling is canonical.
print(f"\n\n{'='*70}")
print("FILLING CHAIN UNIQUENESS AT n=5")
print("="*70)

unique_count = 0
multi_count = 0
for A in all_tournaments_gen(5):
    a1 = enumerate_allowed_paths(A, 5, 1)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a3 = enumerate_allowed_paths(A, 5, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, 5, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0: continue

    om3 = compute_omega_basis(A, 5, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if dim_om3 == 0: continue

    bd3 = build_full_boundary_matrix(a3_list, a2_list)
    im3 = bd3 @ om3

    # rank of ∂_3|Ω_3
    rank3 = np.linalg.matrix_rank(im3, tol=1e-8)
    # ker(∂_3|Ω_3)
    ker3 = dim_om3 - rank3

    if ker3 == 0:
        unique_count += 1
    else:
        multi_count += 1

print(f"Unique filling (ker ∂_3|Ω_3 = 0): {unique_count} tournaments")
print(f"Multiple fillings (ker ∂_3|Ω_3 > 0): {multi_count} tournaments")

print("\nDone.")
