#!/usr/bin/env python3
"""
beta2_multicone_n7.py — Multi-vertex cone at n=7,8 (sampled)

Tests whether 2 cones always suffice.

Author: opus-2026-03-08-S49
"""
import sys, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_cone_image(A, n, v, ap2, ap3_list):
    C = np.zeros((len(ap3_list), len(ap2)))
    for j, path in enumerate(ap2):
        if v in path:
            continue
        if not A[v][path[0]]:
            continue
        cone_path = tuple([v] + list(path))
        if cone_path in ap3_list:
            C[ap3_list.index(cone_path), j] = 1.0
    return C


print("=" * 70)
print("MULTI-VERTEX CONE AT n=7,8 (SAMPLED)")
print("=" * 70)

for n, samples in [(7, 500), (8, 100)]:
    t0 = time.time()
    min_vertices_needed = Counter()
    max_needed = 0

    for trial in range(samples):
        A = build_random_adj(n)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0:
            min_vertices_needed[0] += 1
            continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z2_dim = d2 - rk

        if z2_dim == 0:
            min_vertices_needed[0] += 1
            continue

        z2_om = Vt[rk:, :]
        z2_A2 = om2 @ z2_om.T

        ap3_list = [tuple(p) for p in ap3]
        bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

        z2_pinv = np.linalg.pinv(z2_A2)

        # Compute projected cone fills for all vertices
        projected = {}
        for v in range(n):
            C_v = compute_cone_image(A, n, v, ap2, ap3_list)
            fill_v = bd3 @ C_v @ z2_A2
            projected[v] = z2_pinv @ fill_v

        # Try increasing subset sizes
        from itertools import combinations
        found = False
        for k in range(1, n+1):
            if found:
                break
            for S_set in combinations(range(n), k):
                cols = np.hstack([projected[v] for v in S_set])
                if np.linalg.matrix_rank(cols, tol=1e-8) >= z2_dim:
                    min_vertices_needed[k] += 1
                    max_needed = max(max_needed, k)
                    found = True
                    break

        if not found:
            min_vertices_needed[-1] += 1
            max_needed = max(max_needed, n+1)

        if (trial+1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{samples} ({elapsed:.0f}s) max_needed={max_needed}")

    elapsed = time.time() - t0
    print(f"\nn={n}: {samples} tournaments in {elapsed:.0f}s")
    print(f"  Min vertices needed: {dict(sorted(min_vertices_needed.items()))}")
    print(f"  Max needed: {max_needed}")
    if -1 in min_vertices_needed:
        print(f"  FAILURES: {min_vertices_needed[-1]}")

print("\nDone.")
