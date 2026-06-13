#!/usr/bin/env python3
"""
FAST Betti concentration test — avoids running validation from path_homology_v2.
Tests the conjecture: at most one β_p (p≥1) is nonzero for any tournament.

opus-2026-03-15-S72d
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import random
import time

# ---- Inline the needed functions from path_homology_v2 ----

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(p + 1)]

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0: return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        betti.append(max(0, ker_dim - im_dim))
    return betti

# ---- Tournament generation ----

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

# ---- MAIN ----
print("=" * 70)
print("BETTI CONCENTRATION — FAST VERSION")
print("=" * 70)

# n=5: exhaustive
print("\n--- n=5 EXHAUSTIVE ---")
n = 5; E = n*(n-1)//2
betti_dist = Counter(); multi = []
t0 = time.time()
for bits in range(1 << E):
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1
    nz = [p for p in range(1, len(betti)) if betti[p] > 0]
    if len(nz) > 1: multi.append((bits, betti))
print(f"  Time: {time.time()-t0:.1f}s, Total: {1<<E}")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  Multi-nonzero: {len(multi)}")

# n=6: exhaustive (this is the big one)
print("\n--- n=6 EXHAUSTIVE ---")
n = 6; E = n*(n-1)//2
betti_dist = Counter(); multi = []; nz_count = Counter()
t0 = time.time()
for bits in range(1 << E):
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1
    nz = [p for p in range(1, len(betti)) if betti[p] > 0]
    nz_count[len(nz)] += 1
    if len(nz) > 1: multi.append((bits, betti))
    if (bits+1) % 5000 == 0:
        print(f"  Progress: {bits+1}/{1<<E} ({100*(bits+1)/(1<<E):.1f}%), "
              f"time: {time.time()-t0:.0f}s, multi: {len(multi)}", flush=True)
print(f"  Time: {time.time()-t0:.1f}s, Total: {1<<E}")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  Nonzero count dist: {dict(nz_count)}")
print(f"  Multi-nonzero: {len(multi)}")
if multi:
    for bits, betti in multi[:10]:
        print(f"    VIOLATION: bits={bits}, β={betti}")

# n=7: sample 2000
print("\n--- n=7 SAMPLED (2000) ---")
n = 7; E = n*(n-1)//2
random.seed(42)
betti_dist = Counter(); multi = []; nz_dims = Counter()
t0 = time.time()
for trial in range(2000):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1
    nz = tuple(p for p in range(1, len(betti)) if betti[p] > 0)
    nz_dims[nz] += 1
    if len(nz) > 1: multi.append((bits, betti))
print(f"  Time: {time.time()-t0:.1f}s")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  Nonzero dimensions:")
for dims, cnt in sorted(nz_dims.items()):
    print(f"    p∈{dims}: {cnt}")
print(f"  Multi-nonzero: {len(multi)}")
if multi:
    for bits, betti in multi[:10]:
        print(f"    VIOLATION: bits={bits}, β={betti}")

# n=8: sample 300 (expensive)
print("\n--- n=8 SAMPLED (300) ---")
n = 8; E = n*(n-1)//2
random.seed(99)
betti_dist = Counter(); multi = []; nz_dims = Counter()
t0 = time.time()
for trial in range(300):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    try:
        betti = path_betti_numbers(A, n, max_dim=5)  # up to β₅
        bt = tuple(betti)
        betti_dist[bt] += 1
        nz = tuple(p for p in range(1, len(betti)) if betti[p] > 0)
        nz_dims[nz] += 1
        if len(nz) > 1: multi.append((bits, betti))
    except Exception as e:
        print(f"  Error at trial {trial}: {e}")
print(f"  Time: {time.time()-t0:.1f}s")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  Nonzero dimensions:")
for dims, cnt in sorted(nz_dims.items()):
    print(f"    p∈{dims}: {cnt}")
print(f"  Multi-nonzero: {len(multi)}")
if multi:
    for bits, betti in multi[:10]:
        print(f"    VIOLATION: bits={bits}, β={betti}")

# Euler characteristic check
print("\n" + "=" * 70)
print("EULER CHARACTERISTIC CHECK")
print("=" * 70)
print("If concentration holds, χ = 1 (contractible) or χ = 1 + (-1)^k (sphere S^k)")
print("  Odd k: χ = 0;  Even k: χ = 2")
# Re-scan the n=7 betti dist for Euler char
chi_dist = Counter()
for bt, cnt in betti_dist.items():
    chi = sum((-1)**p * b for p, b in enumerate(bt))
    chi_dist[chi] += cnt
print(f"  n=8: Euler chars: {dict(chi_dist)}")
print("\nDone.")
