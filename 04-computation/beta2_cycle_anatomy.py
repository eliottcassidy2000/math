#!/usr/bin/env python3
"""
β_2 = 0 PROOF EXPLORATION: Anatomy of 2-cycles

For each tournament with nontrivial ker(∂_2|Ω_2), we:
1. Find explicit 2-cycles (elements of ker ∂_2 ∩ Ω_2)
2. Decompose them into TT and cancellation-chain components
3. Find which Ω_3 elements map to them under ∂_3
4. Analyze the local 4-vertex structure that kills each cycle

Goal: find the ALGEBRAIC MECHANISM that forces β_2 = 0.
"""
import numpy as np
from itertools import combinations
import sys, time
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

def analyze_cycles(A, n, verbose=True):
    """Find and decompose 2-cycles."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    a1_set = set(a1_list)

    # Ω_2 basis
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        return None

    # TT classification
    tt_idx = [i for i, p in enumerate(a2_list) if A[p[0]][p[2]] == 1]
    ntt_idx = [i for i, p in enumerate(a2_list) if A[p[0]][p[2]] == 0]

    # ∂_2: A_2 → A_1
    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    # Restrict to Ω_2
    bd2_om = bd2 @ om2

    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        return None

    # Find kernel basis of ∂_2|Ω_2
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    # Kernel = last (dim_om2 - rank2) rows of Vt
    ker_basis_in_om2 = Vt[rank2:, :]  # ker_dim × dim_om2
    # Convert to A_2 coordinates
    ker_basis_in_a2 = (om2 @ ker_basis_in_om2.T).T  # ker_dim × len(a2)

    # Ω_3 basis
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # ∂_3: A_3 → A_2
    bd3 = build_full_boundary_matrix(a3_list, a2_list)
    # ∂_3(Ω_3) in A_2 coords
    if dim_om3 > 0:
        im3_in_a2 = bd3 @ om3  # len(a2) × dim_om3
    else:
        im3_in_a2 = np.zeros((len(a2), 0))

    if verbose:
        print(f"\n  dim(Ω_2) = {dim_om2}, |TT| = {len(tt_idx)}, |non-TT| = {len(ntt_idx)}")
        print(f"  ker(∂_2|Ω_2) = {ker_dim}, dim(Ω_3) = {dim_om3}")

        for i in range(ker_dim):
            cyc = ker_basis_in_a2[i]
            # Normalize for readability
            nz = np.abs(cyc) > 1e-8
            if np.sum(nz) == 0:
                continue

            print(f"\n  2-cycle #{i+1}:")
            # TT part
            tt_part = [(a2_list[j], cyc[j]) for j in tt_idx if abs(cyc[j]) > 1e-8]
            ntt_part = [(a2_list[j], cyc[j]) for j in ntt_idx if abs(cyc[j]) > 1e-8]

            if tt_part:
                print(f"    TT components ({len(tt_part)}):")
                for path, coeff in tt_part:
                    print(f"      {coeff:+.4f} * {path}")
            if ntt_part:
                print(f"    non-TT components ({len(ntt_part)}):")
                for path, coeff in ntt_part:
                    # What's the non-allowed face?
                    a, b, c = path
                    face_str = f"({a},{c})" if A[c][a] else f"?"
                    print(f"      {coeff:+.4f} * {path} [non-allowed face: {face_str}]")

            # Which vertices are involved?
            verts = set()
            for j in range(len(a2_list)):
                if abs(cyc[j]) > 1e-8:
                    verts.update(a2_list[j])
            print(f"    Vertices involved: {sorted(verts)}")

            # Check: is this cycle in im(∂_3|Ω_3)?
            if dim_om3 > 0:
                # Solve im3_in_a2 @ x = cyc
                x, res, _, _ = np.linalg.lstsq(im3_in_a2, cyc, rcond=None)
                err = np.max(np.abs(im3_in_a2 @ x - cyc))
                print(f"    In im(∂_3|Ω_3)? err = {err:.2e}")
                if err < 1e-6:
                    # Show which Ω_3 elements contribute
                    om3_coeffs = x
                    for j in range(dim_om3):
                        if abs(om3_coeffs[j]) > 1e-6:
                            # What is this Ω_3 basis element?
                            elem = om3[:, j]
                            nz_3 = [(a3_list[k], elem[k]) for k in range(len(a3_list)) if abs(elem[k]) > 1e-8]
                            if len(nz_3) <= 6:
                                print(f"      Ω_3 element #{j}: {nz_3}")

    return ker_dim, dim_om2, dim_om3

# ===== Find tournaments with nontrivial ker at n=5 =====
print("=" * 70)
print("ANATOMY OF 2-CYCLES IN TOURNAMENTS")
print("=" * 70)

n = 5
count = 0
found = 0
examples = []
for A in all_tournaments_gen(n):
    count += 1
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(
        [tuple(p) for p in a2], [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim > 0 and found < 3:
        print(f"\n--- Tournament #{count} (ker_dim={ker_dim}) ---")
        # Show adjacency
        for i in range(n):
            nbrs = [j for j in range(n) if A[i][j] == 1]
            print(f"  {i} → {nbrs}")
        analyze_cycles(A, n)
        found += 1

# ===== Statistics on cycle composition =====
print(f"\n\n{'='*70}")
print("STATISTICS: 2-CYCLE STRUCTURE AT n=5")
print("="*70)

pure_tt_cycles = 0
mixed_cycles = 0
total_cycles = 0

for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    total_cycles += ker_dim

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    for i in range(ker_dim):
        cyc = ker_basis[i]
        has_ntt = any(abs(cyc[j]) > 1e-8 for j, p in enumerate(a2_list) if A[p[0]][p[2]] == 0)
        if has_ntt:
            mixed_cycles += 1
        else:
            pure_tt_cycles += 1

print(f"Total 2-cycles (kernel vectors): {total_cycles}")
print(f"  Pure TT: {pure_tt_cycles}")
print(f"  Mixed (TT + non-TT): {mixed_cycles}")

# ===== For each cycle, count vertices involved =====
print(f"\n--- Vertex count distribution ---")
from collections import Counter
vert_dist = Counter()
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    for i in range(ker_dim):
        cyc = ker_basis[i]
        verts = set()
        for j in range(len(a2_list)):
            if abs(cyc[j]) > 1e-8:
                verts.update(a2_list[j])
        vert_dist[len(verts)] += 1

print(f"  #vertices in cycle: count")
for k in sorted(vert_dist.keys()):
    print(f"    {k}: {vert_dist[k]}")

# ===== Key question: is every 2-cycle supported on a 4-vertex subset? =====
print(f"\n--- Can cycles be localized to 4 vertices? ---")
n4_sufficient = 0
n4_insufficient = 0
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    for i in range(ker_dim):
        cyc = ker_basis[i]
        verts = set()
        for j in range(len(a2_list)):
            if abs(cyc[j]) > 1e-8:
                verts.update(a2_list[j])
        if len(verts) <= 4:
            n4_sufficient += 1
        else:
            n4_insufficient += 1

print(f"  Supported on ≤4 vertices: {n4_sufficient}")
print(f"  Need 5 vertices: {n4_insufficient}")

print("\nDone.")
