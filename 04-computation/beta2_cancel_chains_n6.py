#!/usr/bin/env python3
"""Analyze the cancellation 3-chains needed for β_2=0 at n=6.

Focus on the 960 tournaments where DT-only fails.
What do the cancellation 3-chains look like?
Are they formed by pairs of 3-paths sharing a non-A_2 face?
"""
import numpy as np
from itertools import combinations
import sys, time
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

n = 6
print("=" * 70)
print("CANCELLATION 3-CHAINS AT n=6")
print("=" * 70)

# Find one DT-failure tournament and study it in detail
t0 = time.time()
found = 0
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]
    a2_set = set(a2_list)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2

    # DT paths
    dt_paths = []
    for p in a3_list:
        a, b, c, d = tuple(p)
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        if all(f in a2_set for f in faces):
            dt_paths.append(tuple(p))

    a2_idx = {p: i for i, p in enumerate(a2_list)}
    dt_bds = []
    for p in dt_paths:
        bd_vec = np.zeros(len(a2_list))
        a, b, c, d = p
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        for fi, f in enumerate(faces):
            if f in a2_idx:
                bd_vec[a2_idx[f]] += (-1)**fi
        dt_bds.append(bd_vec)

    if not dt_bds:
        continue

    dt_bd_matrix = np.column_stack(dt_bds)
    dt_in_om2, _, _, _ = np.linalg.lstsq(om2, dt_bd_matrix, rcond=None)
    rank_dt = np.linalg.matrix_rank(dt_in_om2, tol=1e-8)

    if rank_dt >= ker_dim:
        continue  # DT sufficient, skip

    found += 1
    if found > 2:
        continue  # Just show first 2

    print(f"\n--- Tournament (DT gap = {ker_dim - rank_dt}) ---")
    for i in range(n):
        out = [j for j in range(n) if A[i][j]]
        print(f"  {i} → {out}")

    out_deg = [sum(A[i]) for i in range(n)]
    t3 = sum(1 for a, b, c in combinations(range(n), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    print(f"  out_deg={out_deg}, t3={t3}")
    print(f"  dim(Ω_2)={dim_om2}, ker(∂_2|Ω_2)={ker_dim}, rank(∂_3|DT)={rank_dt}")

    # Full Ω_3
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    print(f"  dim(Ω_3)={dim_om3}, |DT|={len(dt_paths)}")

    # Find cancellation 3-chains (Ω_3 basis elements that aren't individual DT paths)
    # An Ω_3 basis element is "cancellation" if it has nonzero coefficient on a non-DT path
    dt_set = set(dt_paths)
    a3_idx = {p: i for i, p in enumerate(a3_list)}

    cancel_elements = []
    dt_elements = []
    for j in range(dim_om3):
        col = om3[:, j]
        nz_paths = [(a3_list[k], col[k]) for k in range(len(a3_list)) if abs(col[k]) > 1e-8]
        has_non_dt = any(p not in dt_set for p, c in nz_paths)
        if has_non_dt:
            cancel_elements.append((j, nz_paths))
        else:
            dt_elements.append((j, nz_paths))

    print(f"\n  Ω_3 basis: {len(dt_elements)} DT elements, {len(cancel_elements)} cancellation elements")

    # Show first few cancellation elements
    for ci, (j, nz_paths) in enumerate(cancel_elements[:3]):
        print(f"\n  Cancellation element #{ci+1} ({len(nz_paths)} terms):")
        for p, c in nz_paths:
            is_dt = p in dt_set
            tag = "DT" if is_dt else "non-DT"
            # For non-DT: which face is not in A_2?
            a, b, c_v, d = p
            faces = [(b,c_v,d), (a,c_v,d), (a,b,d), (a,b,c_v)]
            bad_faces = [f for f in faces if f not in a2_set]
            bad_str = f" [bad: {bad_faces}]" if bad_faces else ""
            print(f"    {c:+.4f} * {p} [{tag}]{bad_str}")

    # KEY QUESTION: Do the non-DT paths in cancellation chains come in PAIRS
    # with matching bad faces that cancel?
    print(f"\n  Non-DT path pairs with matching bad faces:")
    for ci, (j, nz_paths) in enumerate(cancel_elements[:3]):
        bad_face_groups = defaultdict(list)
        for p, coeff in nz_paths:
            if p not in dt_set:
                a, b, c_v, d = p
                faces = [(b,c_v,d), (a,c_v,d), (a,b,d), (a,b,c_v)]
                bad = tuple(f for f in faces if f not in a2_set)
                bad_face_groups[bad].append((p, coeff))

        for bad, paths in bad_face_groups.items():
            coeffs = [c for _, c in paths]
            total = sum(coeffs)
            print(f"    Bad face {bad}: {len(paths)} paths, coeff sum = {total:.4f}")

print(f"\n\nTotal DT-failure tournaments found: {found} (in {time.time()-t0:.0f}s)")

# ===== Part 2: Minimal cancellation chain structure =====
print(f"\n\n{'='*70}")
print("CANCELLATION CHAIN MECHANISM")
print("="*70)
print("For each pair of non-DT 3-paths sharing a bad face,")
print("check if their difference is in Ω_3.")

# Take the first failure tournament
found = 0
for A in all_tournaments_gen(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_set = set(tuple(p) for p in a2)
    a3_list = [tuple(p) for p in a3]
    a3_idx = {tuple(p): i for i, p in enumerate(a3_list)}

    dt_paths = []
    non_dt = []
    for p in a3_list:
        a, b, c, d = tuple(p)
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        if all(f in a2_set for f in faces):
            dt_paths.append(tuple(p))
        else:
            non_dt.append(tuple(p))

    # Check DT sufficiency
    a1 = enumerate_allowed_paths(A, n, 1)
    a2_list = [tuple(p) for p in a2]
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue
    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2

    a2_idx = {p: i for i, p in enumerate(a2_list)}
    dt_bds = []
    for p in dt_paths:
        bd_vec = np.zeros(len(a2_list))
        a, b, c, d = p
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        for fi, f in enumerate(faces):
            if f in a2_idx:
                bd_vec[a2_idx[f]] += (-1)**fi
        dt_bds.append(bd_vec)

    if not dt_bds: continue
    dt_bd_matrix = np.column_stack(dt_bds)
    dt_in_om2, _, _, _ = np.linalg.lstsq(om2, dt_bd_matrix, rcond=None)
    rank_dt = np.linalg.matrix_rank(dt_in_om2, tol=1e-8)
    if rank_dt >= ker_dim: continue

    found += 1
    if found > 1: break

    # Group non-DT paths by their bad face
    bad_face_groups = defaultdict(list)
    for p in non_dt:
        a, b, c, d = p
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        bad = tuple(f for f in faces if f not in a2_set)
        bad_face_groups[bad].append(p)

    print(f"\nNon-DT 3-paths grouped by bad face:")
    for bad, paths in sorted(bad_face_groups.items(), key=lambda x: len(x[1]), reverse=True):
        print(f"  Bad face {bad}: {len(paths)} paths")
        for p in paths:
            print(f"    {p}")

    # For each pair in a group, form the difference and check if it's in Ω_3
    print(f"\n  Pairwise differences in Ω_3?")
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    for bad, paths in bad_face_groups.items():
        if len(paths) < 2:
            continue
        for i in range(len(paths)):
            for j in range(i+1, len(paths)):
                diff = np.zeros(len(a3_list))
                diff[a3_idx[paths[i]]] = 1
                diff[a3_idx[paths[j]]] = -1

                # Check if diff is in column span of om3
                x, _, _, _ = np.linalg.lstsq(om3, diff, rcond=None)
                err = np.max(np.abs(om3 @ x - diff))
                in_om3 = err < 1e-8

                if in_om3:
                    print(f"    {paths[i]} - {paths[j]} ∈ Ω_3 ✓ (bad: {bad})")
                else:
                    print(f"    {paths[i]} - {paths[j]} ∉ Ω_3 (err={err:.2e})")

print("\nDone.")
