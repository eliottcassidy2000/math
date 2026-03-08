#!/usr/bin/env python3
"""Verify subtournament boundary decomposition at n=6.
Does im(∂_3|subtournament DT) = ker(∂_2|Ω_2) for ALL 32768 tournaments?"""
import numpy as np
from itertools import combinations
import sys, time
from collections import Counter
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
print(f"SUBTOURNAMENT BOUNDARY DECOMPOSITION AT n={n}")
print("=" * 70)

t0 = time.time()
success = 0
failure = 0
total_with_ker = 0
total_cycles = 0
count = 0
failure_details = []

for A in all_tournaments_gen(n):
    count += 1
    if count % 5000 == 0:
        print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0: continue

    total_with_ker += 1
    total_cycles += ker_dim

    # Collect boundaries of ALL DT paths from 4-vertex subtournaments
    a2_idx = {p: i for i, p in enumerate(a2_list)}
    a2_set = set(a2_list)

    all_sub_bds = []
    for p in a3_list:
        a, b, c, d = tuple(p)
        # Check if this is a DT path (all faces in A_2)
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        if all(f in a2_set for f in faces):
            bd_vec = np.zeros(len(a2_list))
            for fi, f in enumerate(faces):
                bd_vec[a2_idx[f]] += (-1)**fi
            all_sub_bds.append(bd_vec)

    if not all_sub_bds:
        failure += ker_dim
        continue

    sub_bd_matrix = np.column_stack(all_sub_bds)

    # Project to Ω_2 coordinates
    sub_in_om2, _, _, _ = np.linalg.lstsq(om2, sub_bd_matrix, rcond=None)
    err = np.max(np.abs(om2 @ sub_in_om2 - sub_bd_matrix))
    if err > 1e-6:
        pass  # DT boundaries not in Ω_2? This shouldn't happen

    # ker basis
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)

    # Check each cycle
    for ci in range(ker_dim):
        z_om2 = Vt[rank2 + ci, :]
        combined = np.column_stack([sub_in_om2, z_om2.reshape(-1,1)])
        rank_without = np.linalg.matrix_rank(sub_in_om2, tol=1e-8)
        rank_with = np.linalg.matrix_rank(combined, tol=1e-8)

        if rank_with == rank_without:
            success += 1
        else:
            failure += 1
            if len(failure_details) < 5:
                t3 = sum(1 for a, b, c in combinations(range(n), 3)
                         if (A[a][b] and A[b][c] and A[c][a]) or
                            (A[b][a] and A[a][c] and A[c][b]))
                out_deg = sorted([sum(A[i]) for i in range(n)])
                failure_details.append({
                    't3': t3,
                    'out_deg': out_deg,
                    'ker_dim': ker_dim,
                    'rank_DT': rank_without,
                })

t1 = time.time()
print(f"\nDone in {t1-t0:.1f}s")
print(f"\nTotal tournaments with ker > 0: {total_with_ker}")
print(f"Total 2-cycles: {total_cycles}")
print(f"Filled by DT boundaries: {success}")
print(f"NOT filled (gap): {failure}")

if failure_details:
    print(f"\nFailure details:")
    for fd in failure_details:
        print(f"  t3={fd['t3']}, out_deg={fd['out_deg']}, ker={fd['ker_dim']}, rank_DT={fd['rank_DT']}")

# Also test: do ALL Ω_3 boundaries (including cancellation chains) fill everything?
print(f"\n\n{'='*70}")
print("FULL Ω_3 BOUNDARY TEST (should be 0 failures if β_2=0)")
print("="*70)

full_success = 0
full_failure = 0
count = 0
for A in all_tournaments_gen(n):
    count += 1
    if count % 5000 == 0:
        print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0: continue

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd3 = build_full_boundary_matrix(a3_list, a2_list)
    if dim_om3 > 0:
        im3 = bd3 @ om3
        im3_in_om2, _, _, _ = np.linalg.lstsq(om2, im3, rcond=None)
        rank3 = np.linalg.matrix_rank(im3_in_om2, tol=1e-8)
    else:
        rank3 = 0

    if rank3 >= ker_dim:
        full_success += 1
    else:
        full_failure += 1

print(f"Full Ω_3 success: {full_success}, failure: {full_failure}")

print("\nDone.")
