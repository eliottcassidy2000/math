#!/usr/bin/env python3
"""
beta2_dt_cancel_filling.py — Test if DT + cancellation pairs fill Z₂ at n=6

At n=5: DT + single cancellation pairs fill ALL of Z₂ (1024/1024).
Does this extend to n=6?

If yes, this gives a proof path:
1. DT 4-paths are in Ω₃ (all faces in A₂, by tournament completeness)
2. Cancellation pairs (sharing bad face) are in Ω₃ (bad face cancels)
3. Their boundaries span Z₂
4. Therefore im(∂₃) ⊇ Z₂ ⊇ im(∂₃), so im(∂₃) = Z₂, so β₂ = 0

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ======================================================================
print("="*70)
print("DT + CANCELLATION FILLING OF Z₂")
print("="*70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    success = 0
    failure = 0
    total = 0

    for A in all_tournaments(n):
        total += 1
        if total % 5000 == 0:
            print(f"  ... {total}", flush=True)

        a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
        a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
        a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        d_om2 = om2.shape[1] if om2.ndim == 2 else 0
        if d_om2 == 0:
            success += 1
            continue

        # Compute Z₂ = ker(∂₂|Ω₂)
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2
        S = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = sum(s > 1e-8 for s in S)
        z2 = d_om2 - rk2

        if z2 == 0:
            success += 1
            continue

        # Build DT + cancellation space
        a3_idx = {p: i for i, p in enumerate(a3)}

        # DT paths
        dt_indices = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

        # Cancellation pairs: paths sharing a bad face
        bad_face_groups = defaultdict(list)
        for idx, p in enumerate(a3):
            a_, b_, c_, d_ = p
            if not A[a_][c_]:
                bad_face_groups[('acd', a_, c_, d_)].append(idx)
            if not A[b_][d_]:
                bad_face_groups[('abd', a_, b_, d_)].append(idx)

        # Build vectors for DT + cancellation
        vectors = []
        for idx in dt_indices:
            v = np.zeros(len(a3))
            v[idx] = 1
            vectors.append(v)

        for bf, indices in bad_face_groups.items():
            if len(indices) >= 2:
                for j in range(1, len(indices)):
                    v = np.zeros(len(a3))
                    v[indices[0]] = 1
                    v[indices[j]] = -1
                    vectors.append(v)

        if not vectors:
            if z2 == 0:
                success += 1
            else:
                failure += 1
            continue

        V = np.column_stack(vectors)

        # Compute im(∂₃ restricted to this subspace)
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_V = bd3 @ V  # in A₂ coords

        # Project onto Ω₂ and check if it fills Z₂
        coords, _, _, _ = np.linalg.lstsq(om2, bd3_V, rcond=None)

        # rk of boundary map restricted to our subspace
        rk = np.linalg.matrix_rank(coords, tol=1e-8)

        if rk >= z2:
            success += 1
        else:
            failure += 1
            if failure <= 5:
                t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                         if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))
                scores = tuple(sorted(sum(A[i]) for i in range(n)))
                print(f"  FAILURE: rk={rk}, Z₂={z2}, DT={len(dt_indices)}, cancel_pairs={len(vectors)-len(dt_indices)}, t3={t3}, scores={scores}")

    print(f"\n  RESULTS: success={success}, failure={failure}, total={total}")
    if failure == 0:
        print(f"  ✓ DT + cancellation pairs fill ALL of Z₂ at n={n}!")
    else:
        print(f"  ✗ {failure} tournaments need higher-order elements")

print(f"\n{'='*70}")
print("DONE")
print("="*70)
