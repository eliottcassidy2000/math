#!/usr/bin/env python3
"""
beta2_relative_debug.py - Debug the relative H_2 computation

Manual calculation for T = {0->1,1->2,2->0,0->3,1->3,2->3} with v=3 (sink):
  Omega_2 = span{(0,1,3),(1,2,3),(2,0,3)}, dim=3
  Omega_2(T\v) = 0
  Omega_3 = 0

  d_2 on Omega_2 projected to v-arcs:
    (0,1,3) -> (1,3)-(0,3) in v-arcs
    (1,2,3) -> (2,3)-(1,3) in v-arcs
    (2,0,3) -> (0,3)-(2,3) in v-arcs
  Rank = 2, ker = 1

  ker(d_2^rel) = span{(0,1,3)+(1,2,3)+(2,0,3)}, dim=1
  im(d_3^rel) = 0

  H_2^rel = 1 (NOT 0!)

But beta2_relative_homology.py reported H_2^rel = 0 at n=4.

The BUG: the script computed ker(d_2^rel) as (ker d_2 + V_2)/V_2,
but the CORRECT formula is:
  ker(d_2^rel) = {x in Omega_2 : d_2(x) in V_1} / V_2
where V_1 = non-v arcs.

This is LARGER than (ker d_2 + V_2)/V_2 because it includes elements
whose boundary is nonzero but lies entirely in non-v arcs.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
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

def compute_relative_h2_CORRECT(A, n, v, verbose=False):
    """CORRECT relative H_2 computation.

    ker(d_2^rel) = dim{x in Omega_2 : d_2(x) in V_1} - dim(V_2)
    where V_1 = Omega_1(T\\v) = non-v arcs, V_2 = Omega_2(T\\v).

    Key: the old code computed (ker d_2 + V_2)/V_2 which MISSES elements
    whose boundary is nonzero but lands in V_1.
    """
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dim_om2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
    dim_om3 = omega[3].shape[1] if omega[3].ndim == 2 else 0

    if dim_om2 == 0:
        return 0

    # Boundary d_2: A_2 -> A_1, restricted to Omega_2
    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]  # |A_1| x dim_om2

    # V_1 = non-v arcs (indices in A_1)
    v_arc_idx = [i for i in range(len(paths[1])) if v in paths[1][i]]
    nv_arc_idx = [i for i in range(len(paths[1])) if v not in paths[1][i]]

    # The "relative boundary" = projection of d_2 onto v-arcs
    # ker(d_2^rel) = ker of this projection restricted to Omega_2
    M = bd2_om[v_arc_idx, :]  # |v-arcs| x dim_om2
    sv_M = np.linalg.svd(M, compute_uv=False)
    rk_M = int(np.sum(np.abs(sv_M) > 1e-8))
    dim_preimage = dim_om2 - rk_M  # dim{x in Omega_2 : d_2(x) in V_1}

    # V_2 = Omega_2 ∩ no-v paths
    v_path_idx_2 = [i for i in range(len(paths[2])) if v in paths[2][i]]
    if v_path_idx_2:
        R_v2 = omega[2][v_path_idx_2, :]
        sv_rv2 = np.linalg.svd(R_v2, compute_uv=False)
        rk_rv2 = int(np.sum(np.abs(sv_rv2) > 1e-8))
        dim_V2 = dim_om2 - rk_rv2
    else:
        dim_V2 = dim_om2

    ker2_rel = dim_preimage - dim_V2

    if verbose:
        print(f"    dim_om2={dim_om2}, rk_M={rk_M}, dim_preimage={dim_preimage}, dim_V2={dim_V2}")
        print(f"    ker(d_2^rel) = {ker2_rel}")

    if ker2_rel <= 0:
        return 0

    # im(d_3^rel): image of d_3 in quotient Omega_2/V_2
    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]  # |A_2| x dim_om3

        # Express in Omega_2 coordinates
        im3_om2, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)

        # V_2 basis in Omega_2 coordinates
        if v_path_idx_2 and rk_rv2 < dim_om2:
            Uv, Sv, Vtv = np.linalg.svd(omega[2][v_path_idx_2, :], full_matrices=True)
            V2_basis = Vtv[rk_rv2:].T
        elif not v_path_idx_2:
            V2_basis = np.eye(dim_om2)
        else:
            V2_basis = np.zeros((dim_om2, 0))

        if V2_basis.shape[1] > 0:
            combined = np.hstack([im3_om2, V2_basis])
            rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
            im3_rel = rk_combined - V2_basis.shape[1]
        else:
            im3_rel = np.linalg.matrix_rank(im3_om2, tol=1e-8)
    else:
        im3_rel = 0

    if verbose:
        print(f"    dim_om3={dim_om3}, im(d_3^rel)={im3_rel}")
        print(f"    H_2^rel = {max(0, ker2_rel - im3_rel)}")

    return max(0, ker2_rel - im3_rel)


# ===== Test specific case =====
print("=" * 70)
print("DEBUG: CORRECT RELATIVE H_2 COMPUTATION")
print("=" * 70)

# Tournament: 0->1, 1->2, 2->0, 0->3, 1->3, 2->3 (sink at 3)
# bits = 1 + 4 + 8 + 16 + 32 = 61
n = 4
bits = 61
A = build_adj(n, bits)
print(f"\nTournament bits={bits} (n={n}):")
print(f"  Adjacency:")
for i in range(n):
    print(f"    {i}: {[j for j in range(n) if A[i][j]]}")
print(f"  Scores: {[sum(row) for row in A]}")

for v_test in range(n):
    print(f"\n  v={v_test} (d_v={sum(A[v_test])})")
    h2_rel = compute_relative_h2_CORRECT(A, n, v_test, verbose=True)
    print(f"    => H_2^rel = {h2_rel}")

# ===== Exhaustive n=3,4,5 =====
for n_test in [3, 4, 5]:
    print(f"\n{'='*70}")
    print(f"EXHAUSTIVE n={n_test}")
    print("=" * 70)
    n_arcs = n_test*(n_test-1)//2
    total = 1 << n_arcs

    nonzero_count = 0
    nonzero_examples = []

    for bits_t in range(total):
        A_t = build_adj(n_test, bits_t)
        for v_t in range(n_test):
            h2 = compute_relative_h2_CORRECT(A_t, n_test, v_t)
            if h2 > 0:
                nonzero_count += 1
                if len(nonzero_examples) < 5:
                    scores = tuple(sorted(sum(row) for row in A_t))
                    d_v_t = sum(A_t[v_t])
                    nonzero_examples.append((bits_t, v_t, d_v_t, scores, h2))

    print(f"  Total (T, v) pairs: {total * n_test}")
    print(f"  H_2^rel > 0: {nonzero_count}")
    if nonzero_examples:
        for bits_t, v_t, dv, sc, h2 in nonzero_examples:
            print(f"    bits={bits_t}, v={v_t}, d_v={dv}, scores={sc}, H_2^rel={h2}")

# ===== n=6 exhaustive =====
print(f"\n{'='*70}")
print(f"EXHAUSTIVE n=6")
print("=" * 70)
n_test = 6
n_arcs = n_test*(n_test-1)//2
total = 1 << n_arcs
nonzero_count = 0

import time
t0 = time.time()
for bits_t in range(total):
    if bits_t % 5000 == 0 and bits_t > 0:
        dt = time.time() - t0
        print(f"  ... {bits_t}/{total} ({dt:.0f}s), nonzero={nonzero_count}")
    A_t = build_adj(n_test, bits_t)
    for v_t in range(n_test):
        h2 = compute_relative_h2_CORRECT(A_t, n_test, v_t)
        if h2 > 0:
            nonzero_count += 1
            if nonzero_count <= 3:
                scores = tuple(sorted(sum(row) for row in A_t))
                d_v_t = sum(A_t[v_t])
                print(f"    *** bits={bits_t}, v={v_t}, d_v={d_v_t}, scores={scores}, H_2^rel={h2}")

dt = time.time() - t0
print(f"\n  Total (T, v) pairs: {total * n_test}")
print(f"  H_2^rel > 0: {nonzero_count} ({dt:.0f}s)")

print("\nDone.")
