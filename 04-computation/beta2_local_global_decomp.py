#!/usr/bin/env python3
"""
beta2_local_global_decomp.py - Test local-to-global decomposition of 2-cycles

Key question: Can every 2-cycle z in ker(d_2|Omega_2) be written as a
sum of elementary cycles, each supported on a single 4-vertex subset?

An elementary cycle on {a,b,c,d} (transitive subtournament with a<b<c<d
in the total order) has the form:
  d_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

This is the boundary of the unique DT 4-path on this transitive 4-tuple.

If every 2-cycle decomposes this way, then im(d_3|DT) = ker(d_2|Omega_2),
which gives beta_2 = 0.

Strategy: compute a basis for Z_2 = ker(d_2|Omega_2) and check if each
basis vector is in span of DT boundaries.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import defaultdict
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

def find_DT_paths(A, n):
    """Find all DT 4-paths: (a,b,c,d) with a->b->c->d, a->c, b->d."""
    dt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if not A[a][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if A[b][d]:
                        dt.append((a,b,c,d))
    return dt

def compute_DT_boundary_matrix(dt_paths, allowed_2, n):
    """Build the boundary matrix from DT 4-paths to allowed 2-paths."""
    if not dt_paths or not allowed_2:
        return np.zeros((len(allowed_2), 0))

    path_to_idx = {p: i for i, p in enumerate(allowed_2)}
    mat = np.zeros((len(allowed_2), len(dt_paths)))

    for j, (a, b, c, d) in enumerate(dt_paths):
        # d_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in path_to_idx:
                mat[path_to_idx[face], j] += sign

    return mat


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"LOCAL-TO-GLOBAL DECOMPOSITION AT n={n}")
print("=" * 70)

success = 0
failure = 0
trivial = 0  # Z_2 = 0

for bits in range(total):
    A = build_adj(n, bits)

    # Get allowed paths
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    # Get Omega_2 basis
    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        trivial += 1
        continue

    # Compute Z_2 = ker(d_2|Omega_2)
    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd2_om = bd2 @ omega2
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    rk2 = int(np.sum(np.abs(S) > 1e-8))
    Z2_dim = dim_O2 - rk2

    if Z2_dim == 0:
        trivial += 1
        continue

    # Get Z_2 basis (in Omega_2 coordinates, then convert to A_2 coordinates)
    Z2_omega = Vt[rk2:].T  # dim_O2 x Z2_dim
    Z2_A2 = omega2 @ Z2_omega  # |A_2| x Z2_dim

    # Get DT 4-path boundary matrix
    dt_paths = find_DT_paths(A, n)
    DT_bd = compute_DT_boundary_matrix(dt_paths, allowed_2, n)

    # Check: is Z_2 in column span of DT_bd?
    if DT_bd.shape[1] == 0:
        if Z2_dim > 0:
            failure += 1
        continue

    # Stack [DT_bd | Z2_A2] and check if rank increases
    combined = np.hstack([DT_bd, Z2_A2])
    rk_DT = np.linalg.matrix_rank(DT_bd, tol=1e-8)
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)

    if rk_combined == rk_DT:
        success += 1
    else:
        failure += 1
        if failure <= 5:
            scores = tuple(sorted(sum(row) for row in A))
            print(f"  FAILURE at bits={bits}: Z2_dim={Z2_dim}, rk(DT_bd)={rk_DT}, rk(combined)={rk_combined}")
            print(f"    scores={scores}, |DT|={len(dt_paths)}")

print(f"\n  n={n}: success={success}, failure={failure}, trivial={trivial}")
print(f"  Z_2 in span(DT boundaries): {'YES' if failure == 0 else 'NO'}")

# Now n=6
n = 6
print(f"\n{'='*70}")
print(f"LOCAL-TO-GLOBAL DECOMPOSITION AT n={n}")
print(f"{'='*70}")

n_arcs = n*(n-1)//2
total = 1 << n_arcs

success = 0
failure = 0
trivial = 0
deficit_1 = 0

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s)")

    A = build_adj(n, bits)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        trivial += 1
        continue

    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd2_om = bd2 @ omega2
    S_v = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S_v) > 1e-8))
    Z2_dim = dim_O2 - rk2

    if Z2_dim == 0:
        trivial += 1
        continue

    # Get DT boundaries
    dt_paths = find_DT_paths(A, n)
    DT_bd = compute_DT_boundary_matrix(dt_paths, allowed_2, n)

    rk_DT = np.linalg.matrix_rank(DT_bd, tol=1e-8)

    if rk_DT >= Z2_dim:
        # DT alone might suffice - check if Z_2 is in span
        U_dt, S_dt, Vt_dt = np.linalg.svd(bd2_om, full_matrices=True)
        Z2_omega = Vt_dt[rk2:].T
        Z2_A2 = omega2 @ Z2_omega

        combined = np.hstack([DT_bd, Z2_A2])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)

        if rk_combined == rk_DT:
            success += 1
        else:
            deficit_1 += 1
    else:
        deficit_1 += 1

dt = time.time() - t0
print(f"\n  n={n} ({dt:.0f}s):")
print(f"  DT alone fills Z_2: {success}/{total} ({100*success/total:.1f}%)")
print(f"  DT deficit: {deficit_1}/{total} ({100*deficit_1/total:.1f}%)")
print(f"  Trivial (Z_2=0): {trivial}/{total} ({100*trivial/total:.1f}%)")

# But the FULL Omega_3 always fills Z_2 (beta_2 = 0).
# The question is: does the DT-only part already span Z_2?

print("\nDone.")
