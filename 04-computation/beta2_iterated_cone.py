#!/usr/bin/env python3
"""
beta2_iterated_cone.py — Residual analysis of partial cone

For non-source tournaments: what does the residual r = z - d3(c_v(z_+)) look like?
Is it in Omega_2? In Z_2? In B_2?

Author: opus-2026-03-08-S49
"""
import sys, time
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

print("=" * 70)
print("RESIDUAL OF PARTIAL CONE")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

r_in_omega2 = Counter()
r_in_z2 = Counter()
r_in_b2 = Counter()

count = 0
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    v = scores.index(max(scores))
    if scores[v] == n-1:
        continue  # Skip source

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk
    if z2_dim == 0:
        continue

    z2_A2 = om2 @ Vt[rk:, :].T
    ap3_list = [tuple(p) for p in ap3]
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    for j in range(z2_dim):
        z = z2_A2[:, j]
        cv_z = np.zeros(len(ap3_list))
        for i, path in enumerate(ap2):
            if abs(z[i]) < 1e-10 or v in path or not A[v][path[0]]:
                continue
            cone_path = tuple([v] + list(path))
            if cone_path in ap3_list:
                cv_z[ap3_list.index(cone_path)] += z[i]

        d3_cv = bd3 @ cv_z
        r = z - d3_cv
        r_norm = np.linalg.norm(r)

        if r_norm < 1e-8:
            r_in_omega2[True] += 1
            r_in_z2[True] += 1
            r_in_b2[True] += 1
            continue

        # Is r in Omega_2?
        r_coords = np.linalg.lstsq(om2, r, rcond=None)[0]
        in_om2 = np.linalg.norm(r - om2 @ r_coords) < 1e-6
        r_in_omega2[in_om2] += 1

        if in_om2:
            # Is r in Z_2?
            d2_r = bd2 @ r
            # Need d2_r in Omega_1 and equal to 0
            d2_r_in_om1 = np.linalg.lstsq(om1, d2_r, rcond=None)[0]
            d2_r_check = np.linalg.norm(bd2 @ r)
            # Actually just check if boundary is 0
            # But r might not be in Omega_2 properly. Let me project.
            r_om2_coords = np.linalg.lstsq(om1, bd2 @ (om2 @ r_coords), rcond=None)[0]
            in_z2 = np.linalg.norm(r_om2_coords) < 1e-6
            r_in_z2[in_z2] += 1

            if in_z2 and d3 > 0:
                # Is r in B_2?
                B2 = bd3 @ om3
                sigma = np.linalg.lstsq(B2, r, rcond=None)[0]
                in_b2 = np.linalg.norm(B2 @ sigma - r) < 1e-6
                r_in_b2[in_b2] += 1
            else:
                r_in_b2[False] += 1
        else:
            r_in_z2['N/A'] += 1
            r_in_b2['N/A'] += 1

    count += 1

print(f"\nn=5: {count} non-source tournaments analyzed")
print(f"  r in Omega_2: {dict(r_in_omega2)}")
print(f"  r in Z_2: {dict(r_in_z2)}")
print(f"  r in B_2: {dict(r_in_b2)}")

# KEY: if residual is ALWAYS in B_2, then z = d3(c_v(z_+)) + r
# and both terms are in B_2, so z in B_2. QED.

print("\nDone.")
