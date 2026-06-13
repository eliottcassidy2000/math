#!/usr/bin/env python3
"""
beta2_missing_dir.py — Find the Z₂ direction that DT boundaries miss at n=6

960 tournaments at n=6 have DT rank = 9 but Z₂ = 10.
Find the missing 2-cycle and understand its structure.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print("=" * 70)
print("MISSING Z₂ DIRECTION ANALYSIS")
print("=" * 70)

bits = 154
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(pairs):
    if (bits >> idx) & 1: A[i][j] = 1
    else: A[j][i] = 1

print(f"\nT#{bits}:")
for i in range(n):
    out = [j for j in range(n) if A[i][j]]
    print(f"  {i} → {out}")

ap0 = enumerate_allowed_paths(A, n, 0)
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

om1 = compute_omega_basis(A, n, 1, ap1, ap0)
om2 = compute_omega_basis(A, n, 2, ap2, ap1)
om3 = compute_omega_basis(A, n, 3, ap3, ap2)

d2 = om2.shape[1]
d3 = om3.shape[1]

# Z₂ basis
bd2 = build_full_boundary_matrix(ap2, ap1)
bd2_om = bd2 @ om2
coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
rk2 = int(sum(s > 1e-8 for s in S2))
z2_basis_om2 = Vt2[rk2:]  # (z2_dim) × d2, rows are Z₂ basis in Ω₂ coords
z2_dim = z2_basis_om2.shape[0]
print(f"\nZ₂ dim = {z2_dim}, Ω₂ dim = {d2}, Ω₃ dim = {d3}")

# DT paths
ap3_t = [tuple(p) for p in ap3]
dt_idx = [i for i, (a,b,c,d) in enumerate(ap3_t) if A[a][c] and A[b][d]]
ndt_idx = [i for i in range(len(ap3_t)) if i not in dt_idx]

# Boundary in Ω₂ coords
bd3 = build_full_boundary_matrix(ap3, ap2)

# DT boundaries in Ω₂ coords
dt_bd_A2 = bd3[:, dt_idx]
dt_bd_om2 = np.linalg.lstsq(om2, dt_bd_A2, rcond=None)[0]

# DT boundaries projected into Z₂ space
dt_in_z2 = z2_basis_om2 @ dt_bd_om2  # z2_dim × #DT
print(f"DT image in Z₂: shape {dt_in_z2.shape}")

# Find missing direction: left null space complement
U_dt, S_dt, _ = np.linalg.svd(dt_in_z2, full_matrices=True)
rk_dt = int(sum(s > 1e-8 for s in S_dt))
print(f"DT rank in Z₂: {rk_dt}")

# Missing direction in Z₂ coords
missing_z2_coords = U_dt[:, rk_dt:]  # z2_dim × (z2_dim - rk_dt)
print(f"Missing directions: {missing_z2_coords.shape[1]}")

# Convert to Ω₂ coords
missing_om2 = z2_basis_om2.T @ missing_z2_coords  # d2 × 1

# Convert to A₂ coords
missing_a2 = om2 @ missing_om2  # |A₂| × 1

print(f"\nMissing Z₂ vector (in A₂ coords):")
for col in range(missing_a2.shape[1]):
    v = missing_a2[:, col]
    terms = [(tuple(ap2[i]), round(v[i], 5)) for i in range(len(v)) if abs(v[i]) > 1e-8]
    for path, coeff in terms:
        a, b, c = path
        tt = A[a][c]
        print(f"  {coeff:+.5f} * ({a},{b},{c}) [{'TT' if tt else 'NT'}]")

# Verify it's a cycle
print(f"\nVerify ∂₂ of missing vector = 0:")
bd2_missing = bd2 @ missing_a2
for col in range(bd2_missing.shape[1]):
    residual = np.max(np.abs(bd2_missing[:, col]))
    print(f"  ||∂₂(missing)||_inf = {residual:.2e}")

# Now: which Ω₃ element kills this cycle?
print(f"\n{'='*70}")
print("WHICH Ω₃ ELEMENT FILLS THE GAP?")
print("=" * 70)

# Full Ω₃ boundary in Z₂ space
full_bd_om2 = np.linalg.lstsq(om2, bd3 @ om3, rcond=None)[0]
full_in_z2 = z2_basis_om2 @ full_bd_om2  # z2_dim × d3

# Project missing direction onto full Ω₃ image
proj = missing_z2_coords.T @ full_in_z2  # 1 × d3
print(f"Missing direction's projection onto Ω₃ basis:")
for j in range(d3):
    if abs(proj[0, j]) > 1e-8:
        print(f"  Ω₃ basis e{j}: {proj[0,j]:+.5f}")
        # Show this Ω₃ element
        v = om3[:, j]
        terms = [(ap3_t[i], round(v[i], 4)) for i in range(len(v)) if abs(v[i]) > 1e-8]
        dt_count = sum(1 for p, _ in terms if A[p[0]][p[2]] and A[p[1]][p[3]])
        ndt_count = len(terms) - dt_count
        print(f"    {len(terms)} terms ({dt_count} DT, {ndt_count} non-DT)")

# Check: what types of non-DT paths appear?
print(f"\n{'='*70}")
print("CLASSIFYING NON-DT PATHS")
print("=" * 70)

# A non-DT path (a,b,c,d) has a→b→c→d, but NOT (a→c AND b→d)
# Types:
#   ST (semi-transitive-1): a→c, NOT b→d (i.e., d→b)
#   TS (semi-transitive-2): NOT a→c (i.e., c→a), b→d
#   OT (anti-transitive): NOT a→c AND NOT b→d
for i, (a,b,c,d) in enumerate(ap3_t):
    if i in dt_idx:
        continue
    ac = A[a][c]
    bd = A[b][d]
    if ac and not bd:
        typ = "ST1(a→c, d→b)"
    elif not ac and bd:
        typ = "ST2(c→a, b→d)"
    elif not ac and not bd:
        typ = "OT(c→a, d→b)"
    else:
        typ = "??"  # shouldn't happen

# Count types across ALL failures
print(f"\n{'='*70}")
print("NON-DT TYPE DISTRIBUTION ACROSS FAILURES")
print("=" * 70)

from collections import Counter
type_dist = Counter()
ndt_in_om3 = Counter()

for bits_test in range(1 << m):
    A_t = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits_test >> idx) & 1: A_t[i][j] = 1
        else: A_t[j][i] = 1

    ap3_l = enumerate_allowed_paths(A_t, n, 3)
    if not ap3_l:
        continue

    ap2_l = enumerate_allowed_paths(A_t, n, 2)
    ap1_l = enumerate_allowed_paths(A_t, n, 1)
    ap0_l = enumerate_allowed_paths(A_t, n, 0)
    if not ap2_l:
        continue

    om1_l = compute_omega_basis(A_t, n, 1, ap1_l, ap0_l)
    om2_l = compute_omega_basis(A_t, n, 2, ap2_l, ap1_l)
    d2_l = om2_l.shape[1] if om2_l.ndim == 2 and om2_l.shape[0] > 0 else 0
    if d2_l == 0:
        continue

    bd2_l = build_full_boundary_matrix(ap2_l, ap1_l)
    bd2_om_l = bd2_l @ om2_l
    coords2_l = np.linalg.lstsq(om1_l, bd2_om_l, rcond=None)[0]
    rk2_l = np.linalg.matrix_rank(coords2_l, tol=1e-8)
    z2_l = d2_l - rk2_l

    if z2_l == 0:
        continue

    ap3_t_l = [tuple(p) for p in ap3_l]
    dt_idx_l = [i for i, (a,b,c,d) in enumerate(ap3_t_l)
                if A_t[a][c] and A_t[b][d]]

    bd3_l = build_full_boundary_matrix(ap3_l, ap2_l)
    dt_bd_l = bd3_l[:, dt_idx_l] if dt_idx_l else np.zeros((len(ap2_l), 0))
    dt_om2_l = np.linalg.lstsq(om2_l, dt_bd_l, rcond=None)[0] if dt_idx_l else np.zeros((d2_l, 0))
    rk_dt_l = np.linalg.matrix_rank(dt_om2_l, tol=1e-8) if dt_idx_l else 0

    is_fail = rk_dt_l < z2_l
    if is_fail:
        # Count non-DT types
        for i, (a,b,c,d) in enumerate(ap3_t_l):
            if i in dt_idx_l:
                continue
            ac = A_t[a][c]
            bd = A_t[b][d]
            if ac and not bd:
                type_dist['ST1'] += 1
            elif not ac and bd:
                type_dist['ST2'] += 1
            else:
                type_dist['OT'] += 1

        # Check if non-DT paths are in Ω₃
        om3_l = compute_omega_basis(A_t, n, 3, ap3_l, ap2_l)
        d3_l = om3_l.shape[1] if om3_l.ndim == 2 and om3_l.shape[0] > 0 else 0
        if d3_l > 0:
            for col in range(d3_l):
                v = om3_l[:, col]
                has_dt = any(abs(v[j]) > 1e-8 for j in dt_idx_l)
                has_ndt = any(abs(v[j]) > 1e-8 for j in range(len(ap3_t_l)) if j not in dt_idx_l)
                if has_ndt and not has_dt:
                    ndt_in_om3['pure_ndt'] += 1
                elif has_ndt and has_dt:
                    ndt_in_om3['mixed'] += 1
                else:
                    ndt_in_om3['pure_dt'] += 1

print(f"Non-DT type distribution in failures: {dict(type_dist)}")
print(f"Ω₃ basis composition in failures: {dict(ndt_in_om3)}")

print("\nDone.")
