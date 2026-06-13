#!/usr/bin/env python3
"""
beta2_cone_omega_preserve.py — Does the cone c_v map Ω₂ to Ω₃?

The chain homotopy identity ∂₃ ∘ c_v + c_v ∘ ∂₂ = id works on A_* (allowed paths).
If c_v also maps Ω_p to Ω_{p+1}, then for z ∈ Z₂ ⊆ Ω₂:
  ∂₃(c_v(z)) + c_v(∂₂(z)) = z
  ∂₃(c_v(z)) = z   (since ∂₂(z)=0)
and c_v(z) ∈ Ω₃, so z ∈ B₂. Hence β₂ = 0.

Test: for v = source vertex, does c_v: Ω₂ → Ω₃?
More generally, for any v, does c_v: Z₂ → Ω₃? (weaker but sufficient)

CRITICAL: also test whether ∂₃(c_v(z)) = z holds at the A₂ level
(without needing c_v(z) ∈ Ω₃).

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
print("CONE OMEGA PRESERVATION TEST")
print("=" * 70)

# Test 1: Does chain homotopy ∂₃(c_v(z)) = z hold at A₂ level for source v?
print("\nTest 1: Chain homotopy ∂₃(c_v(z)) = z (A₂ level, source v)")

n = 5
m = n*(n-1)//2
total = 1 << m

homotopy_works = 0
homotopy_fails = 0
omega_preserve = 0
omega_no_preserve = 0
z2_trivial = 0

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    # Find source-like vertex (max out-degree)
    v = scores.index(max(scores))

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    if d2 == 0:
        z2_trivial += 1
        continue

    # Z₂ in Ω₂ coords
    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        z2_trivial += 1
        continue

    z2_om = Vt[rk:, :]
    z2_A2 = om2 @ z2_om.T  # |A₂| × z2_dim

    ap2_list = [tuple(p) for p in ap2]
    ap3_list = [tuple(p) for p in ap3]

    # For each Z₂ basis vector: compute c_v(z) and ∂₃(c_v(z))
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    all_homotopy_ok = True
    all_in_omega3 = True

    for j in range(z2_dim):
        z = z2_A2[:, j]

        # c_v(z) in A₃ coordinates
        cv_z = np.zeros(len(ap3_list))
        for i, path in enumerate(ap2):
            if abs(z[i]) < 1e-10:
                continue
            if v in path:
                continue
            if not A[v][path[0]]:
                continue
            cone_path = tuple([v] + list(path))
            if cone_path in ap3_list:
                cv_z[ap3_list.index(cone_path)] += z[i]

        # ∂₃(c_v(z)) in A₂ coordinates
        d3_cv = bd3 @ cv_z

        # Check if ∂₃(c_v(z)) = z at A₂ level
        diff = d3_cv - z
        if np.linalg.norm(diff) > 1e-6:
            all_homotopy_ok = False

        # Check if c_v(z) ∈ Ω₃
        d3 = dim_om(om3)
        if d3 > 0:
            coords = np.linalg.lstsq(om3, cv_z, rcond=None)[0]
            resid = cv_z - om3 @ coords
            if np.linalg.norm(resid) > 1e-6:
                all_in_omega3 = False
        elif np.linalg.norm(cv_z) > 1e-6:
            all_in_omega3 = False

    if all_homotopy_ok:
        homotopy_works += 1
    else:
        homotopy_fails += 1

    if all_in_omega3:
        omega_preserve += 1
    else:
        omega_no_preserve += 1

print(f"\nn={n}: {total} tournaments")
print(f"  Z₂ trivial: {z2_trivial}")
print(f"  Chain homotopy ∂₃(c_v(z))=z works (max outdegree v): {homotopy_works}")
print(f"  Chain homotopy fails: {homotopy_fails}")
print(f"  c_v(Z₂) ⊆ Ω₃: {omega_preserve}")
print(f"  c_v(Z₂) ⊄ Ω₃: {omega_no_preserve}")

# Test 2: For SOURCE vertices (d+=n-1)
print(f"\n{'='*70}")
print("Test 2: SOURCE VERTEX (d+=n-1) CONE")
print("=" * 70)

n = 5
src_homotopy_ok = 0
src_homotopy_fail = 0
src_omega_ok = 0
src_omega_fail = 0
has_source = 0
no_source = 0

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    # Find actual source (d+=n-1)
    sources = [i for i in range(n) if scores[i] == n-1]
    if not sources:
        no_source += 1
        continue

    has_source += 1
    v = sources[0]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
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

    all_ok = True
    all_in_om3 = True

    for j in range(z2_dim):
        z = z2_A2[:, j]

        cv_z = np.zeros(len(ap3_list))
        for i, path in enumerate(ap2):
            if abs(z[i]) < 1e-10:
                continue
            if v in path:
                continue
            # Source: v→path[0] always true
            cone_path = tuple([v] + list(path))
            if cone_path in ap3_list:
                cv_z[ap3_list.index(cone_path)] += z[i]

        d3_cv = bd3 @ cv_z
        diff = d3_cv - z
        if np.linalg.norm(diff) > 1e-6:
            all_ok = False

        d3 = dim_om(om3)
        if d3 > 0:
            coords = np.linalg.lstsq(om3, cv_z, rcond=None)[0]
            resid = cv_z - om3 @ coords
            if np.linalg.norm(resid) > 1e-6:
                all_in_om3 = False
        elif np.linalg.norm(cv_z) > 1e-6:
            all_in_om3 = False

    if all_ok:
        src_homotopy_ok += 1
    else:
        src_homotopy_fail += 1

    if all_in_om3:
        src_omega_ok += 1
    else:
        src_omega_fail += 1

print(f"\nn={n}: {has_source} with source, {no_source} without")
print(f"  Source cone: homotopy ok={src_homotopy_ok}, fail={src_homotopy_fail}")
print(f"  Source cone: c_v(Z₂) ⊆ Ω₃ = {src_omega_ok}, ⊄ Ω₃ = {src_omega_fail}")

# Test 3: CRITICAL - does the chain homotopy at A level imply β₂=0?
# If ∂₃(c_v(z)) = z in A₂, and z ∈ Ω₂, then z ∈ im(∂₃|_{A₃}).
# But we need z ∈ im(∂₃|_{Ω₃}). Is im(∂₃|_{A₃}) ∩ Ω₂ = im(∂₃|_{Ω₃})?
print(f"\n{'='*70}")
print("Test 3: im(∂₃|A₃) ∩ Ω₂ vs im(∂₃|Ω₃)")
print("=" * 70)

n = 5
match_count = 0
mismatch_count = 0

for bits in range(total):
    A = build_adj(n, bits)

    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)

    if d2 == 0:
        continue

    bd3_full = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    # im(∂₃|Ω₃) in A₂: bd3_full @ om3
    if d3 > 0:
        im_omega3 = bd3_full @ om3  # columns give image
        # Project to Ω₂: om2^+ @ im_omega3
        im_om3_in_om2 = np.linalg.lstsq(om2, im_omega3, rcond=None)[0]
        rk_omega = np.linalg.matrix_rank(im_om3_in_om2, tol=1e-8)
    else:
        rk_omega = 0

    # im(∂₃|A₃) in A₂: bd3_full (all columns)
    if len(ap3) > 0:
        # Project ∂₃(A₃) to Ω₂
        im_A3 = bd3_full  # all of ∂₃(A₃)
        im_A3_in_om2 = np.linalg.lstsq(om2, im_A3, rcond=None)[0]
        # Residual check: is im_A3 in Ω₂?
        resid = im_A3 - om2 @ im_A3_in_om2
        # Only keep columns that project into Ω₂
        in_omega2_mask = np.array([np.linalg.norm(resid[:, j]) < 1e-6 for j in range(im_A3.shape[1])])
        im_A3_in_om2_filtered = im_A3_in_om2[:, in_omega2_mask]
        rk_A3_filtered = np.linalg.matrix_rank(im_A3_in_om2_filtered, tol=1e-8) if im_A3_in_om2_filtered.shape[1] > 0 else 0
    else:
        rk_A3_filtered = 0

    if rk_omega == rk_A3_filtered:
        match_count += 1
    else:
        mismatch_count += 1

print(f"\nn={n}: im(∂₃|Ω₃)∩Ω₂ == im(∂₃|A₃)∩Ω₂?")
print(f"  Match: {match_count}, Mismatch: {mismatch_count}")

# Test 4: Does ∂₃(A₃) ∩ Ω₂ = ∂₃(Ω₃)?
# (Without restricting to Z₂)
print(f"\nTest 4: rank of ∂₃ restricted to Ω₃ vs projected from A₃")
print(f"  rk_omega == rk_A3_filtered for all tournaments: {mismatch_count == 0}")

print("\nDone.")
