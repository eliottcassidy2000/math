#!/usr/bin/env python3
"""
beta2_proof_path_n7.py — Verify both proof ingredients at n=7,8

Ingredient 1 (HYP-269): Z₂(Ω) ⊆ im(∂₃|A₃) for all tournaments
Ingredient 2 (HYP-270): im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃)

Together: Z₂(Ω) ⊆ im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃) = B₂(Ω), so β₂ = 0.

Author: opus-2026-03-08-S49
"""
import sys, time, random
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

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


print("=" * 70)
print("PROOF PATH VERIFICATION AT n=7,8")
print("=" * 70)

for n, samples in [(7, 500), (8, 100)]:
    t0 = time.time()

    hyp269_ok = 0
    hyp269_fail = 0
    hyp270_ok = 0
    hyp270_fail = 0

    for trial in range(samples):
        A = build_random_adj(n)

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
            hyp269_ok += 1
            hyp270_ok += 1
            continue

        # Z₂(Ω)
        bd2 = build_full_boundary_matrix(ap2, ap1)
        d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z2_dim = d2 - rk

        if z2_dim == 0:
            hyp269_ok += 1
            hyp270_ok += 1
            continue

        z2_A2 = om2 @ Vt[rk:, :].T

        # ∂₃: A₃ → A₂
        bd3_full = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

        # HYP-269: Z₂(Ω) ⊆ im(∂₃|A₃)?
        if bd3_full.shape[1] > 0:
            all_in = True
            for j in range(z2_dim):
                z = z2_A2[:, j]
                sigma = np.linalg.lstsq(bd3_full, z, rcond=None)[0]
                resid = np.linalg.norm(bd3_full @ sigma - z)
                if resid > 1e-6:
                    all_in = False
                    break
            if all_in:
                hyp269_ok += 1
            else:
                hyp269_fail += 1
        else:
            hyp269_fail += 1

        # HYP-270: im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃)?
        if d3 > 0 and bd3_full.shape[1] > 0:
            # im(∂₃|Ω₃): project bd3 @ om3 to Ω₂ coords
            im_omega = np.linalg.lstsq(om2, bd3_full @ om3, rcond=None)[0]
            rk_omega = np.linalg.matrix_rank(im_omega, tol=1e-8)

            # im(∂₃|A₃) ∩ Ω₂: project all bd3 columns to Ω₂, keep those in Ω₂
            im_A3 = bd3_full
            im_A3_om2 = np.linalg.lstsq(om2, im_A3, rcond=None)[0]
            resid = im_A3 - om2 @ im_A3_om2
            in_omega2 = np.array([np.linalg.norm(resid[:, j]) < 1e-6 for j in range(resid.shape[1])])
            if in_omega2.any():
                rk_A3 = np.linalg.matrix_rank(im_A3_om2[:, in_omega2], tol=1e-8)
            else:
                rk_A3 = 0

            if rk_omega == rk_A3:
                hyp270_ok += 1
            else:
                hyp270_fail += 1
        elif d3 == 0 and bd3_full.shape[1] > 0:
            # Ω₃ = 0 but A₃ ≠ 0: check if any ∂₃(A₃) lands in Ω₂
            im_A3 = bd3_full
            im_A3_om2 = np.linalg.lstsq(om2, im_A3, rcond=None)[0]
            resid = im_A3 - om2 @ im_A3_om2
            in_omega2 = np.array([np.linalg.norm(resid[:, j]) < 1e-6 for j in range(resid.shape[1])])
            rk_A3 = np.linalg.matrix_rank(im_A3_om2[:, in_omega2], tol=1e-8) if in_omega2.any() else 0
            if rk_A3 == 0:
                hyp270_ok += 1
            else:
                hyp270_fail += 1
        else:
            hyp270_ok += 1

        if (trial+1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{samples} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {samples} tournaments in {elapsed:.0f}s")
    print(f"  HYP-269 (Z₂ ⊆ im∂₃A₃): ok={hyp269_ok}, fail={hyp269_fail}")
    print(f"  HYP-270 (im∂₃A₃∩Ω₂=im∂₃Ω₃): ok={hyp270_ok}, fail={hyp270_fail}")

    if hyp269_fail == 0 and hyp270_fail == 0:
        print(f"  *** BOTH ingredients verified for all n={n} samples! ***")

print("\nDone.")
