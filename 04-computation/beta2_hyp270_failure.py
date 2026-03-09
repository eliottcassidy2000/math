#!/usr/bin/env python3
"""
beta2_hyp270_failure.py — Investigate HYP-270 failures at n=7

HYP-270: im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃)
This FAILS at n=7 for some tournaments. Is it related to β₃ > 0?

Also verify HYP-270 at n=6 (exhaustive) and characterize failures at n=7.

Author: opus-2026-03-08-S49
"""
import sys, time, random
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

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_betti(A, n, max_p=None):
    """Compute β_0,...,β_{max_p} for tournament A."""
    if max_p is None:
        max_p = n-1
    ap = {}
    om = {}
    for p in range(max_p+2):
        ap[p] = enumerate_allowed_paths(A, n, p)

    om[0] = np.eye(n)
    for p in range(1, max_p+2):
        if ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0,0))

    betti = []
    for p in range(max_p+1):
        dp = dim_om(om[p]) if p > 0 else n
        dp1 = dim_om(om[p+1])

        if dp == 0:
            betti.append(0)
            continue

        if p == 0:
            bd_p = np.zeros((1, n))  # trivial
            rk_p = 0
        else:
            bd_p = build_full_boundary_matrix(ap[p], ap[p-1])
            d_mat = np.linalg.lstsq(om[p-1] if p > 1 else np.eye(n),
                                     bd_p @ om[p], rcond=None)[0]
            rk_p = np.linalg.matrix_rank(d_mat, tol=1e-8)

        ker_p = dp - rk_p

        if dp1 > 0:
            bd_p1 = build_full_boundary_matrix(ap[p+1], ap[p])
            d_mat1 = np.linalg.lstsq(om[p], bd_p1 @ om[p+1], rcond=None)[0]
            rk_p1 = np.linalg.matrix_rank(d_mat1, tol=1e-8)
        else:
            rk_p1 = 0

        betti.append(ker_p - rk_p1)

    return betti


# Test 1: HYP-270 at n=6 (exhaustive)
print("=" * 70)
print("HYP-270 at n=6 (EXHAUSTIVE)")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

hyp270_ok = 0
hyp270_fail = 0

for bits in range(total):
    A = build_adj(n, bits)

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap1 = enumerate_allowed_paths(A, n, 1)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    bd3_full = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    if d3 > 0 and bd3_full.shape[1] > 0:
        im_omega = np.linalg.lstsq(om2, bd3_full @ om3, rcond=None)[0]
        rk_omega = np.linalg.matrix_rank(im_omega, tol=1e-8)

        im_A3_om2 = np.linalg.lstsq(om2, bd3_full, rcond=None)[0]
        resid = bd3_full - om2 @ im_A3_om2
        in_omega2 = np.array([np.linalg.norm(resid[:, j]) < 1e-6 for j in range(resid.shape[1])])
        rk_A3 = np.linalg.matrix_rank(im_A3_om2[:, in_omega2], tol=1e-8) if in_omega2.any() else 0

        if rk_omega == rk_A3:
            hyp270_ok += 1
        else:
            hyp270_fail += 1
    else:
        hyp270_ok += 1

    if (bits+1) % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s) ok={hyp270_ok} fail={hyp270_fail}")

elapsed = time.time() - t0
print(f"\nn=6: {total} tournaments in {elapsed:.0f}s")
print(f"  HYP-270 ok: {hyp270_ok}, fail: {hyp270_fail}")

# Test 2: HYP-270 failures at n=7 — correlate with β₃
print(f"\n{'='*70}")
print("HYP-270 FAILURES AT n=7: CORRELATION WITH β₃")
print("=" * 70)

n = 7
samples = 1000
t0 = time.time()

fail_b3 = Counter()
ok_b3 = Counter()

for trial in range(samples):
    A = build_random_adj(n)

    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap4 = enumerate_allowed_paths(A, n, 4)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    bd3_full = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    # HYP-270 check
    if d3 > 0 and bd3_full.shape[1] > 0:
        im_omega = np.linalg.lstsq(om2, bd3_full @ om3, rcond=None)[0]
        rk_omega = np.linalg.matrix_rank(im_omega, tol=1e-8)

        im_A3_om2 = np.linalg.lstsq(om2, bd3_full, rcond=None)[0]
        resid = bd3_full - om2 @ im_A3_om2
        in_omega2 = np.array([np.linalg.norm(resid[:, j]) < 1e-6 for j in range(resid.shape[1])])
        rk_A3 = np.linalg.matrix_rank(im_A3_om2[:, in_omega2], tol=1e-8) if in_omega2.any() else 0

        hyp270_holds = (rk_omega == rk_A3)
    else:
        hyp270_holds = True

    # Compute β₃
    om4 = compute_omega_basis(A, n, 4, ap4, ap3) if ap4 else np.zeros((0,0))
    d4 = dim_om(om4)

    # ker(∂₃|Ω₃)
    if d3 > 0:
        bd3_om = np.linalg.lstsq(om2, bd3_full @ om3, rcond=None)[0]
        rk3 = np.linalg.matrix_rank(bd3_om, tol=1e-8)
        ker3 = d3 - rk3
    else:
        ker3 = 0

    # im(∂₄|Ω₄)
    if d4 > 0 and d3 > 0:
        bd4 = build_full_boundary_matrix(ap4, ap3)
        bd4_om = np.linalg.lstsq(om3, bd4 @ om4, rcond=None)[0]
        rk4 = np.linalg.matrix_rank(bd4_om, tol=1e-8)
    else:
        rk4 = 0

    b3 = ker3 - rk4

    if hyp270_holds:
        ok_b3[b3] += 1
    else:
        fail_b3[b3] += 1

    if (trial+1) % 200 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{samples} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=7: {samples} tournaments in {elapsed:.0f}s")
print(f"  HYP-270 OK by β₃: {dict(sorted(ok_b3.items()))}")
print(f"  HYP-270 FAIL by β₃: {dict(sorted(fail_b3.items()))}")

total_ok = sum(ok_b3.values())
total_fail = sum(fail_b3.values())
print(f"  Total OK: {total_ok}, Total FAIL: {total_fail}")

if 1 in fail_b3:
    pct = 100*fail_b3[1]/(ok_b3.get(1,0)+fail_b3.get(1,0)) if (ok_b3.get(1,0)+fail_b3.get(1,0)) > 0 else 0
    print(f"  β₃=1: {pct:.1f}% have HYP-270 failure")
if 0 in fail_b3:
    pct = 100*fail_b3[0]/(ok_b3.get(0,0)+fail_b3.get(0,0)) if (ok_b3.get(0,0)+fail_b3.get(0,0)) > 0 else 0
    print(f"  β₃=0: {pct:.1f}% have HYP-270 failure")

print("\nDone.")
