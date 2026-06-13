#!/usr/bin/env python3
"""
beta2_dt_coverage_n6.py — Verify DT boundaries span Z₂ at n=6

At n=5, DT path boundaries span Z₂ in ALL 1024 tournaments.
This extends the verification to n=6 (32768 tournaments).

If true, this confirms the "tetrahedral decomposition" proof approach:
β₂=0 because DT tetrahedra tile all 2-cycles.

Also analyzes the STRUCTURE of the coverage more carefully:
- How many DT paths are needed?
- What's the minimal covering set?
- Does the rank of DT boundaries always equal dim(Z₂)?

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
import time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def check_dt_coverage(A, n):
    """Check if DT path boundaries span Z₂ for tournament A.

    Returns (z2_dim, dt_rank, num_dt, sufficient).
    """
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    if not ap2:
        return (0, 0, 0, True)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    if d2 == 0:
        return (0, 0, 0, True)

    # Compute Z₂ = ker(∂₂: Ω₂ → Ω₁)
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2  # Image of Ω₂ basis in A₁ coords

    # Rank of ∂₂ restricted to Ω₂ (in Ω₁ coords)
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2 = d2 - rk2

    if z2 == 0:
        return (0, 0, 0, True)

    if not ap3:
        return (z2, 0, 0, z2 == 0)

    # Identify DT paths: (a,b,c,d) with a→c AND b→d
    ap3_tuples = [tuple(p) for p in ap3]
    dt_indices = [i for i, (a,b,c,d) in enumerate(ap3_tuples)
                  if A[a][c] and A[b][d]]
    num_dt = len(dt_indices)

    if num_dt == 0:
        return (z2, 0, 0, z2 == 0)

    # DT boundaries in A₂ space
    bd3 = build_full_boundary_matrix(ap3, ap2)
    dt_boundaries = bd3[:, dt_indices]  # |A₂| × num_dt

    # Project into Ω₂ coords
    dt_in_om2 = np.linalg.lstsq(om2, dt_boundaries, rcond=None)[0]
    rk_dt = np.linalg.matrix_rank(dt_in_om2, tol=1e-8)

    return (z2, rk_dt, num_dt, rk_dt >= z2)


def check_all_om3_coverage(A, n):
    """Check if ALL Ω₃ boundaries (not just DT) span Z₂."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    if not ap2:
        return (0, 0, True)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0, 0))

    d2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    d3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0
    if d2 == 0:
        return (0, 0, True)

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2 = d2 - rk2

    if z2 == 0 or d3 == 0:
        return (z2, 0, z2 == 0)

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    bd3_in_om2 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
    b2 = np.linalg.matrix_rank(bd3_in_om2, tol=1e-8)

    return (z2, b2, z2 == b2)  # β₂ = z2 - b2 = 0 iff z2 == b2


# ===== MAIN =====
print("=" * 70)
print("DT COVERAGE OF Z₂ AT n=6")
print("=" * 70)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)
total_T = 1 << m
print(f"n={n}, {total_T} tournaments, {m} edges")

dt_sufficient = Counter()
z2_dist = Counter()
dt_excess = Counter()  # rk_dt - z2 when sufficient
failures = []
t0 = time.time()

for bits in range(total_T):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    z2, rk_dt, num_dt, suff = check_dt_coverage(A, n)
    z2_dist[z2] += 1
    dt_sufficient[suff] += 1

    if suff and z2 > 0:
        dt_excess[rk_dt - z2] += 1

    if not suff:
        scores = sorted(sum(A[i]) for i in range(n))
        failures.append((bits, scores, z2, rk_dt, num_dt))
        if len(failures) <= 5:
            print(f"  FAIL: T#{bits} scores={scores}, Z₂={z2}, DT_rank={rk_dt}, #DT={num_dt}")

    if bits % 2000 == 0 and bits > 0:
        elapsed = time.time() - t0
        rate = bits / elapsed
        eta = (total_T - bits) / rate
        print(f"  ... {bits}/{total_T} ({elapsed:.0f}s, ETA {eta:.0f}s), "
              f"fails={len(failures)}")

elapsed = time.time() - t0
print(f"\nCompleted in {elapsed:.1f}s")

print(f"\nZ₂ dimension distribution:")
for val, count in sorted(z2_dist.items()):
    print(f"  Z₂ = {val}: {count}")

print(f"\nDT coverage:")
print(f"  Sufficient: {dt_sufficient[True]}/{total_T}")
print(f"  Insufficient: {dt_sufficient.get(False, 0)}/{total_T}")

if dt_excess:
    print(f"\nDT rank excess over Z₂ (when sufficient and Z₂>0):")
    for val, count in sorted(dt_excess.items()):
        print(f"  excess {val}: {count}")

if failures:
    print(f"\n{len(failures)} failures — analyzing:")
    for bits, scores, z2, rk_dt, num_dt in failures[:10]:
        print(f"  T#{bits} scores={scores}: Z₂={z2}, DT_rank={rk_dt}, #DT={num_dt}")

    # For first failure, check if ALL Ω₃ boundaries cover Z₂
    print("\n  Checking if full Ω₃ covers Z₂ in failure cases...")
    for bits, scores, z2, rk_dt, num_dt in failures[:5]:
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        z2_full, b2_full, full_ok = check_all_om3_coverage(A, n)
        print(f"    T#{bits}: Z₂={z2_full}, B₂={b2_full}, β₂={z2_full-b2_full}, "
              f"full Ω₃ covers Z₂: {full_ok}")
else:
    print("\nDT BOUNDARIES SPAN Z₂ IN ALL TOURNAMENTS! ✓")
    print("This confirms the tetrahedral decomposition approach.")

print("\nDone.")
