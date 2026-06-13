#!/usr/bin/env python3
"""
beta2_arcflip_n7_sample.py - Sample arc-flip preservation at n=7

Test HYP-220 (surplus >= |drop|) at n=7 via random sampling.
Also test whether the tight-case stability (HYP-221) extends to n=7.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
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

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_surplus(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dim_O2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    dim_O3 = omega_basis[3].shape[1] if omega_basis[3].ndim == 2 else 0

    if dim_O2 == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed[1] if 1 in allowed else [])
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    # im(d_3) = rk(d_3|Omega_3)
    if dim_O3 == 0 or 3 not in allowed or not allowed[3]:
        im_d3 = 0
    else:
        bd3 = build_full_boundary_matrix(allowed[3], allowed[2] if 2 in allowed else [])
        bd3_om = bd3 @ omega_basis[3]
        S_v = np.linalg.svd(bd3_om, compute_uv=False)
        im_d3 = int(np.sum(np.abs(S_v) > 1e-8))

    surplus = dim_O3 - Z2
    beta2 = Z2 - im_d3
    return surplus, beta2, dim_O2, dim_O3, Z2

def flip_arc(A, u, v, n):
    B = [row[:] for row in A]
    B[u][v] = 0
    B[v][u] = 1
    return B


n = 7
print("=" * 70)
print(f"ARC-FLIP PRESERVATION AT n={n} (SAMPLING)")
print("=" * 70)

N_SAMPLES = 500
N_FLIPS_PER = 5  # random flips per tournament
total_flips = 0
violations = 0
safe_count = 0
min_surplus = float('inf')
min_surplus_after = float('inf')
surplus_dist = {}

t0 = time.time()

for trial in range(N_SAMPLES):
    if trial % 100 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s, {total_flips} flips, {violations} violations)")

    A = random_tournament(n)
    surplus, beta2, _, _, _ = compute_surplus(A, n)

    if beta2 != 0:
        print(f"  WARNING: beta_2 = {beta2} for trial {trial}!")
        violations += 1
        continue

    min_surplus = min(min_surplus, surplus)
    surplus_dist[surplus] = surplus_dist.get(surplus, 0) + 1

    # Try random flips
    for _ in range(N_FLIPS_PER):
        # Pick a random arc to flip
        u = random.randint(0, n-1)
        v = random.randint(0, n-1)
        while u == v or A[u][v] == 0:
            u = random.randint(0, n-1)
            v = random.randint(0, n-1)

        B = flip_arc(A, u, v, n)
        surplus_flip, beta2_flip, _, _, _ = compute_surplus(B, n)
        total_flips += 1

        min_surplus_after = min(min_surplus_after, surplus_flip)

        if beta2_flip != 0:
            print(f"  VIOLATION! trial={trial}, flip ({u},{v}): surplus {surplus}->{surplus_flip}, beta2={beta2_flip}")
            violations += 1

        if surplus >= abs(surplus_flip - surplus) or surplus_flip >= 0:
            safe_count += 1

dt = time.time() - t0
print(f"\n  Samples: {N_SAMPLES}, total flips: {total_flips}")
print(f"  Violations (beta_2 != 0 after flip): {violations}")
print(f"  Min surplus (original): {min_surplus}")
print(f"  Min surplus (after flip): {min_surplus_after}")
print(f"  Safe (surplus >= |drop|): {safe_count}/{total_flips}")
print(f"  Time: {dt:.1f}s")
print(f"\n  Surplus distribution:")
for s in sorted(surplus_dist.keys()):
    print(f"    surplus={s}: {surplus_dist[s]} ({100*surplus_dist[s]/N_SAMPLES:.1f}%)")

# Additional: try to find tight tournaments (low surplus)
print(f"\n{'='*70}")
print(f"HUNTING FOR TIGHT TOURNAMENTS (low surplus)")
print(f"{'='*70}")

N_HUNT = 2000
min_found = float('inf')
tight_count = 0

for trial in range(N_HUNT):
    A = random_tournament(n)
    surplus, beta2, dim_O2, dim_O3, Z2 = compute_surplus(A, n)
    if surplus < min_found:
        min_found = surplus
        print(f"  New minimum surplus: {surplus} (O2={dim_O2}, O3={dim_O3}, Z2={Z2})")
    if surplus <= 3:
        tight_count += 1

print(f"\n  Minimum surplus found: {min_found}")
print(f"  Tight (surplus<=3): {tight_count}/{N_HUNT}")

print(f"\nDone. Total time: {time.time()-t0:.1f}s")
