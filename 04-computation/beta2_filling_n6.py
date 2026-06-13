#!/usr/bin/env python3
"""
beta2_filling_n6.py - Z_2 cycle structure at n=6 (sampled)

Check whether the patterns from n=5 extend:
1. Do ALL Z_2 cycles use all n vertices?
2. Is surplus = dim(Omega_3) - dim(Z_2) >= 0 always?
3. What is the distribution of dim(Z_2) and dim(Z_3)?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
from collections import Counter
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


def analyze_tournament(A, n):
    """Compute key quantities for beta_2 analysis."""
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

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(5)}

    rks = {}
    bd_om = {}
    for p in [2, 3]:
        if dims[p] > 0 and len(paths.get(p-1, [])) > 0 and dims[p-1] > 0:
            bd = build_full_boundary_matrix(paths[p], paths[p-1])
            raw = bd @ omega[p]
            coords, _, _, _ = np.linalg.lstsq(omega[p-1], raw, rcond=None)
            bd_om[p] = coords
            rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)
        else:
            bd_om[p] = None
            rks[p] = 0

    if dims[1] > 0:
        bd1 = build_full_boundary_matrix(paths[1], paths[0])
        rks[1] = np.linalg.matrix_rank(bd1 @ omega[1], tol=1e-8)
    else:
        rks[1] = 0

    dimZ2 = dims[2] - rks.get(2, 0)
    dimZ3 = dims[3] - rks.get(3, 0)
    beta1 = (dims[1] - rks.get(1, 0)) - rks.get(2, 0)
    beta2 = dimZ2 - rks.get(3, 0)
    surplus = dims[3] - dimZ2

    # Cycle support analysis (only for small n or few cycles)
    cycle_min_support = n  # default
    if dimZ2 > 0 and bd_om.get(2) is not None and n <= 6:
        U, S, Vt = np.linalg.svd(bd_om[2], full_matrices=True)
        rk2 = int(np.sum(np.abs(S) > 1e-8))
        Z2_coords = Vt[rk2:].T

        min_support = n
        for j in range(Z2_coords.shape[1]):
            z_paths = omega[2] @ Z2_coords[:, j]
            verts = set()
            for i, c in enumerate(z_paths):
                if abs(c) > 1e-8:
                    verts.update(paths[2][i])
            min_support = min(min_support, len(verts))
        cycle_min_support = min_support

    return {
        'dims': dims, 'rks': rks, 'dimZ2': dimZ2, 'dimZ3': dimZ3,
        'beta1': beta1, 'beta2': beta2, 'surplus': surplus,
        'min_cycle_support': cycle_min_support
    }


random.seed(42)
print("=" * 70)
print("Z_2 CYCLE STRUCTURE AT n=6 (sampled)")
print("=" * 70)

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs
samples = 5000

dimZ2_dist = Counter()
dimZ3_dist = Counter()
surplus_dist = Counter()
beta2_dist = Counter()
beta1_dist = Counter()
support_dist = Counter()

t0 = time.time()
for s in range(samples):
    if s % 1000 == 0 and s > 0:
        dt = time.time() - t0
        print(f"  {s}/{samples} ({dt:.0f}s)")

    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    data = analyze_tournament(A, n)

    dimZ2_dist[data['dimZ2']] += 1
    dimZ3_dist[data['dimZ3']] += 1
    surplus_dist[data['surplus']] += 1
    beta2_dist[data['beta2']] += 1
    beta1_dist[data['beta1']] += 1
    support_dist[data['min_cycle_support']] += 1

dt = time.time() - t0
print(f"  Done in {dt:.0f}s")

print(f"\ndim(Z_2) distribution: {dict(sorted(dimZ2_dist.items()))}")
print(f"dim(Z_3) distribution: {dict(sorted(dimZ3_dist.items()))}")
print(f"surplus distribution: {dict(sorted(surplus_dist.items()))}")
print(f"beta_2 distribution: {dict(sorted(beta2_dist.items()))}")
print(f"beta_1 distribution: {dict(sorted(beta1_dist.items()))}")
print(f"min cycle support sizes: {dict(sorted(support_dist.items()))}")

if all(v == 0 for v in beta2_dist.keys() if beta2_dist[v] > 0):
    print("\nbeta_2 = 0 for ALL sampled tournaments!")

if all(s >= 0 for s in surplus_dist.keys()):
    print("surplus >= 0 for ALL sampled tournaments!")

# Also check: does dim(Z_2) have a formula?
print(f"\n{'='*70}")
print("FORMULA CHECK: dim(Z_2) = dim(Om_2) - rk(d_2) = dim(Om_2) - C(n-1,2) + beta_1")
print("=" * 70)

z2_by_dimOm2_b1 = Counter()
for s in range(min(samples, 1000)):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    data = analyze_tournament(A, n)
    z2_by_dimOm2_b1[(data['dims'][2], data['beta1'], data['dimZ2'])] += 1

print("Checking dim(Z_2) = dim(Om_2) - C(n-1,2) + beta_1:")
Z1 = (n-1)*(n-2)//2
errors = 0
for key, count in z2_by_dimOm2_b1.items():
    dimOm2, b1, dimZ2 = key
    predicted = dimOm2 - Z1 + b1
    if predicted != dimZ2:
        errors += 1
        print(f"  ERROR: dim(Om2)={dimOm2}, b1={b1}, dim(Z2)={dimZ2}, predicted={predicted}")
print(f"  Errors: {errors}")

# Exhaustive check for specific score sequences at n=6
print(f"\n{'='*70}")
print("EXHAUSTIVE n=5 cross-check")
print("=" * 70)

n5 = 5
n5_arcs = n5*(n5-1)//2
total5 = 1 << n5_arcs

surplus5 = Counter()
support5 = Counter()

for bits in range(total5):
    A = build_adj(n5, bits)
    data = analyze_tournament(A, n5)
    surplus5[data['surplus']] += 1
    support5[data['min_cycle_support']] += 1

print(f"n=5 surplus: {dict(sorted(surplus5.items()))}")
print(f"n=5 cycle support: {dict(sorted(support5.items()))}")

print("\nDone.")
