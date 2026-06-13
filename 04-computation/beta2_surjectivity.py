#!/usr/bin/env python3
"""
beta2_surjectivity.py - Analyze WHY d3: Om3 → Z_2 is surjective

For beta2 = 0, we need: every 2-cycle z ∈ Z_2 has a preimage in Om3.

Strategy: For each Z_2 basis element, find which Om3 elements map to it.
Analyze the STRUCTURE of the surjection, not just its existence.

Key questions:
1. Is every Z_2 element a boundary of DD paths alone? (DT sufficiency)
2. What is the "minimal" Om3 covering Z_2?
3. Is there a combinatorial/structural explanation?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
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


def analyze_surjectivity(A, n):
    """Decompose the surjection d3: Om3 → Z_2."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    d2 = om2.shape[1] if om2.ndim == 2 else 0
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    if d2 == 0:
        return {'dim_om2': 0, 'dim_z2': 0, 'dim_om3': 0, 'beta2': 0}

    # Build d2: A_2 → A_1 and d3: A_3 → A_2
    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    # d2 restricted to Om2
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk_d2 = sum(s > 1e-8 for s in S2)
    dim_z2 = d2 - rk_d2  # dim ker(d2|_{Om2})

    if dim_z2 == 0:
        return {'dim_om2': d2, 'dim_z2': 0, 'dim_om3': d3, 'beta2': 0}

    # d3 restricted to Om3 (image lands in A_2)
    if d3 > 0:
        bd3_om = bd3 @ om3  # |A_2| x d3 matrix

        # Rank = dim(im(d3))
        S3 = np.linalg.svd(bd3_om, compute_uv=False)
        rk_d3 = sum(s > 1e-8 for s in S3)
    else:
        bd3_om = np.zeros((len(a2), 0))
        rk_d3 = 0

    beta2 = dim_z2 - rk_d3

    # Now decompose: which part of Z_2 is filled by DD vs non-DD?
    # Classify 3-paths
    dd_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

    if dd_idx and d3 > 0:
        # DD part of Om3: project Om3 onto DD paths
        # Each DD path is individually in Om3, so we can take them directly
        dd_in_a3 = set(dd_idx)

        # Build d3 restricted to DD paths only
        bd3_dd = bd3[:, dd_idx]  # boundary of DD paths
        S_dd = np.linalg.svd(bd3_dd, compute_uv=False)
        rk_dd = sum(s > 1e-8 for s in S_dd)
    else:
        rk_dd = 0

    return {
        'dim_om2': d2,
        'dim_z2': dim_z2,
        'dim_om3': d3,
        'rk_d3': rk_d3,
        'beta2': beta2,
        'n_dd': len(dd_idx),
        'rk_dd_boundary': rk_dd,
        'dd_fills_z2': rk_dd >= dim_z2,
    }


# ============================================================
# TEST 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("SURJECTIVITY ANALYSIS: n=5")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

stats = Counter()
dd_fill_stats = Counter()
beta2_nonzero = 0

for bits in range(total):
    A = build_adj(n, bits)
    r = analyze_surjectivity(A, n)

    if r['beta2'] != 0:
        beta2_nonzero += 1
        print(f"  BETA2 != 0: bits={bits}, {r}")

    if r['dim_z2'] > 0:
        stats[r['dd_fills_z2']] += 1
        dd_fill_stats[(r['dim_z2'], r['rk_dd_boundary'])] += 1

print(f"\nn=5: {beta2_nonzero}/{total} with beta2!=0")
print(f"DD-fills-Z_2: {dict(stats)}")
print(f"(dim_Z_2, rk_DD_boundary) distribution:")
for k in sorted(dd_fill_stats.keys()):
    print(f"  {k}: {dd_fill_stats[k]}")


# ============================================================
# TEST 2: n=6 exhaustive
# ============================================================
print(f"\n{'='*70}")
print("SURJECTIVITY ANALYSIS: n=6")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

stats6 = Counter()
dd_fill6 = Counter()
beta2_nonzero6 = 0
deficit_examples = []

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), beta2!=0: {beta2_nonzero6}")

    A = build_adj(n, bits)
    r = analyze_surjectivity(A, n)

    if r['beta2'] != 0:
        beta2_nonzero6 += 1
        if beta2_nonzero6 <= 3:
            print(f"  BETA2 != 0: bits={bits}, {r}")

    if r['dim_z2'] > 0:
        stats6[r['dd_fills_z2']] += 1
        dd_fill6[(r['dim_z2'], r['rk_dd_boundary'])] += 1
        if not r['dd_fills_z2']:
            deficit = r['dim_z2'] - r['rk_dd_boundary']
            if len(deficit_examples) < 10:
                deficit_examples.append((bits, r))

dt = time.time() - t0
print(f"\nn=6: Done in {dt:.0f}s")
print(f"beta2!=0: {beta2_nonzero6}/{total}")
print(f"DD-fills-Z_2: {dict(stats6)}")
print(f"\n(dim_Z_2, rk_DD_boundary) distribution:")
for k in sorted(dd_fill6.keys()):
    print(f"  {k}: {dd_fill6[k]}")

if deficit_examples:
    print(f"\nDD-deficit examples (DD doesn't fill Z_2):")
    for bits, r in deficit_examples[:5]:
        A = build_adj(6, bits)
        scores = tuple(sorted([sum(row) for row in A]))
        deficit = r['dim_z2'] - r['rk_dd_boundary']
        print(f"  bits={bits}, scores={scores}, Z_2={r['dim_z2']}, DD_rk={r['rk_dd_boundary']}, "
              f"deficit={deficit}, Om3={r['dim_om3']}, d3_rk={r['rk_d3']}, beta2={r['beta2']}")


# ============================================================
# TEST 3: Boundary structure - what do non-DD elements contribute?
# ============================================================
print(f"\n{'='*70}")
print("NON-DD CONTRIBUTION TO FILLING")
print("=" * 70)

# At n=6, for DD-deficit cases, what part of Z_2 do non-DD Om3 elements fill?
n = 6
for bits, r in deficit_examples[:3]:
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    dd_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    nondd_idx = [i for i, p in enumerate(a3) if not (A[p[0]][p[2]] and A[p[1]][p[3]])]

    # Get Om3 basis
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    # Analyze which Om3 elements are "pure DD"
    n_pure_dd = 0
    n_mixed = 0
    for col in range(d3):
        v = om3[:, col]
        has_dd = any(abs(v[i]) > 1e-8 for i in dd_idx) if dd_idx else False
        has_nondd = any(abs(v[i]) > 1e-8 for i in nondd_idx) if nondd_idx else False
        if has_dd and not has_nondd:
            n_pure_dd += 1
        elif has_nondd:
            n_mixed += 1

    scores = tuple(sorted([sum(row) for row in A]))
    print(f"\n  bits={bits}, scores={scores}")
    print(f"  |DD|={len(dd_idx)}, |non-DD|={len(nondd_idx)}")
    print(f"  dim(Om3)={d3}: pure_DD={n_pure_dd}, mixed/pure_nonDD={n_mixed}")

    # Classify non-DD types
    type_02 = [i for i in nondd_idx if not A[a3[i][0]][a3[i][2]] and A[a3[i][1]][a3[i][3]]]
    type_13 = [i for i in nondd_idx if A[a3[i][0]][a3[i][2]] and not A[a3[i][1]][a3[i][3]]]
    type_nn = [i for i in nondd_idx if not A[a3[i][0]][a3[i][2]] and not A[a3[i][1]][a3[i][3]]]
    print(f"  non-DD breakdown: 02={len(type_02)}, 13={len(type_13)}, NN={len(type_nn)}")


# ============================================================
# TEST 4: n=7 sampled
# ============================================================
print(f"\n{'='*70}")
print("SURJECTIVITY ANALYSIS: n=7 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 7
n_samples = 5000

beta2_nonzero7 = 0
dd_fill_fail7 = 0
max_z2 = 0

t0 = time.time()
for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    r = analyze_surjectivity(A, n)

    if r['beta2'] != 0:
        beta2_nonzero7 += 1

    if r['dim_z2'] > 0:
        if not r.get('dd_fills_z2', True):
            dd_fill_fail7 += 1
        max_z2 = max(max_z2, r['dim_z2'])

    if trial % 1000 == 999:
        dt = time.time() - t0
        print(f"  ... {trial+1}/{n_samples} ({dt:.0f}s), beta2!=0: {beta2_nonzero7}, DD_fail: {dd_fill_fail7}")

dt = time.time() - t0
print(f"\nn=7: {n_samples} samples in {dt:.0f}s")
print(f"beta2!=0: {beta2_nonzero7}/{n_samples}")
print(f"DD doesn't fill Z_2: {dd_fill_fail7}/{n_samples}")
print(f"Max dim(Z_2): {max_z2}")


# ============================================================
# TEST 5: n=8 sampled
# ============================================================
print(f"\n{'='*70}")
print("SURJECTIVITY ANALYSIS: n=8 (sampled)")
print("=" * 70)

n = 8
n_samples = 1000

beta2_nonzero8 = 0
dd_fill_fail8 = 0

t0 = time.time()
for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    r = analyze_surjectivity(A, n)

    if r['beta2'] != 0:
        beta2_nonzero8 += 1

    if r['dim_z2'] > 0 and not r.get('dd_fills_z2', True):
        dd_fill_fail8 += 1

dt = time.time() - t0
print(f"\nn=8: {n_samples} samples in {dt:.0f}s")
print(f"beta2!=0: {beta2_nonzero8}/{n_samples}")
print(f"DD doesn't fill Z_2: {dd_fill_fail8}/{n_samples}")


print("\n\nDone.")
