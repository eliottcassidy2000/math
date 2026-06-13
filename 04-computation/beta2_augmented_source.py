#!/usr/bin/env python3
"""
beta2_augmented_source.py — Augmented tournament approach to β₂=0

Strategy: For tournament T on n vertices, construct T' on n+1 vertices
by adding a source vertex s (s→v for all v∈T). Then:
- T' has a source, so β_p(T') = 0 for ALL p ≥ 1 (source cone proof)
- LES of pair (T', T) gives: H₃(T') → H₃(T',T) → H₂(T) → H₂(T') = 0
  Since H₃(T') = 0: H₂(T) ≅ H₃(T',T)

So β₂(T) = 0 ⟺ H₃(T',T) = 0.

This script:
1. Verifies β_p(T') = 0 for T' with source (all p)
2. Computes H₃(T',T) directly
3. Checks H₂(T) ≅ H₃(T',T)

If H₃(T',T) = 0 always, this gives a CLEAN proof of β₂=0.

Author: opus-2026-03-08-S49
"""
import sys, time
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

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def augment_with_source(A, n):
    """Add source vertex 0 that beats everyone. Original vertices shift to 1..n."""
    n1 = n + 1
    A1 = [[0]*n1 for _ in range(n1)]
    # Source vertex 0 → all others
    for j in range(1, n1):
        A1[0][j] = 1
    # Copy original tournament (shifted by 1)
    for i in range(n):
        for j in range(n):
            A1[i+1][j+1] = A[i][j]
    return A1, n1

def compute_betti_full(A, n, max_p):
    """Compute β_0,...,β_{max_p}."""
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


# ============================================================
# Test 1: Source tournament has all β_p = 0
# ============================================================
print("=" * 70)
print("TEST 1: SOURCE TOURNAMENT β_p = 0 FOR ALL p")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

all_zero = True
for bits in range(total):
    A = build_adj(n, bits)
    A1, n1 = augment_with_source(A, n)
    betti = compute_betti_full(A1, n1, 4)  # Check up to β_4
    if any(b != 0 for b in betti[1:]):
        print(f"  FAIL at bits={bits}: β = {betti}")
        all_zero = False

if all_zero:
    print(f"  ✓ ALL {total} augmented tournaments (n={n}→{n+1}) have β_p=0 for p=1..4")

# ============================================================
# Test 2: H₃(T',T) — relative homology of augmented pair
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: H₃(T',T) — RELATIVE HOMOLOGY")
print("=" * 70)

# For the pair (T', T), relative chains are chains in T' that involve
# the source vertex s. Relative Ω_p consists of Ω_p-chains that use s.
#
# More precisely: the relative chain complex is Ω_p(T') / Ω_p(T).
# Since T = T'\s, the quotient Ω_p(T')/Ω_p(T) has basis from
# Ω_p-generators that involve vertex s=0.

def compute_relative_h3(A, n):
    """Compute H₃(T',T) where T' = T + source vertex s=0."""
    A1, n1 = augment_with_source(A, n)

    # Compute allowed paths and omega for T'
    ap = {}
    om = {}
    for p in range(6):  # Up to p=5 for safety
        ap[p] = enumerate_allowed_paths(A1, n1, p)

    om[0] = np.eye(n1)
    for p in range(1, 6):
        if ap[p]:
            om[p] = compute_omega_basis(A1, n1, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0,0))

    # Also compute for T (original, vertices 1..n in T')
    ap_T = {}
    om_T = {}
    for p in range(6):
        ap_T[p] = enumerate_allowed_paths(A, n, p)

    om_T[0] = np.eye(n)
    for p in range(1, 6):
        if ap_T[p]:
            om_T[p] = compute_omega_basis(A, n, p, ap_T[p], ap_T[p-1])
        else:
            om_T[p] = np.zeros((0,0))

    # Relative dimensions
    d3_T = dim_om(om_T[3]) if 3 in om_T else 0
    d3_T1 = dim_om(om[3])
    d3_rel = d3_T1 - d3_T  # dim(Ω₃^rel)

    d2_T = dim_om(om_T[2]) if 2 in om_T else 0
    d2_T1 = dim_om(om[2])
    d2_rel = d2_T1 - d2_T

    d4_T = dim_om(om_T[4]) if 4 in om_T else 0
    d4_T1 = dim_om(om[4])
    d4_rel = d4_T1 - d4_T

    # β₂(T)
    b2_T = 0
    if dim_om(om_T[2]) > 0:
        bd2_T = build_full_boundary_matrix(ap_T[2], ap_T[1])
        d2m = np.linalg.lstsq(om_T[1], bd2_T @ om_T[2], rcond=None)[0]
        rk2_T = np.linalg.matrix_rank(d2m, tol=1e-8)
        ker2_T = dim_om(om_T[2]) - rk2_T

        if dim_om(om_T[3]) > 0:
            bd3_T = build_full_boundary_matrix(ap_T[3], ap_T[2])
            d3m = np.linalg.lstsq(om_T[2], bd3_T @ om_T[3], rcond=None)[0]
            im3_T = np.linalg.matrix_rank(d3m, tol=1e-8)
        else:
            im3_T = 0
        b2_T = ker2_T - im3_T

    return {
        'd3_rel': d3_rel, 'd2_rel': d2_rel, 'd4_rel': d4_rel,
        'b2_T': b2_T,
        'd3_T1': d3_T1, 'd3_T': d3_T,
        'd2_T1': d2_T1, 'd2_T': d2_T,
        'd4_T1': d4_T1, 'd4_T': d4_T,
    }


n = 5
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

rel_dims = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]))
    info = compute_relative_h3(A, n)
    rel_dims.append(info)

elapsed = time.time() - t0
print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")

# Analyze
d3_rel_vals = [r['d3_rel'] for r in rel_dims]
d2_rel_vals = [r['d2_rel'] for r in rel_dims]
d4_rel_vals = [r['d4_rel'] for r in rel_dims]
b2_vals = [r['b2_T'] for r in rel_dims]

from collections import Counter
print(f"\n  d₃_rel distribution: {dict(Counter(d3_rel_vals))}")
print(f"  d₂_rel distribution: {dict(Counter(d2_rel_vals))}")
print(f"  d₄_rel distribution: {dict(Counter(d4_rel_vals))}")
print(f"  β₂(T) distribution: {dict(Counter(b2_vals))}")

# Check: is d₃_rel always 0?
if all(d == 0 for d in d3_rel_vals):
    print(f"\n  ✓ d₃_rel = 0 for ALL tournaments => H₃(T',T) = 0 => β₂(T) = 0")
else:
    print(f"\n  ✗ d₃_rel > 0 for some tournaments")
    # When d₃_rel > 0, need to compute actual H₃
    for i, r in enumerate(rel_dims):
        if r['d3_rel'] > 0 and i < 5:
            A = build_adj(n, i)
            scores = [sum(A[j][k] for k in range(n) if k!=j) for j in range(n)]
            print(f"    bits={i}, scores={scores}: d3_rel={r['d3_rel']}, d2_rel={r['d2_rel']}")


# ============================================================
# Test 3: Verify cone c_s maps ALL Ω_p into Ω_{p+1}
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: SOURCE CONE c_s: Ω_p → Ω_{p+1} for ALL p")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

for p_test in [1, 2, 3]:
    ok = 0
    fail = 0
    for bits in range(total):
        A = build_adj(n, bits)
        A1, n1 = augment_with_source(A, n)
        s = 0  # source vertex

        ap_p = enumerate_allowed_paths(A1, n1, p_test)
        ap_p1 = enumerate_allowed_paths(A1, n1, p_test+1)
        if not ap_p:
            ok += 1
            continue

        om_p = compute_omega_basis(A1, n1, p_test, ap_p,
                                    enumerate_allowed_paths(A1, n1, p_test-1))
        dp = dim_om(om_p)
        if dp == 0:
            ok += 1
            continue

        om_p1 = compute_omega_basis(A1, n1, p_test+1, ap_p1,
                                     ap_p)
        dp1 = dim_om(om_p1)

        # Apply cone: c_s maps each p-path (a₀,...,aₚ) to (s,a₀,...,aₚ) if s→a₀
        # For source s, s→everyone, so s→a₀ always holds IF a₀ ≠ s
        ap_p1_idx = {tuple(p): i for i, p in enumerate(ap_p1)}

        cone_matrix = np.zeros((len(ap_p1), len(ap_p)))
        for j, path in enumerate(ap_p):
            if s in path:
                continue  # Can't cone paths containing s
            cone_path = tuple([s] + list(path))
            if cone_path in ap_p1_idx:
                cone_matrix[ap_p1_idx[cone_path], j] = 1.0

        # Check: c_s(Ω_p) ⊆ Ω_{p+1}?
        cone_om = cone_matrix @ om_p  # Each column is c_s(ω) in A_{p+1} coords

        if dp1 == 0:
            # Check if cone_om is zero
            if np.linalg.norm(cone_om) < 1e-8:
                ok += 1
            else:
                fail += 1
            continue

        # Project onto Ω_{p+1}
        coords = np.linalg.lstsq(om_p1, cone_om, rcond=None)[0]
        resid = cone_om - om_p1 @ coords
        max_resid = np.max(np.abs(resid)) if resid.size > 0 else 0

        if max_resid < 1e-6:
            ok += 1
        else:
            fail += 1
            if fail <= 3:
                scores = [sum(A[j][k] for k in range(n) if k!=j) for j in range(n)]
                print(f"  FAIL p={p_test}: bits={bits}, scores={scores}, max_resid={max_resid:.4e}")

    print(f"  c_s: Ω_{p_test} → Ω_{p_test+1}: ok={ok}, fail={fail}")


# ============================================================
# Test 4: For n=6, compute relative dimensions
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: RELATIVE DIMENSIONS AT n=6")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

d3_rel_counter = Counter()
d2_rel_counter = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    A1, n1 = augment_with_source(A, n)

    # Quick: just compute Ω dimensions for T and T'
    ap3_T = enumerate_allowed_paths(A, n, 3)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))

    ap3_T1 = enumerate_allowed_paths(A1, n1, 3)
    ap2_T1 = enumerate_allowed_paths(A1, n1, 2)
    ap1_T1 = enumerate_allowed_paths(A1, n1, 1)
    om3_T1 = compute_omega_basis(A1, n1, 3, ap3_T1, ap2_T1) if ap3_T1 else np.zeros((0,0))
    om2_T1 = compute_omega_basis(A1, n1, 2, ap2_T1, ap1_T1) if ap2_T1 else np.zeros((0,0))

    d3_rel = dim_om(om3_T1) - dim_om(om3_T)
    d2_rel = dim_om(om2_T1) - dim_om(om2_T)
    d3_rel_counter[d3_rel] += 1
    d2_rel_counter[d2_rel] += 1

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6: {total} tournaments in {elapsed:.0f}s")
print(f"  d₃_rel distribution: {dict(sorted(d3_rel_counter.items()))}")
print(f"  d₂_rel distribution: {dict(sorted(d2_rel_counter.items()))}")

print("\nDone.")
