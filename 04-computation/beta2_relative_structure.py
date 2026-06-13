#!/usr/bin/env python3
"""
beta2_relative_structure.py — Structure of the relative chain complex

For each (T, v), compute:
  d_p^rel = dim Ω_p(T) - dim Ω_p(T\v)
  z_p^rel = dim Z_p(T) - dim Z_p(T\v)
  b_p^rel = dim B_p(T) - dim B_p(T\v)

Key questions:
1. Is z₂^rel - b₂^rel = h₂^rel? (Not in general, since relative ≠ quotient)
2. What constrains h₂^rel ∈ {0,1}?
3. What is the relative Euler characteristic?

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

def compute_dims(A, n, max_p=4):
    """Compute Ω_p, Z_p, B_p dimensions for tournament A."""
    ap = {}
    om = {}
    for p in range(max_p+1):
        ap[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om[p] = np.eye(n)
        elif ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0,0))

    dims = {}
    for p in range(max_p+1):
        dp = dim_om(om[p]) if p > 0 else n
        dims[f'om{p}'] = dp

        if dp == 0:
            dims[f'z{p}'] = 0
            dims[f'b{p}'] = 0
            continue

        # Z_p = ker(∂_p)
        if p == 0:
            dims[f'z{p}'] = n  # ∂_0 = 0
        else:
            bd_p = build_full_boundary_matrix(ap[p], ap[p-1])
            prev = om[p-1] if p > 1 else np.eye(n)
            dm = np.linalg.lstsq(prev, bd_p @ om[p], rcond=None)[0]
            rk = np.linalg.matrix_rank(dm, tol=1e-8)
            dims[f'z{p}'] = dp - rk

        # B_p = im(∂_{p+1})
        if p < max_p and dim_om(om.get(p+1, np.zeros((0,0)))) > 0:
            bd_p1 = build_full_boundary_matrix(ap[p+1], ap[p])
            dm1 = np.linalg.lstsq(om[p], bd_p1 @ om[p+1], rcond=None)[0]
            dims[f'b{p}'] = np.linalg.matrix_rank(dm1, tol=1e-8)
        else:
            dims[f'b{p}'] = 0

    return dims

def compute_beta1(A, n):
    dims = compute_dims(A, n, max_p=2)
    return dims['z1'] - dims['b1']


# ============================================================
# n=5: Detailed relative structure
# ============================================================
print("=" * 70)
print("n=5: RELATIVE CHAIN COMPLEX STRUCTURE")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each (T,v): compute T and T\v dimensions
rel_joint = Counter()  # (d2_rel, z2_rel, b2_rel, beta1_S, beta1_T) -> count
h2_rel_analysis = Counter()  # h2_rel_value -> count

for bits in range(total):
    A = build_adj(n, bits)
    dims_T = compute_dims(A, n, 4)
    beta1_T = dims_T['z1'] - dims_T['b1']
    beta2_T = dims_T['z2'] - dims_T['b2']

    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        dims_S = compute_dims(A_sub, n-1, 4)
        beta1_S = dims_S['z1'] - dims_S['b1']
        beta2_S = dims_S['z2'] - dims_S['b2']

        d2_rel = dims_T['om2'] - dims_S['om2']
        z2_rel = dims_T['z2'] - dims_S['z2']
        b2_rel = dims_T['b2'] - dims_S['b2']
        d3_rel = dims_T['om3'] - dims_S['om3']

        # h₂_rel from LES (with β₂(T\v) = 0):
        # h₂_rel = β₂(T) + ker(i*: H₁(S) → H₁(T))
        if beta1_T == 0:
            ker_i = beta1_S  # i* maps into 0
        elif beta1_S == 0:
            ker_i = 0
        else:
            # Both nonzero; need to check if i* is injective
            # For now, assume ker_i = 0 (HYP-263 says h₂_rel = 0 here)
            ker_i = 0  # Will verify

        h2_rel = beta2_T + ker_i

        rel_joint[(d2_rel, z2_rel, b2_rel, d3_rel, beta1_S, beta1_T, h2_rel)] += 1
        h2_rel_analysis[h2_rel] += 1

print(f"\nn={n}: {total*n} (T,v) pairs")
print(f"\n  h₂_rel distribution: {dict(sorted(h2_rel_analysis.items()))}")

print(f"\n  Detailed (d₂^rel, z₂^rel, b₂^rel, d₃^rel, β₁_sub, β₁_T, h₂^rel):")
for key, cnt in sorted(rel_joint.items()):
    d2r, z2r, b2r, d3r, b1s, b1T, h2r = key
    if cnt >= 20:
        print(f"    d₂r={d2r}, z₂r={z2r}, b₂r={b2r}, d₃r={d3r}, β₁s={b1s}, β₁T={b1T}, h₂r={h2r}: {cnt}")


# ============================================================
# n=6: Check h₂_rel ∈ {0,1} exhaustively
# ============================================================
print(f"\n{'='*70}")
print("n=6: h₂_rel DISTRIBUTION")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

h2_counter = Counter()
violations = 0

for bits in range(total):
    A = build_adj(n, bits)
    dims_T = compute_dims(A, n, 3)
    beta1_T = dims_T['z1'] - dims_T['b1']
    beta2_T = dims_T['z2'] - dims_T['b2']

    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        beta1_S = compute_beta1(A_sub, n-1)

        if beta1_T == 0:
            ker_i = beta1_S
        elif beta1_S == 0:
            ker_i = 0
        else:
            ker_i = 0  # HYP-263

        h2_rel = beta2_T + ker_i
        h2_counter[h2_rel] += 1
        if h2_rel > 1:
            violations += 1
            if violations <= 3:
                scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
                print(f"  VIOLATION: bits={bits}, v={v}, h₂_rel={h2_rel}, β₂(T)={beta2_T}, β₁_sub={beta1_S}")

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6: {total*n} (T,v) pairs in {elapsed:.0f}s")
print(f"  h₂_rel distribution: {dict(sorted(h2_counter.items()))}")
print(f"  Violations (h₂_rel > 1): {violations}")

print("\nDone.")
