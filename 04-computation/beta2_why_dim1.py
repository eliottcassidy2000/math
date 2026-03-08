#!/usr/bin/env python3
"""
beta2_why_dim1.py — WHY is h₂(T,T\\v) always 0 or 1?

Understanding the structural reason for dim H₂(T,T\\v) ≤ 1.

The relative chain complex: Ω_p(T)/Ω_p(T\\v).
We need to understand:
1. What determines dim Ω₂(T)/Ω₂(T\\v)?
2. What is the relative boundary map?
3. Why is the kernel of ∂₂^rel mod im(∂₃^rel) at most 1-dimensional?

Key quantities:
- dim Ω₂(T) = |TT triples| - |non-trivial Ω₂ conditions|
- dim Ω₂(T\\v) = same for T\\v
- d_rel = dim Ω₂(T) - dim Ω₂(T\\v) (not exactly, need to account for embedding)

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
print("WHY IS dim H₂(T,T\\v) ≤ 1?")
print("=" * 70)

# Part 1: Collect all dimensions of the relative chain complex
n = 5
m = n*(n-1)//2

print(f"\n--- n={n}: Complete dimension analysis ---\n")

# For each interior v with h2_rel > 0, collect all relevant dimensions
dims_data = []

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap0_T = enumerate_allowed_paths(A, n, 0)
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap3_T = enumerate_allowed_paths(A, n, 3)

        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
        om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
        om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))

        d = {
            'bits': bits, 'v': v, 'dv': dv,
            'd1_T': dim_om(om1_T), 'd2_T': dim_om(om2_T), 'd3_T': dim_om(om3_T),
            'd1_sub': dim_om(om1_sub), 'd2_sub': dim_om(om2_sub), 'd3_sub': dim_om(om3_sub),
            'A2_T': len(ap2_T), 'A2_sub': len(ap2_sub) if ap2_sub else 0,
            'A3_T': len(ap3_T) if ap3_T else 0, 'A3_sub': len(ap3_sub) if ap3_sub else 0,
            'beta1_sub': 0,
        }

        # Relative dimensions
        ap2_T_list = [tuple(p) for p in ap2_T]
        if d['d2_T'] == 0:
            continue

        if ap2_sub and d['d2_sub'] > 0:
            embed = np.zeros((len(ap2_T_list), d['d2_sub']))
            for j in range(d['d2_sub']):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d['d2_T'], 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        d['rk_embed2'] = rk_phi
        d['d_rel2'] = d['d2_T'] - rk_phi

        if d['d_rel2'] == 0:
            continue

        # Relative ∂₂
        coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]
        ap1_T_list = [tuple(p) for p in ap1_T]
        if d['d1_sub'] > 0:
            embed1 = np.zeros((len(ap1_T_list), d['d1_sub']))
            for j in range(d['d1_sub']):
                for k, path_sub in enumerate(ap1_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap1_T_list:
                        embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
            psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
        else:
            psi = np.zeros((d['d1_T'], 0))

        rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
        d['rk_embed1'] = rk_psi
        d['d_rel1'] = d['d1_T'] - rk_psi

        if rk_psi > 0:
            U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
            R = U_psi[:, rk_psi:]
        else:
            R = np.eye(d['d1_T'])

        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True) if rk_phi > 0 else (np.eye(d['d2_T']), None, None)
        Q = U_phi[:, rk_phi:] if rk_phi > 0 else np.eye(d['d2_T'])

        coords2_rel = R.T @ coords2_T @ Q
        rk_d2_rel = np.linalg.matrix_rank(coords2_rel, tol=1e-8)
        d['rk_d2_rel'] = rk_d2_rel
        d['z2_rel'] = d['d_rel2'] - rk_d2_rel

        if d['z2_rel'] == 0:
            continue

        # Relative ∂₃
        if d['d3_T'] > 0:
            coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
            d3_proj = Q.T @ coords3_T
            if d['d3_sub'] > 0:
                ap3_T_list = [tuple(p) for p in ap3_T]
                embed3 = np.zeros((len(ap3_T_list), d['d3_sub']))
                for j in range(d['d3_sub']):
                    for k, path_sub in enumerate(ap3_sub):
                        path_T = tuple(remap[x] for x in path_sub)
                        if path_T in ap3_T_list:
                            embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
                chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
                rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
                if rk_chi > 0:
                    U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                    d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
                else:
                    d3_proj_rel = d3_proj
                d['rk_embed3'] = rk_chi
            else:
                d3_proj_rel = d3_proj
                d['rk_embed3'] = 0
        else:
            d3_proj_rel = np.zeros((d['d_rel2'], 0))
            d['rk_embed3'] = 0

        d['d_rel3'] = d['d3_T'] - d['rk_embed3']
        rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
        d['rk_d3_rel'] = rk_d3_rel
        d['b2_rel'] = rk_d3_rel
        d['h2_rel'] = d['z2_rel'] - rk_d3_rel

        # β₁(T\v)
        bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
        rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
        z1 = d['d1_sub'] - rk_d1
        if d['d2_sub'] > 0:
            bd2_om = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
            b1 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
        else:
            b1 = 0
        d['beta1_sub'] = z1 - b1

        dims_data.append(d)

# Analyze the collected dimensions
print(f"Total interior (T,v) pairs with h₂_rel > 0: {len(dims_data)}")

# Print dimension table
print(f"\n  {'d_rel2':>6} {'d_rel1':>6} {'d_rel3':>6} | {'rk_∂2':>5} {'z2_rel':>6} {'b2_rel':>6} {'h2_rel':>6} | {'β1':>3} | count")
print(f"  {'─'*6} {'─'*6} {'─'*6}   {'─'*5} {'─'*6} {'─'*6} {'─'*6}   {'─'*3}   {'─'*5}")

dim_combos = Counter()
for d in dims_data:
    key = (d['d_rel2'], d['d_rel1'], d.get('d_rel3', 0), d['rk_d2_rel'], d['z2_rel'], d['b2_rel'], d['h2_rel'], d['beta1_sub'])
    dim_combos[key] += 1

for key, cnt in sorted(dim_combos.items()):
    d_rel2, d_rel1, d_rel3, rk, z2, b2, h2, b1 = key
    print(f"  {d_rel2:>6} {d_rel1:>6} {d_rel3:>6} | {rk:>5} {z2:>6} {b2:>6} {h2:>6} | {b1:>3} | {cnt}")


# Part 2: Look at the LES more carefully
print(f"\n{'='*70}")
print("Part 2: Full LES analysis")
print("=" * 70)

# The LES: ... → H₂(T\v) → H₂(T) → H₂(T,T\v) →^δ H₁(T\v) → H₁(T) → ...
# We know β₂(T\v) = 0 (by induction for n-1=4).
# So: 0 → H₂(T) → H₂(T,T\v) →^δ H₁(T\v) → H₁(T) → H₁(T,T\v) → ...
#
# Exactness at H₂(T,T\v): ker(δ) = im(j*) = im(H₂(T) → H₂(T,T\v))
# Since β₂(T\v) = 0, j* is injective: H₂(T) ↪ H₂(T,T\v)
# So ker(δ) ≅ H₂(T).
# If δ is injective, ker(δ) = 0, so H₂(T) = 0.
#
# From the table: when h2_rel = 1, we need δ ≠ 0.
# Equivalently: the single generator of H₂(T,T\v) maps to a nonzero element of H₁(T\v).
#
# But H₁(T\v) is also 1-dimensional. So δ is a linear map R → R.
# If δ = 0, then H₂(T) = H₂(T,T\v) = R, contradicting β₂ = 0.
# So actually, the β₂=0 fact IMPLIES δ is injective!
# But we're trying to PROVE β₂=0 from δ-injectivity...
#
# The circular argument breaks if we can prove h2_rel ≤ 1 and
# find a structural reason why δ ≠ 0.

# Alternatively: by the LES, β₂(T) = dim ker(δ).
# If h2_rel = 0 for all interior v: β₂(T) = 0 directly from h2_rel = 0
# (because H₂(T) ↪ H₂(T,T\v) = 0).
# If h2_rel = 1 for some interior v: need δ ≠ 0.

# How many tournaments have h2_rel > 0 for ALL interior vertices?
# vs h2_rel = 0 for SOME interior vertex (which would give β₂=0 automatically)?

print("\nHow often is h2_rel = 0 for some interior v?")
tours_with_zero_v = 0
tours_with_all_nonzero_v = 0

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    interior_v = [v for v in range(n) if 1 <= scores[v] <= n-2]
    if not interior_v:
        continue

    has_zero = False
    for v in interior_v:
        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap0_T = enumerate_allowed_paths(A, n, 0)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)

        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
        d2_T = dim_om(om2_T)
        if d2_T == 0:
            has_zero = True
            break

        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        d2_sub = dim_om(om2_sub)

        ap2_T_list = [tuple(p) for p in ap2_T]
        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        d_rel = d2_T - rk_phi
        if d_rel == 0:
            has_zero = True
            break

    if has_zero:
        tours_with_zero_v += 1
    else:
        tours_with_all_nonzero_v += 1

print(f"  Tournaments with SOME interior v giving h₂_rel=0 (or d_rel=0): {tours_with_zero_v}")
print(f"  Tournaments where ALL interior v give d_rel > 0: {tours_with_all_nonzero_v}")
print(f"  (Total with interior vertices: {tours_with_zero_v + tours_with_all_nonzero_v})")

if tours_with_zero_v == tours_with_zero_v + tours_with_all_nonzero_v:
    print("  ✓ EVERY tournament has SOME interior v with d_rel=0!")
    print("  => H₂(T) ↪ H₂(T,T\\v) = 0, so β₂=0 WITHOUT needing δ!")
elif tours_with_all_nonzero_v > 0:
    print(f"  Need δ-injectivity for {tours_with_all_nonzero_v} tournaments")


print(f"\n{'='*70}")
print("Part 3: Check at n=6 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 6
m = n*(n-1)//2
total = 1 << m
samples = random.sample(range(total), 2000)

tours_zero_6 = 0
tours_nonzero_6 = 0
t0 = time.time()

for idx, bits in enumerate(samples):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    interior_v = [v for v in range(n) if 1 <= scores[v] <= n-2]

    if not interior_v:
        continue

    has_zero = False
    for v in interior_v[:3]:  # check up to 3
        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
        d2_T = dim_om(om2_T)
        if d2_T == 0:
            has_zero = True
            break

        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        d2_sub = dim_om(om2_sub)

        ap2_T_list = [tuple(p) for p in ap2_T]
        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        d_rel = d2_T - rk_phi
        if d_rel == 0:
            has_zero = True
            break

    if has_zero:
        tours_zero_6 += 1
    else:
        tours_nonzero_6 += 1

    if (idx + 1) % 500 == 0:
        print(f"  {idx+1}/2000 ({time.time()-t0:.0f}s)")

print(f"\nn=6 (2000 samples):")
print(f"  Has SOME interior v with d_rel=0: {tours_zero_6}")
print(f"  ALL interior v have d_rel > 0: {tours_nonzero_6}")

if tours_nonzero_6 == 0:
    print("  ✓ EVERY tournament has SOME interior v with d_rel=0!")
    print("  => β₂=0 WITHOUT needing δ!")
else:
    print(f"  Need δ for {tours_nonzero_6} tournaments")


print("\nDone.")
