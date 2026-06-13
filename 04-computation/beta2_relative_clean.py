#!/usr/bin/env python3
"""
beta2_relative_clean.py - Clean analysis of relative complex C_*(T, T\\v)

For the pair (T, T\\v), compute the relative chain complex and look for
dimensional formulas that explain H_2^rel = 0.

The key idea: if we can express dim(ker d_2^rel) and dim(im d_3^rel) in
terms of simple invariants (d_v, n, local structure), we might see WHY
they're always equal.

The relative complex has:
  Omega_p^rel(T, T\\v) = Omega_p(T) / Omega_p(T\\v)

For p=0: dim = 1 (just vertex v)
For p=1: dim = n-1 (all arcs involving v; for tournament, this is always n-1)

Strategy: use the beta2_relative_homology.py approach (known correct) to compute
H_2^rel, and also compute:
  - dim(Omega_2^rel)
  - dim(ker d_2^rel) = dim(Z_2^rel)
  - dim(im d_3^rel)
  - dim(Omega_3^rel)

as functions of the local structure at v.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
from collections import Counter, defaultdict
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

def random_tournament(n):
    return random.getrandbits(n*(n-1)//2)

def compute_relative_complex(A, n, v):
    """Compute the relative complex (T, T\\v) correctly.

    Returns dict with dims and ranks.
    Uses the CORRECT computation from beta2_relative_homology.py.
    """
    # Get paths and Omega bases for p=0,1,2,3
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

    result = {}

    # For each p, compute:
    # dim(Omega_p) = omega[p].shape[1]
    # dim(Omega_p(T\\v)) = dim(Omega_p) ∩ no-v = omega[p].shape[1] - rank(R_v)
    # dim(Omega_p^rel) = rank(R_v) where R_v = omega[p] restricted to v-paths
    for p in range(5):
        dim_om = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_om == 0:
            result[p] = {'dim_om': 0, 'dim_sub': 0, 'dim_rel': 0}
            continue

        v_idx = [i for i in range(len(paths[p])) if v in paths[p][i]]
        if not v_idx:
            result[p] = {'dim_om': dim_om, 'dim_sub': dim_om, 'dim_rel': 0}
            continue

        R_v = omega[p][v_idx, :]
        sv = np.linalg.svd(R_v, compute_uv=False)
        rank_rv = int(np.sum(np.abs(sv) > 1e-8))

        result[p] = {
            'dim_om': dim_om,
            'dim_sub': dim_om - rank_rv,
            'dim_rel': rank_rv,
        }

    # Now compute boundary ranks in the FULL complex
    for p in [1, 2, 3, 4]:
        dim_src = result[p]['dim_om']
        dim_tgt = result[p-1]['dim_om']
        if dim_src == 0 or dim_tgt == 0:
            result[p]['rk_bd'] = 0
            continue

        bd = build_full_boundary_matrix(paths[p], paths[p-1])
        bd_om = bd @ omega[p]  # |A_{p-1}| x dim_src

        # Express in Omega_{p-1} coordinates
        im_coords, _, _, _ = np.linalg.lstsq(omega[p-1], bd_om, rcond=None)
        result[p]['rk_bd'] = np.linalg.matrix_rank(im_coords, tol=1e-8)

    # Compute H_2^rel using the method from beta2_relative_homology.py
    # This is the CORRECT computation.

    dim_om2 = result[2]['dim_om']
    dim_om3 = result[3]['dim_om']

    if dim_om2 == 0:
        result['H2_rel'] = 0
        result['ker2_rel'] = 0
        result['im3_rel'] = 0
        return result

    # Build boundary d_2 in Omega coords
    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]

    # ker(d_2) in Omega_2 coords
    U2, S2, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk_d2 = int(np.sum(np.abs(S2) > 1e-8))
    ker2_dim = dim_om2 - rk_d2

    # V_p = Omega_p ∩ no-v (subspace of Omega_p)
    def get_V_basis(p):
        dim = result[p]['dim_om']
        if dim == 0: return np.zeros((dim, 0))
        v_idx = [i for i in range(len(paths[p])) if v in paths[p][i]]
        if not v_idx: return np.eye(dim)  # all of Omega_p is no-v
        R = omega[p][v_idx, :]
        U, S, Vt = np.linalg.svd(R, full_matrices=True)
        rk = int(np.sum(np.abs(S) > 1e-8))
        if rk >= dim: return np.zeros((dim, 0))
        return Vt[rk:].T

    V2 = get_V_basis(2)
    dim_V2 = V2.shape[1]

    # ker(d_2) ∩ V_2
    if ker2_dim == 0:
        result['H2_rel'] = 0
        result['ker2_rel'] = 0
        result['im3_rel'] = 0
        return result

    ker2_basis = Vt2[rk_d2:].T  # dim_om2 x ker2_dim

    if dim_V2 > 0:
        combined = np.hstack([ker2_basis, V2])
        rk_comb = np.linalg.matrix_rank(combined, tol=1e-8)
        inter_dim = ker2_dim + dim_V2 - rk_comb
    else:
        inter_dim = 0

    ker2_rel = ker2_dim - inter_dim
    result['ker2_rel'] = ker2_rel

    if ker2_rel == 0:
        result['H2_rel'] = 0
        result['im3_rel'] = 0
        return result

    # im(d_3^rel) in quotient Omega_2/V_2
    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        im3_coords, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)

        if dim_V2 > 0:
            combined_im = np.hstack([im3_coords, V2])
            rk_im_comb = np.linalg.matrix_rank(combined_im, tol=1e-8)
            im3_rel = rk_im_comb - dim_V2
        else:
            im3_rel = np.linalg.matrix_rank(im3_coords, tol=1e-8)
    else:
        im3_rel = 0

    result['im3_rel'] = im3_rel
    result['H2_rel'] = max(0, ker2_rel - im3_rel)

    return result


# ===== n=5 exhaustive =====
print("=" * 70)
print("CLEAN RELATIVE COMPLEX ANALYSIS")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print(f"\n--- n={n} exhaustive ({total} tournaments) ---")

# Group by (d_v, dim_rel_2, ker2_rel, im3_rel, H2_rel)
stats = Counter()
dim_rel_by_dv = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    for v_test in range(n):
        d_v = sum(A[v_test])
        res = compute_relative_complex(A, n, v_test)

        key = (d_v, res[2]['dim_rel'], res.get('ker2_rel', 0),
               res.get('im3_rel', 0), res.get('H2_rel', 0))
        stats[key] += 1

        # Track dim_rel_p as function of d_v
        dims = (res[1]['dim_rel'], res[2]['dim_rel'], res[3]['dim_rel'])
        dim_rel_by_dv[d_v].append(dims)

print(f"\n  (d_v, dim_O2_rel, ker_d2_rel, im_d3_rel, H2_rel): count")
for key in sorted(stats.keys()):
    d_v, dr2, k2r, i3r, h2r = key
    print(f"    d_v={d_v}: dim_O2_rel={dr2}, ker={k2r}, im={i3r}, H2_rel={h2r}, count={stats[key]}")

print(f"\n  FORMULAS for dim(Omega_p^rel) vs d_v at n={n}:")
for d_v in sorted(dim_rel_by_dv):
    all_dims = dim_rel_by_dv[d_v]
    unique_dims = Counter(all_dims)
    print(f"    d_v={d_v}: {dict(unique_dims)}")

# ===== n=6 exhaustive =====
n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print(f"\n--- n={n} exhaustive ({total} tournaments) ---")

stats_6 = Counter()
dim_rel_by_dv_6 = defaultdict(Counter)
h2_violations = 0

t0 = time.time()
for bits in range(total):
    if bits % 5000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), violations={h2_violations}")

    A = build_adj(n, bits)
    for v_test in range(n):
        d_v = sum(A[v_test])
        res = compute_relative_complex(A, n, v_test)

        h2r = res.get('H2_rel', 0)
        if h2r > 0:
            h2_violations += 1

        key = (d_v, res[2]['dim_rel'], res.get('ker2_rel', 0),
               res.get('im3_rel', 0), h2r)
        stats_6[key] += 1

        dims = (res[1]['dim_rel'], res[2]['dim_rel'], res[3]['dim_rel'])
        dim_rel_by_dv_6[d_v][dims] += 1

dt = time.time() - t0
print(f"\n  Completed in {dt:.0f}s. H2_rel violations: {h2_violations}")

print(f"\n  (d_v, dim_O2_rel, ker_d2_rel, im_d3_rel, H2_rel): count")
for key in sorted(stats_6.keys()):
    d_v, dr2, k2r, i3r, h2r = key
    marker = " ***" if h2r > 0 else ""
    print(f"    d_v={d_v}: dim_O2_rel={dr2}, ker={k2r}, im={i3r}, H2_rel={h2r}, count={stats_6[key]}{marker}")

print(f"\n  FORMULAS for dim(Omega_p^rel) vs d_v at n={n}:")
for d_v in sorted(dim_rel_by_dv_6):
    unique_dims = dim_rel_by_dv_6[d_v]
    total_count = sum(unique_dims.values())
    print(f"    d_v={d_v}: {len(unique_dims)} patterns, total={total_count}")
    for dims, cnt in sorted(unique_dims.items(), key=lambda x: -x[1])[:5]:
        print(f"      {dims}: {cnt} ({100*cnt/total_count:.1f}%)")

print("\nDone.")
