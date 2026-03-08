#!/usr/bin/env python3
"""Cone construction proof of beta_2 = 0 for tournaments.

KEY INSIGHT: For a source vertex v (v->w for all w!=v), the cone map
  s_v(sigma) = (v, sigma_0, ..., sigma_p)
gives a chain homotopy: d s + s d = id on the subcomplex of T\v.

Claim: s_v maps Omega_p elements (not using v) to Omega_{p+1}.
Proof: For u in Omega_p with no v, d(s_v(u)) = u - s_v(du).
  u in A_p. s_v(du): each term (v,tau) with tau in A_{p-1}, v->tau_0 (source).
  So d(s_v(u)) in A_p, meaning s_v(u) in Omega_{p+1}.

Then if z in ker d_p, d(s_v(z)) = z. So z is exact.

BUT: general tournaments DON'T have source vertices!

TEST:
1. Does cone from source work? (verified at n=5)
2. Do 2-cycles ever use ALL n vertices? (need to avoid v for cone)
3. For non-source v, does the cone still work on paths not using v?
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    boundary_coeffs
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def cone_map_vec(v, z, a2_t, a3_t, a3_idx, A, n):
    """Apply cone from v to vector z (on a2_t paths not using v).
    Returns w (vector on a3_t), or None if cone fails."""
    w = np.zeros(len(a3_t))
    for i, path in enumerate(a2_t):
        if abs(z[i]) < 1e-10:
            continue
        if v in path:
            return None
        if A[v][path[0]] != 1:
            return None
        cp = (v,) + tuple(path)
        if cp in a3_idx:
            w[a3_idx[cp]] += z[i]
        else:
            return None
    return w

# ===== Part 1: Verify cone from source =====
print("=" * 70)
print("PART 1: CONE FROM SOURCE VERTEX KILLS beta_2")
print("=" * 70)

n = 5
source_works = 0
has_source = 0
for A in all_tournaments_gen(n):
    source = None
    for v in range(n):
        if sum(A[v]) == n - 1:
            source = v
            break
    if source is None:
        continue

    has_source += 1
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]
    a3_idx = {p: i for i, p in enumerate(a3_t)}

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        source_works += 1
        continue

    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    _, S2, Vt2 = np.linalg.svd(bd2_om)
    r2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - r2
    if ker_dim == 0:
        source_works += 1
        continue

    ker_basis = om2 @ Vt2[r2:].T
    bd3 = build_full_boundary_matrix(a3_t, a2_t)

    all_ok = True
    for k in range(ker_dim):
        z = ker_basis[:, k]
        w = cone_map_vec(source, z, a2_t, a3_t, a3_idx, A, n)
        if w is None:
            all_ok = False
            break
        dw = bd3 @ w
        residual = np.max(np.abs(dw - z))
        if residual > 1e-6:
            all_ok = False
            break

    if all_ok:
        source_works += 1

print(f"  Tournaments with source vertex: {has_source}")
print(f"  Source cone fills all: {source_works}/{has_source}")

# ===== Part 2: Do 2-cycles avoid a vertex? =====
print(f"\n{'='*70}")
print("PART 2: DO 2-CYCLES AVOID AT LEAST ONE VERTEX?")
print("=" * 70)

for n in [5, 6]:
    avoids_stats = Counter()
    total_cycles = 0

    for tidx, A in enumerate(all_tournaments_gen(n)):
        if n == 6 and tidx % 5000 == 0:
            print(f"  n={n}: ... {tidx}", flush=True)

        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a1_t = [tuple(p) for p in a1]
        a2_t = [tuple(p) for p in a2]

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        if dim_om2 == 0:
            continue

        bd2 = build_full_boundary_matrix(a2_t, a1_t)
        bd2_om = bd2 @ om2
        _, S2, Vt2 = np.linalg.svd(bd2_om)
        r2 = sum(s > 1e-8 for s in S2)
        ker_dim = dim_om2 - r2
        if ker_dim == 0:
            continue

        ker_basis = om2 @ Vt2[r2:].T
        total_cycles += ker_dim

        for k in range(ker_dim):
            z = ker_basis[:, k]
            support_verts = set()
            for i in range(len(a2_t)):
                if abs(z[i]) > 1e-10:
                    support_verts.update(a2_t[i])
            n_avoids = n - len(support_verts)
            avoids_stats[n_avoids] += 1

    print(f"\n  n={n}: {total_cycles} total 2-cycles")
    print(f"  Vertices avoided per cycle:")
    for k in sorted(avoids_stats):
        pct = 100 * avoids_stats[k] / total_cycles
        print(f"    Avoids {k}: {avoids_stats[k]} ({pct:.1f}%)")

# ===== Part 3: Multi-vertex cone =====
print(f"\n{'='*70}")
print("PART 3: MULTI-VERTEX CONE (n=5)")
print("=" * 70)
print("For each 2-cycle z, can we find v not in support(z)")
print("such that v -> first vertex of every path in support(z)?")

n = 5
multi_cone_success = 0
multi_cone_fail = 0

for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]
    a3_idx = {p: i for i, p in enumerate(a3_t)}

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        multi_cone_success += 1
        continue

    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    _, S2, Vt2 = np.linalg.svd(bd2_om)
    r2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - r2
    if ker_dim == 0:
        multi_cone_success += 1
        continue

    ker_basis = om2 @ Vt2[r2:].T
    bd3 = build_full_boundary_matrix(a3_t, a2_t)

    # Collect all cone fillings
    filling_images = []
    for v in range(n):
        for k in range(ker_dim):
            z = ker_basis[:, k]
            w = cone_map_vec(v, z, a2_t, a3_t, a3_idx, A, n)
            if w is not None:
                dw = bd3 @ w
                if np.max(np.abs(dw - z)) < 1e-6:
                    filling_images.append(z.copy())

    if not filling_images:
        multi_cone_fail += 1
        continue

    fill_matrix = np.column_stack(filling_images)
    combined = np.hstack([fill_matrix, ker_basis])
    rf = np.linalg.matrix_rank(fill_matrix, tol=1e-8)
    rc = np.linalg.matrix_rank(combined, tol=1e-8)

    if rc == rf:
        multi_cone_success += 1
    else:
        multi_cone_fail += 1

print(f"  n={n}: Multi-cone success: {multi_cone_success}")
print(f"  n={n}: Multi-cone failure: {multi_cone_fail}")

# ===== Part 4: Does the cone formula d(s_v(z)) = z actually hold? =====
print(f"\n{'='*70}")
print("PART 4: VERIFY CONE FORMULA d(s_v(z)) = z - s_v(dz)")
print("=" * 70)
print("For z in ker d_2, dz = 0, so d(s_v(z)) should = z.")
print("This only works if z is supported on paths NOT using v")
print("AND v -> path[0] for every path in support.")

n = 5
for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx > 10:
        break

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]
    a3_idx = {p: i for i, p in enumerate(a3_t)}

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    _, S2, Vt2 = np.linalg.svd(bd2_om)
    r2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - r2
    if ker_dim == 0:
        continue

    ker_basis = om2 @ Vt2[r2:].T
    bd3 = build_full_boundary_matrix(a3_t, a2_t)

    print(f"\n  T#{tidx}: ker_dim={ker_dim}")
    for k in range(min(ker_dim, 2)):
        z = ker_basis[:, k]
        support_verts = set()
        for i in range(len(a2_t)):
            if abs(z[i]) > 1e-10:
                support_verts.update(a2_t[i])

        for v in range(n):
            if v in support_verts:
                continue
            w = cone_map_vec(v, z, a2_t, a3_t, a3_idx, A, n)
            if w is not None:
                dw = bd3 @ w
                print(f"    cycle {k}, cone v={v}: |dw-z|={np.max(np.abs(dw-z)):.2e}")
                break

print("\nDone.")
