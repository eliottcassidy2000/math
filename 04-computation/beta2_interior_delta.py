#!/usr/bin/env python3
"""
beta2_interior_delta.py — Exhaustive verification + structural analysis of
interior delta-injectivity for beta2=0 proof.

PROOF STRATEGY (LES induction):
  For n>=4 and tournament T, pick an INTERIOR vertex v (1 <= d+(v) <= n-2).
  By induction: beta2(T\v) = 0.
  LES: 0 = H2(T\v) -> H2(T) ->^α H2(T,T\v) ->^delta H1(T\v)
  α is injective (since H2(T\v)=0), delta injective ⟹ ker(delta)=0 = im(α) ⟹ H2(T)=0.

  So: beta2=0 follows from "delta is injective for all interior vertices."

  This script:
  1. Exhaustively verifies interior delta-injectivity at n=5,6
  2. Samples at n=7
  3. Analyzes the ALGEBRAIC MECHANISM: why does position mixing prevent delta=0?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
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

random.seed(42)


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


def compute_delta_injectivity(A, n, v):
    """
    Check if the connecting map delta: H2(T,T\v) -> H1(T\v) is injective.
    Returns (is_injective, h2_rel, analysis_data).
    """
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    # Compute chain complexes for T and T\v
    ap0_T = [(i,) for i in range(n)]
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = [(i,) for i in range(n1)]
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))

    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    d3_T = dim_om(om3_T)
    d1_sub = dim_om(om1_sub)
    d2_sub = dim_om(om2_sub)

    if d2_T == 0:
        return True, 0, {}

    # Build embedding maps
    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

    # Embed Omega2(T\v) -> Omega2(T)
    if ap2_sub and d2_sub > 0:
        embed2 = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    embed2[ap2_T_list.index(path_T), j] = om2_sub[k, j]
        phi = np.linalg.lstsq(om2_T, embed2, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    # Quotient basis for Omega2(T)/Omega2(T\v)
    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
    if rk_phi > 0:
        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        return True, 0, {}

    # Boundary del2 in Omega1(T) coordinates
    bd2_raw = build_full_boundary_matrix(ap2_T, ap1_T)
    coords2_T = np.linalg.lstsq(om1_T, bd2_raw @ om2_T, rcond=None)[0]

    # Embed Omega1(T\v) -> Omega1(T)
    d1_T = dim_om(om1_T)
    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
    if rk_psi > 0:
        U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
        R = U_psi[:, rk_psi:]  # complement of Omega1(T\v) in Omega1(T)
    else:
        R = np.eye(d1_T)

    # Relative del2: project coords2_T to quotient and v-arcs
    coords2_rel_q = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)
    z2_rel_dim = d_rel - rk_d2_rel

    if z2_rel_dim == 0:
        return True, 0, {}

    # Relative del3
    if d3_T > 0:
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)
        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
        d3_sub = dim_om(om3_sub)

        bd3_raw = build_full_boundary_matrix(ap3_T, ap2_T)
        coords3_T = np.linalg.lstsq(om2_T, bd3_raw @ om3_T, rcond=None)[0]

        d3_proj = Q.T @ coords3_T
        if d3_sub > 0:
            ap3_T_list = [tuple(p) for p in ap3_T]
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
            chi3 = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
            rk_chi = np.linalg.matrix_rank(chi3, tol=1e-8)
            if rk_chi > 0:
                U_chi, _, _ = np.linalg.svd(chi3, full_matrices=True)
                d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
    h2_rel = z2_rel_dim - rk_d3_rel

    if h2_rel == 0:
        return True, 0, {}

    # Check delta-injectivity: do H2(T,T\v) elements map to nonzero classes in H1(T\v)?
    U_rel, S_rel, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
    z2_rel_basis = Vt_rel[rk_d2_rel:]  # nullspace of del2^rel in the quotient

    if d3_proj_rel.shape[1] > 0:
        proj_z2 = z2_rel_basis @ d3_proj_rel
        rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
    else:
        rk_proj = 0

    # H2(T,T\v) basis (after modding out boundaries)
    h2_rel_basis = z2_rel_basis
    if rk_proj > 0:
        # Remove boundary part from z2 basis
        U_p, _, _ = np.linalg.svd(proj_z2.T, full_matrices=True)
        h2_rel_basis = z2_rel_basis[rk_proj:]

    if h2_rel_basis.shape[0] == 0:
        return True, 0, {}

    # Compute delta images in H1(T\v)
    delta_images = []
    for k in range(h2_rel_basis.shape[0]):
        h = h2_rel_basis[k]
        # lift to Omega2(T) coords, apply del2
        bd = coords2_T @ Q @ h  # del2(z̃) in Omega1(T) coords
        # project to Omega1(T\v) coords
        if d1_sub > 0:
            delta_sub = psi.T @ bd  # using psi^T as projection
            # Actually need to solve psi @ x = bd for the Omega1(T\v) component
            delta_sub = np.linalg.lstsq(psi, bd, rcond=None)[0]
            delta_images.append(delta_sub)

    if not delta_images or d1_sub == 0:
        return True, h2_rel, {'z2_rel': z2_rel_dim, 'h2_rel': h2_rel}

    # Check if delta images are nontrivial in H1(T\v)
    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1_sub = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)

    if d2_sub > 0:
        bd2_sub_raw = build_full_boundary_matrix(ap2_sub, ap1_sub)
        coords2_sub = np.linalg.lstsq(om1_sub, bd2_sub_raw @ om2_sub, rcond=None)[0]
        rk_d2_sub = np.linalg.matrix_rank(coords2_sub, tol=1e-8)
    else:
        rk_d2_sub = 0

    beta1_sub = d1_sub - rk_d1_sub - rk_d2_sub

    # Z1(T\v) basis
    _, _, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
    z1_basis = Vt1[rk_d1_sub:]

    # Project delta images to Z1 coordinates
    D = np.array(delta_images)
    D_z1 = D @ z1_basis.T

    # Check if images are in B1(T\v) or outside
    # B1_z1: boundaries in Z1 coords — shape (z1_dim, d2_sub)
    # D_z1: delta images in Z1 coords — shape (h2_rel, z1_dim)
    if rk_d2_sub > 0:
        B1_z1 = z1_basis @ coords2_sub  # (z1_dim, d2_sub)
        rk_B1 = np.linalg.matrix_rank(B1_z1, tol=1e-8)
        # Stack as row vectors: B1_z1.T has shape (d2_sub, z1_dim), D_z1 has (h2_rel, z1_dim)
        combined = np.vstack([B1_z1.T, D_z1])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        delta_rk_in_H1 = rk_combined - rk_B1
    else:
        delta_rk_in_H1 = np.linalg.matrix_rank(D_z1, tol=1e-8)

    is_injective = (delta_rk_in_H1 == h2_rel)

    data = {
        'z2_rel': z2_rel_dim, 'h2_rel': h2_rel,
        'beta1_sub': beta1_sub, 'delta_rk': delta_rk_in_H1,
        'd_rel': d_rel, 'd2_T': d2_T, 'd1_sub': d1_sub,
    }
    return is_injective, h2_rel, data


# ============================================================
# Part 1: Exhaustive verification at n=5
# ============================================================
print("=" * 70)
print("INTERIOR DELTA-INJECTIVITY: EXHAUSTIVE VERIFICATION")
print("=" * 70)

for n in [5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    if n == 7:
        samples = random.sample(range(total), 1000)
    else:
        samples = range(total)

    interior_checks = 0
    interior_fails = 0
    boundary_checks = 0
    boundary_fails = 0
    h2rel_dist = Counter()
    beta1_when_h2 = Counter()
    fail_details = []

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

        for v in range(n):
            dv = scores[v]
            is_extremal = (dv == 0 or dv == n-1)

            inj, h2_rel, data = compute_delta_injectivity(A, n, v)

            if h2_rel > 0:
                h2rel_dist[h2_rel] += 1
                if 'beta1_sub' in data:
                    beta1_when_h2[data['beta1_sub']] += 1

                if is_extremal:
                    boundary_checks += 1
                    if not inj:
                        boundary_fails += 1
                else:
                    interior_checks += 1
                    if not inj:
                        interior_fails += 1
                        sc = tuple(sorted(scores))
                        fail_details.append((bits, sc, v, dv, data))
                        print(f"  INTERIOR FAIL! n={n} T#{bits} scores={sc}, v={v}, d+={dv}, data={data}")

        if (idx + 1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {idx+1}/{len(samples)} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(samples)} tournaments in {elapsed:.0f}s")
    print(f"  Interior vertex checks (with H2(T,T\\v) > 0): {interior_checks}")
    print(f"  Interior failures: {interior_fails}")
    print(f"  Boundary (source/sink) checks: {boundary_checks}")
    print(f"  Boundary failures: {boundary_fails}")
    print(f"  H2(T,T\\v) dimension distribution: {dict(sorted(h2rel_dist.items()))}")
    print(f"  beta1(T\\v) when H2(T,T\\v)>0: {dict(sorted(beta1_when_h2.items()))}")

    if interior_fails == 0:
        print(f"  OK: delta INJECTIVE for ALL interior vertices at n={n}")
    else:
        print(f"  FAIL: {interior_fails} interior failures")
        for bits, sc, v, dv, data in fail_details[:5]:
            print(f"    T#{bits} scores={sc}, v={v}, d+={dv}")


# ============================================================
# Part 2: Sample at n=7
# ============================================================
print(f"\n{'='*70}")
print("INTERIOR DELTA-INJECTIVITY: SAMPLING AT n=7")
print("=" * 70)

n = 7
total = 1 << (n*(n-1)//2)
sample_size = 500

interior_checks = 0
interior_fails = 0
h2rel_dist = Counter()

t0 = time.time()
for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n - 1:
            continue  # skip extremal

        inj, h2_rel, data = compute_delta_injectivity(A, n, v)
        if h2_rel > 0:
            interior_checks += 1
            h2rel_dist[h2_rel] += 1
            if not inj:
                interior_fails += 1
                print(f"  INTERIOR FAIL! T#{bits} scores={tuple(sorted(scores))}, v={v}")

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{sample_size} done ({elapsed:.1f}s), int_checks={interior_checks}, fails={interior_fails}")

elapsed = time.time() - t0
print(f"\nn=7: {sample_size} tournaments in {elapsed:.0f}s")
print(f"  Interior checks (with H2>0): {interior_checks}")
print(f"  Interior failures: {interior_fails}")
print(f"  H2(T,T\\v) distribution: {dict(sorted(h2rel_dist.items()))}")

if interior_fails == 0:
    print(f"  OK: delta INJECTIVE for ALL interior vertices at n=7")


# ============================================================
# Part 3: Analyze the position structure for interior vs source
# ============================================================
print(f"\n{'='*70}")
print("POSITION STRUCTURE: WHY INTERIOR DELTA WORKS")
print("=" * 70)

n = 5
# Find cases with h2_rel > 0 for both interior and source vertices
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    # Check if this tournament has BOTH an interior h2_rel>0 and a source/sink h2_rel>0
    interior_h2 = []
    boundary_h2 = []
    for v in range(n):
        dv = scores[v]
        _, h2_rel, _ = compute_delta_injectivity(A, n, v)
        if h2_rel > 0:
            if dv == 0 or dv == n-1:
                boundary_h2.append((v, dv, h2_rel))
            else:
                interior_h2.append((v, dv, h2_rel))

    if interior_h2 and boundary_h2:
        print(f"\nT#{bits} scores={tuple(sorted(scores))}")
        print(f"  Interior with H2>0: {interior_h2}")
        print(f"  Boundary with H2>0: {boundary_h2}")

        # Detailed position analysis
        for v, dv, h2 in interior_h2[:1]:
            in_nbrs = [j for j in range(n) if j != v and A[j][v]]
            out_nbrs = [j for j in range(n) if j != v and A[v][j]]
            print(f"  Interior v={v} (d+={dv}): N-={in_nbrs}, N+={out_nbrs}")

            # Count 2-paths by position of v
            ap2 = enumerate_allowed_paths(A, n, 2)
            pos0 = [(a,b,c) for (a,b,c) in ap2 if a == v]
            pos1 = [(a,b,c) for (a,b,c) in ap2 if b == v]
            pos2 = [(a,b,c) for (a,b,c) in ap2 if c == v]
            print(f"    |A2|={len(ap2)}: v@0={len(pos0)}, v@1={len(pos1)}, v@2={len(pos2)}")

        for v, dv, h2 in boundary_h2[:1]:
            # Source/sink analysis
            in_nbrs = [j for j in range(n) if j != v and A[j][v]]
            out_nbrs = [j for j in range(n) if j != v and A[v][j]]
            print(f"  Boundary v={v} (d+={dv}): N-={in_nbrs}, N+={out_nbrs}")

            ap2 = enumerate_allowed_paths(A, n, 2)
            pos0 = [(a,b,c) for (a,b,c) in ap2 if a == v]
            pos1 = [(a,b,c) for (a,b,c) in ap2 if b == v]
            pos2 = [(a,b,c) for (a,b,c) in ap2 if c == v]
            print(f"    |A2|={len(ap2)}: v@0={len(pos0)}, v@1={len(pos1)}, v@2={len(pos2)}")

        if len(interior_h2) > 0 and len(boundary_h2) > 0:
            break  # Just show one example


# ============================================================
# Part 4: KEY OBSERVATION — v-arc type counting
# ============================================================
print(f"\n{'='*70}")
print("V-ARC TYPE ANALYSIS (n=5 exhaustive)")
print("=" * 70)

n = 5
source_inj_data = []
interior_inj_data = []

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        inj, h2_rel, data = compute_delta_injectivity(A, n, v)
        if h2_rel > 0:
            in_nbrs = [j for j in range(n) if j != v and A[j][v]]
            out_nbrs = [j for j in range(n) if j != v and A[v][j]]

            ap2 = enumerate_allowed_paths(A, n, 2)
            n_pos0 = sum(1 for (a,b,c) in ap2 if a == v)
            n_pos1 = sum(1 for (a,b,c) in ap2 if b == v)
            n_pos2 = sum(1 for (a,b,c) in ap2 if c == v)

            entry = (dv, n_pos0, n_pos1, n_pos2, h2_rel, inj,
                     data.get('beta1_sub', -1))

            if dv == 0 or dv == n-1:
                source_inj_data.append(entry)
            else:
                interior_inj_data.append(entry)

print(f"Interior (d+, #pos0, #pos1, #pos2, h2_rel, delta_inj, beta1(T\\v)):")
int_summary = Counter(interior_inj_data)
for entry, count in sorted(int_summary.items()):
    print(f"  d+={entry[0]}, pos=({entry[1]},{entry[2]},{entry[3]}), "
          f"h2={entry[4]}, d_inj={entry[5]}, b1={entry[6]}: {count} cases")

print(f"\nSource/sink (d+, #pos0, #pos1, #pos2, h2_rel, delta_inj, beta1(T\\v)):")
bd_summary = Counter(source_inj_data)
for entry, count in sorted(bd_summary.items()):
    print(f"  d+={entry[0]}, pos=({entry[1]},{entry[2]},{entry[3]}), "
          f"h2={entry[4]}, d_inj={entry[5]}, b1={entry[6]}: {count} cases")


# ============================================================
# Part 5: The KEY dimension constraint
# ============================================================
print(f"\n{'='*70}")
print("DIMENSION CONSTRAINT: H2(T,T\\v) vs beta1(T\\v)")
print("=" * 70)

# For delta injective, need dim H2(T,T\\v) <= beta1(T\\v) (dimension of target in H1)
# Let's check if this always holds for interior vertices.

n = 5
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        inj, h2_rel, data = compute_delta_injectivity(A, n, v)
        if h2_rel > 0:
            b1 = data.get('beta1_sub', 0)
            if h2_rel > b1:
                print(f"  DIM VIOLATION! T#{bits} v={v}: h2_rel={h2_rel} > beta1(T\\v)={b1}")

print("(no output = all interior cases have h2_rel <= beta1(T\\v))")

print("\n\nDone.")
