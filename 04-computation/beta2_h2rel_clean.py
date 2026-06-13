#!/usr/bin/env python3
"""
beta2_h2rel_clean.py — Clean H₂^rel computation in A-coordinates

The previous compute_h2_rel() had 426 incorrect cases (gave 0 when LES demands 1).
Root cause: lstsq coordinate transformations lose precision.

Fix: work entirely in A-coordinates (raw path vector spaces). No lstsq needed.

The quotient complex R_p = Ω_p(T) / Ω_p(T\v) has:
  ker(∂₂^R) = {x ∈ Ω₂(T) : ∂₂(x) ∈ Ω₁(T\v)} / Ω₂(T\v)
  im(∂₃^R) = (∂₃(Ω₃(T)) + Ω₂(T\v)) / Ω₂(T\v)

Everything lives in A_p(T) coordinates — integer vectors, no basis changes needed.

Author: opus-2026-03-08-S45
"""
import sys
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


def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B


def get_omega_in_A_coords(A, n, p):
    """Return Ω_p basis as columns in A_p coordinate space."""
    ap = enumerate_allowed_paths(A, n, p)
    if p == 0:
        return np.eye(n), ap
    if not ap:
        return np.zeros((0, 0)), ap
    ap_prev = enumerate_allowed_paths(A, n, p - 1)
    om = compute_omega_basis(A, n, p, ap, ap_prev)
    if om.ndim != 2 or om.shape[0] == 0:
        return np.zeros((len(ap), 0)), ap
    return om, ap  # |A_p| × dim_Ω_p


def embed_Tv_in_T(allowed_T, allowed_Tv, omega_Tv):
    """Embed Ω_p(T\v) into A_p(T) coordinates.

    T\v paths use same vertex labels (just missing v).
    Returns matrix of shape |A_p(T)| × dim_Ω_p(T\v).
    """
    if omega_Tv.ndim != 2 or omega_Tv.shape[1] == 0:
        return np.zeros((len(allowed_T), 0))

    T_list = [tuple(x) for x in allowed_T]
    Tv_list = [tuple(x) for x in allowed_Tv]
    T_idx = {p: i for i, p in enumerate(T_list)}

    # Inclusion map: A_p(T\v) → A_p(T)
    incl = np.zeros((len(T_list), len(Tv_list)))
    for j, path in enumerate(Tv_list):
        if path in T_idx:
            incl[T_idx[path], j] = 1

    # Ω_p(T\v) in A_p(T) coords
    return incl @ omega_Tv


def compute_h2_rel_clean(A, n, v):
    """Compute H₂(T, T\v) cleanly in A-coordinates.

    Returns (h2_rel, debug_info).
    """
    B = delete_vertex(A, n, v)

    # Ω bases in A-coordinates
    om_T = {}
    ap_T = {}
    om_Tv = {}
    ap_Tv = {}

    for p in range(5):
        om_T[p], ap_T[p] = get_omega_in_A_coords(A, n, p)
    for p in range(4):
        om_Tv[p], ap_Tv[p] = get_omega_in_A_coords(B, n - 1, p)

    def dim(M):
        return M.shape[1] if M.ndim == 2 and M.shape[0] > 0 else 0

    d2T = dim(om_T[2])
    d3T = dim(om_T[3])
    d2Tv = dim(om_Tv[2])

    if d2T == 0:
        return 0, {}

    # Embed Ω_p(T\v) into A_p(T) coordinates
    emb = {}
    for p in [1, 2]:
        emb[p] = embed_Tv_in_T(ap_T[p], ap_Tv[p], om_Tv[p])

    # Boundary map ∂₂: A₂(T) → A₁(T) (integer matrix)
    bd2_A = build_full_boundary_matrix(ap_T[2], ap_T[1])

    # ∂₂ restricted to Ω₂(T): bd2_A @ om_T[2]
    # This gives the image of each Ω₂ basis vector in A₁(T) coords
    bd2_om = bd2_A @ om_T[2]  # |A₁(T)| × d2T

    # Embed Ω₁(T\v) in A₁(T): emb[1] columns
    # ker(∂₂^R): x ∈ Ω₂(T) with ∂₂(x) ∈ Ω₁(T\v) (embedded in A₁(T))
    # x = om_T[2] @ c for coefficient vector c
    # ∂₂(x) = bd2_om @ c
    # Need: bd2_om @ c ∈ col_span(emb[1])
    # i.e., bd2_om @ c = emb[1] @ y for some y
    # Null space of [bd2_om | -emb[1]] gives all (c, y)

    emb1 = emb[1]  # |A₁(T)| × d1Tv

    if emb1.shape[1] > 0:
        M = np.hstack([bd2_om, -emb1])
    else:
        # ∂₂(x) must be 0 (since Ω₁(T\v) = 0 means we need ∂₂(x) = 0)
        M = bd2_om

    # Null space via SVD
    U, S_vals, Vt = np.linalg.svd(M, full_matrices=True)
    tol = max(M.shape) * max(S_vals, default=0) * np.finfo(float).eps * 100
    if tol < 1e-10:
        tol = 1e-10
    rk_M = int(sum(s > tol for s in S_vals))

    if rk_M < Vt.shape[0]:
        null_space = Vt[rk_M:].T  # (d2T + d1Tv) × null_dim
    else:
        null_space = np.zeros((Vt.shape[1], 0))

    # Extract c-part (coefficients in Ω₂(T))
    c_part = null_space[:d2T, :]

    # These c vectors give elements of Ω₂(T) in the preimage.
    # The actual elements in A₂(T) are: om_T[2] @ c_part
    preimage_A = om_T[2] @ c_part  # |A₂(T)| × null_dim

    # Now ker(∂₂^R) = span(preimage_A) / span(emb[2])
    # dim(ker) = dim(span(preimage_A ∪ emb[2])) - dim(span(emb[2]))
    # But actually: dim(ker) = dim(preimage_A) - dim(preimage_A ∩ emb[2])
    # Which equals: rank([emb[2] | preimage_A]) - rank(emb[2])

    emb2 = emb[2]  # |A₂(T)| × d2Tv

    if emb2.shape[1] > 0 and preimage_A.shape[1] > 0:
        combined = np.hstack([emb2, preimage_A])
    elif preimage_A.shape[1] > 0:
        combined = preimage_A
    else:
        combined = emb2

    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
    rk_emb2 = np.linalg.matrix_rank(emb2, tol=1e-8) if emb2.shape[1] > 0 else 0
    dim_ker_R2 = rk_combined - rk_emb2

    # im(∂₃^R) = (∂₃(Ω₃(T)) + Ω₂(T\v)) / Ω₂(T\v)
    if d3T > 0:
        bd3_A = build_full_boundary_matrix(ap_T[3], ap_T[2])
        bd3_om = bd3_A @ om_T[3]  # |A₂(T)| × d3T — image of Ω₃ under ∂₃

        if emb2.shape[1] > 0:
            combined_3 = np.hstack([emb2, bd3_om])
        else:
            combined_3 = bd3_om
        rk_combined_3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
        dim_im_R3 = rk_combined_3 - rk_emb2
    else:
        dim_im_R3 = 0

    h2_rel = dim_ker_R2 - dim_im_R3

    debug = {
        'd2T': d2T, 'd3T': d3T, 'd2Tv': d2Tv,
        'dim_ker_R2': dim_ker_R2, 'dim_im_R3': dim_im_R3,
        'rk_preimage': np.linalg.matrix_rank(preimage_A, tol=1e-8) if preimage_A.shape[1] > 0 else 0,
        'null_dim': null_space.shape[1],
    }

    return h2_rel, debug


def compute_betti_1(A, n):
    """Compute β₁ only (fast)."""
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap2 = enumerate_allowed_paths(A, n, 2)

    if not ap1:
        return 0

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    d1 = om1.shape[1] if om1.ndim == 2 else 0

    # rank(∂₁)
    bd1 = build_full_boundary_matrix(ap1, ap0)
    bd1_om = bd1 @ om1
    om0 = np.eye(n)
    coords1, _, _, _ = np.linalg.lstsq(om0, bd1_om, rcond=None)
    S1 = np.linalg.svd(coords1, compute_uv=False)
    rk1 = int(sum(s > 1e-8 for s in S1))

    # rank(∂₂)
    if ap2:
        om2 = compute_omega_basis(A, n, 2, ap2, ap1)
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
            S2 = np.linalg.svd(coords2, compute_uv=False)
            rk2 = int(sum(s > 1e-8 for s in S2))
        else:
            rk2 = 0
    else:
        rk2 = 0

    return (d1 - rk1) - rk2


# ===== MAIN =====
print("=" * 70)
print("CLEAN H₂^rel COMPUTATION (A-coordinates only)")
print("=" * 70)

# First test on the specific case that should have H₂^rel = 1
n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

# Test case: bits=9, v=1 (known problem case)
A_test = [[0] * n for _ in range(n)]
for idx, (i, j) in enumerate(pairs):
    if (9 >> idx) & 1:
        A_test[i][j] = 1
    else:
        A_test[j][i] = 1

v_test = 1
print(f"\nTest case: bits=9, v={v_test}")
for i in range(n):
    arcs = [j for j in range(n) if A_test[i][j]]
    print(f"  {i} -> {arcs}")

B_test = delete_vertex(A_test, n, v_test)
b1_T = compute_betti_1(A_test, n)
b1_Tv = compute_betti_1(B_test, n - 1)
print(f"  beta1(T)={b1_T}, beta1(T\\v)={b1_Tv}")
print(f"  LES demands H2_rel >= {max(0, b1_Tv - b1_T)}")

h2r, dbg = compute_h2_rel_clean(A_test, n, v_test)
print(f"  H2_rel (clean) = {h2r}")
print(f"  Debug: {dbg}")

# Now systematic check
print(f"\n{'='*70}")
print(f"SYSTEMATIC: ALL n=5 (T,v) pairs")
print("=" * 70)

m = len(pairs)
h2_dist = Counter()
les_check = Counter()
bugs = 0
total = 0

for bits in range(1 << m):
    A = [[0] * n for _ in range(n)]
    for idx, (i, j) in enumerate(pairs):
        if (bits >> idx) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    b1_T = compute_betti_1(A, n)

    for v in range(n):
        B = delete_vertex(A, n, v)
        b1_Tv = compute_betti_1(B, n - 1)

        h2r, _ = compute_h2_rel_clean(A, n, v)
        h2_dist[h2r] += 1
        total += 1

        # LES check
        les_min = max(0, b1_Tv - b1_T)
        if h2r < les_min:
            bugs += 1
            les_check['violation'] += 1
            if bugs <= 3:
                print(f"  BUG: bits={bits}, v={v}, h2_rel={h2r}, "
                      f"b1_T={b1_T}, b1_Tv={b1_Tv}, les_min={les_min}")
        else:
            les_check['ok'] += 1

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}, bugs={bugs}")

print(f"\nH2_rel distribution ({total} pairs):")
for val, count in sorted(h2_dist.items()):
    print(f"  H2_rel = {val}: {count}")

print(f"\nLES consistency: {les_check}")
if bugs == 0:
    print("ALL cases consistent with LES!")

    # Cross-check: H2_rel should equal dim(ker(i_*: H1(Tv) -> H1(T)))
    # When beta1(T)=0 and beta1(Tv)=1, H2_rel should be exactly 1.
    # When beta1(T)=beta1(Tv), H2_rel could be 0.
    print(f"\nBreakdown by (beta1_T, beta1_Tv):")
    breakdown = Counter()
    for bits in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        b1_T = compute_betti_1(A, n)
        for v in range(n):
            B = delete_vertex(A, n, v)
            b1_Tv = compute_betti_1(B, n - 1)
            h2r, _ = compute_h2_rel_clean(A, n, v)
            breakdown[(b1_T, b1_Tv, h2r)] += 1

    for key, count in sorted(breakdown.items()):
        print(f"  b1_T={key[0]}, b1_Tv={key[1]}, H2_rel={key[2]}: {count}")
else:
    print(f"STILL {bugs} violations!")

print("\nDone.")
