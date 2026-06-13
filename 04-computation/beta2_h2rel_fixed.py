#!/usr/bin/env python3
"""
beta2_h2rel_fixed.py — Fixed H₂^rel computation

BUG FOUND: When deleting vertex v from the middle (not the last vertex),
T\v vertices get renumbered 0..n-2. Path (a,b,c) in T\v local numbering
corresponds to different global vertex indices. The embedding function
was NOT applying this renumbering, causing paths to not be found in T.

FIX: Build a local→global vertex map and apply it when embedding.

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
    """Delete vertex v from graph A, returning (n-1)×(n-1) adjacency."""
    B = []
    for i in range(n):
        if i == v:
            continue
        row = []
        for j in range(n):
            if j == v:
                continue
            row.append(A[i][j])
        B.append(row)
    return B


def local_to_global_map(n, v):
    """Map from T\v local vertex indices to T global indices.

    After deleting v from {0,...,n-1}, remaining vertices are
    renumbered 0,...,n-2. This returns the inverse mapping.
    """
    mapping = []
    for i in range(n):
        if i != v:
            mapping.append(i)
    return mapping  # mapping[local] = global


def embed_Tv_paths_in_T(allowed_T, allowed_Tv, omega_Tv, n, v):
    """Embed Ω_p(T\v) into A_p(T) coordinates with correct vertex renumbering.

    Returns matrix of shape |A_p(T)| × dim_Ω_p(T\v).
    """
    if omega_Tv.ndim != 2 or omega_Tv.shape[1] == 0:
        return np.zeros((len(allowed_T), 0))

    loc2glob = local_to_global_map(n, v)

    T_list = [tuple(x) for x in allowed_T]
    Tv_list = [tuple(x) for x in allowed_Tv]
    T_idx = {p: i for i, p in enumerate(T_list)}

    # Inclusion map: A_p(T\v) → A_p(T) with renumbering
    incl = np.zeros((len(T_list), len(Tv_list)))
    for j, path_local in enumerate(Tv_list):
        path_global = tuple(loc2glob[k] for k in path_local)
        if path_global in T_idx:
            incl[T_idx[path_global], j] = 1
        else:
            # This should never happen for vertex deletion in tournaments
            print(f"  WARNING: path {path_local} (global {path_global}) NOT in T!")

    # Ω_p(T\v) in A_p(T) coordinates
    return incl @ omega_Tv


def compute_h2_rel_fixed(A, n, v):
    """Compute H₂(T, T\v) with correct vertex renumbering."""
    B = delete_vertex(A, n, v)

    # Compute Ω bases
    om_T = {}
    ap_T = {}
    om_Tv = {}
    ap_Tv = {}

    for p in range(5):
        ap_T[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om_T[p] = np.eye(n)
        elif ap_T[p]:
            om_T[p] = compute_omega_basis(A, n, p, ap_T[p], ap_T[p - 1])
        else:
            om_T[p] = np.zeros((0, 0))

    for p in range(4):
        ap_Tv[p] = enumerate_allowed_paths(B, n - 1, p)
        if p == 0:
            om_Tv[p] = np.eye(n - 1)
        elif ap_Tv[p]:
            om_Tv[p] = compute_omega_basis(B, n - 1, p, ap_Tv[p], ap_Tv[p - 1])
        else:
            om_Tv[p] = np.zeros((0, 0))

    def dim_om(om):
        return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

    d2T = dim_om(om_T[2])
    d3T = dim_om(om_T[3])
    d2Tv = dim_om(om_Tv[2])

    if d2T == 0:
        return 0

    # Embed Ω_p(T\v) into A_p(T) with correct renumbering
    emb1 = embed_Tv_paths_in_T(ap_T[1], ap_Tv[1], om_Tv[1], n, v)
    emb2 = embed_Tv_paths_in_T(ap_T[2], ap_Tv[2], om_Tv[2], n, v)

    # Boundary map ∂₂: A₂(T) → A₁(T) restricted to Ω₂(T)
    bd2_A = build_full_boundary_matrix(ap_T[2], ap_T[1])
    bd2_om = bd2_A @ om_T[2]  # |A₁(T)| × d2T

    # ker(∂₂^R): {x ∈ Ω₂(T) : ∂₂(x) ∈ Ω₁(T\v)} / Ω₂(T\v)
    # x = om_T[2] @ c, ∂₂(x) = bd2_om @ c
    # Need bd2_om @ c ∈ col_span(emb1)
    if emb1.shape[1] > 0:
        M = np.hstack([bd2_om, -emb1])
    else:
        M = bd2_om

    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    tol = 1e-8
    rk_M = int(sum(s > tol for s in S))
    if rk_M < Vt.shape[0]:
        null_space = Vt[rk_M:].T
    else:
        null_space = np.zeros((Vt.shape[1], 0))

    c_part = null_space[:d2T, :]
    preimage_A = om_T[2] @ c_part

    # ker(∂₂^R) = dim(span(preimage ∪ emb2)) - dim(span(emb2))
    if emb2.shape[1] > 0 and preimage_A.shape[1] > 0:
        combined = np.hstack([emb2, preimage_A])
    elif preimage_A.shape[1] > 0:
        combined = preimage_A
    else:
        combined = emb2

    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8) if combined.shape[1] > 0 else 0
    rk_emb2 = np.linalg.matrix_rank(emb2, tol=1e-8) if emb2.shape[1] > 0 else 0
    dim_ker_R2 = rk_combined - rk_emb2

    # im(∂₃^R) = (∂₃(Ω₃(T)) + Ω₂(T\v)) / Ω₂(T\v)
    if d3T > 0:
        bd3_A = build_full_boundary_matrix(ap_T[3], ap_T[2])
        bd3_om = bd3_A @ om_T[3]
        if emb2.shape[1] > 0:
            combined_3 = np.hstack([emb2, bd3_om])
        else:
            combined_3 = bd3_om
        rk_combined_3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
        dim_im_R3 = rk_combined_3 - rk_emb2
    else:
        dim_im_R3 = 0

    return dim_ker_R2 - dim_im_R3


def compute_betti_1(A, n):
    """Compute β₁ only."""
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    if not ap1:
        return 0
    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    d1 = om1.shape[1] if om1.ndim == 2 else 0

    bd1 = build_full_boundary_matrix(ap1, enumerate_allowed_paths(A, n, 0))
    bd1_om = bd1 @ om1
    S1 = np.linalg.svd(np.linalg.lstsq(np.eye(n), bd1_om, rcond=None)[0], compute_uv=False)
    rk1 = int(sum(s > 1e-8 for s in S1))

    if ap2:
        om2 = compute_omega_basis(A, n, 2, ap2, ap1)
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            S2 = np.linalg.svd(np.linalg.lstsq(om1, bd2_om, rcond=None)[0], compute_uv=False)
            rk2 = int(sum(s > 1e-8 for s in S2))
        else:
            rk2 = 0
    else:
        rk2 = 0

    return (d1 - rk1) - rk2


# ===== MAIN =====
print("=" * 70)
print("FIXED H₂^rel: vertex renumbering corrected")
print("=" * 70)

n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
m = len(pairs)

# Quick test on known problem case
print("\nTest: bits=9, v=1")
A_test = [[0] * n for _ in range(n)]
for idx, (i, j) in enumerate(pairs):
    if (9 >> idx) & 1:
        A_test[i][j] = 1
    else:
        A_test[j][i] = 1

h2r = compute_h2_rel_fixed(A_test, n, 1)
b1_T = compute_betti_1(A_test, n)
B_test = delete_vertex(A_test, n, 1)
b1_Tv = compute_betti_1(B_test, n - 1)
print(f"  b1_T={b1_T}, b1_Tv={b1_Tv}, H2_rel={h2r}")
if b1_Tv > b1_T and h2r >= 1:
    print("  FIXED! ✓")
elif b1_Tv > b1_T and h2r < 1:
    print("  STILL BROKEN ✗")
else:
    print("  (no LES constraint)")

# Systematic check
print(f"\nSystematic check: all n={n} (T,v) pairs")
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

        h2r = compute_h2_rel_fixed(A, n, v)
        h2_dist[h2r] += 1
        total += 1

        les_min = max(0, b1_Tv - b1_T)
        if h2r < les_min:
            bugs += 1
            if bugs <= 3:
                print(f"  BUG: bits={bits}, v={v}, h2_rel={h2r}, "
                      f"b1_T={b1_T}, b1_Tv={b1_Tv}")
        elif h2r >= les_min:
            les_check['ok'] += 1

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}, bugs={bugs}")

print(f"\nH₂^rel distribution ({total} pairs):")
for val, count in sorted(h2_dist.items()):
    print(f"  H₂^rel = {val}: {count}")

print(f"\nLES violations: {bugs}")
if bugs == 0:
    print("ALL CONSISTENT WITH LES! ✓")

    # Detailed breakdown
    print(f"\nBreakdown by (β₁(T), β₁(T\\v)):")
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
            h2r = compute_h2_rel_fixed(A, n, v)
            breakdown[(b1_T, b1_Tv, h2r)] += 1

    for key, count in sorted(breakdown.items()):
        print(f"  β₁(T)={key[0]}, β₁(T\\v)={key[1]}, H₂^rel={key[2]}: {count}")

    # KEY CHECK: Is H₂^rel ALWAYS = max(0, β₁(T\v) - β₁(T))?
    print(f"\nIs H₂^rel = max(0, β₁(T\\v) - β₁(T)) always?")
    exact_match = 0
    mismatch = 0
    for bits in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        b1_T = compute_betti_1(A, n)
        for v_idx in range(n):
            B = delete_vertex(A, n, v_idx)
            b1_Tv = compute_betti_1(B, n - 1)
            h2r = compute_h2_rel_fixed(A, n, v_idx)
            expected = max(0, b1_Tv - b1_T)
            if h2r == expected:
                exact_match += 1
            else:
                mismatch += 1
                if mismatch <= 5:
                    print(f"  MISMATCH: bits={bits}, v={v_idx}, H₂^rel={h2r}, "
                          f"expected={expected}")

    print(f"  Exact match: {exact_match}/{total}")
    print(f"  Mismatch: {mismatch}/{total}")
else:
    print(f"Still {bugs} violations — deeper issue.")

print("\nDone.")
