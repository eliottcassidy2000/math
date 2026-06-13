"""
istar_b3b4_fast.py — Find and analyze beta_3=beta_4=1 tournaments quickly

Two-phase: first find them with quick Betti check, then do full i_* analysis.
Uses the same seed kind-pasteur used in their S48 run.

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    full_chain_complex_modp,
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_istar_p3(A, n, v, max_p=7):
    """Compute rank(i_*^3) specifically. Faster than all-degrees version."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    data_T = full_chain_complex_modp(A, n, max_p)
    data_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    bT = data_T['bettis']
    bTv = data_Tv['bettis']

    if bTv.get(3, 0) == 0 or bT.get(3, 0) == 0:
        return 0, bT, bTv

    # Compute rank(i_*^3) using the standard method
    p = 3
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    paths_T_p = ap_T.get(p, [])
    paths_Tv_p = ap_Tv.get(p, [])
    paths_T_pp1 = ap_T.get(p + 1, [])
    paths_Tv_pm1 = ap_Tv.get(p - 1, [])

    # Omega basis for T at degree 3
    P_T, na_rows_T, na_cols_T = _build_constraint_matrix(ap_T, p, PRIME)
    if P_T is not None:
        _, nbasis_T = _gauss_nullbasis_modp(P_T, na_rows_T, na_cols_T, PRIME)
        omega_basis_T = np.array(nbasis_T, dtype=np.int64) if nbasis_T else None
    else:
        omega_basis_T = np.eye(len(paths_T_p), dtype=np.int64)

    if omega_basis_T is None:
        return 0, bT, bTv

    # Omega basis for T\v at degree 3
    P_Tv, na_rows_Tv, na_cols_Tv = _build_constraint_matrix(ap_Tv, p, PRIME)
    if P_Tv is not None:
        _, nbasis_Tv = _gauss_nullbasis_modp(P_Tv, na_rows_Tv, na_cols_Tv, PRIME)
        omega_basis_Tv = np.array(nbasis_Tv, dtype=np.int64) if nbasis_Tv else None
    else:
        omega_basis_Tv = np.eye(len(paths_Tv_p), dtype=np.int64)

    if omega_basis_Tv is None:
        return 0, bT, bTv

    # im(d_4^T) in A_3(T) coords
    idx_Tp = {tuple(q): i for i, q in enumerate(paths_T_p)}
    if paths_T_pp1:
        bd_T_pp1 = np.zeros((len(paths_T_p), len(paths_T_pp1)), dtype=np.int64)
        for j, sigma in enumerate(paths_T_pp1):
            for sign, face in boundary_faces(sigma):
                face_t = tuple(face)
                if face_t in idx_Tp:
                    bd_T_pp1[idx_Tp[face_t], j] = (bd_T_pp1[idx_Tp[face_t], j] + sign) % PRIME

        P_Tpp1, nrpp1, ncpp1 = _build_constraint_matrix(ap_T, p + 1, PRIME)
        if P_Tpp1 is not None:
            rpp1, nb_pp1 = _gauss_nullbasis_modp(P_Tpp1, nrpp1, ncpp1, PRIME)
            ob_pp1 = np.array(nb_pp1, dtype=np.int64) if nb_pp1 else None
        else:
            ob_pp1 = np.eye(len(paths_T_pp1), dtype=np.int64)

        if ob_pp1 is not None:
            im_dp1 = bd_T_pp1 @ ob_pp1.T % PRIME
        else:
            im_dp1 = np.zeros((len(paths_T_p), 1), dtype=np.int64)
    else:
        im_dp1 = np.zeros((len(paths_T_p), 1), dtype=np.int64)

    # Embed Omega_3(T\v) into A_3(T)
    embed_matrix = np.zeros((len(paths_T_p), len(paths_Tv_p)), dtype=np.int64)
    for j, path_tv in enumerate(paths_Tv_p):
        embedded = embed_path(path_tv)
        if embedded in idx_Tp:
            embed_matrix[idx_Tp[embedded], j] = 1

    embedded_omega_Tv = embed_matrix @ omega_basis_Tv.T % PRIME

    # Z_3(T\v) in Omega coords
    idx_Tvpm1 = {tuple(q): i for i, q in enumerate(paths_Tv_pm1)}
    bd_Tv_p = np.zeros((len(paths_Tv_pm1), len(paths_Tv_p)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_p):
        for sign, face in boundary_faces(sigma):
            face_t = tuple(face)
            if face_t in idx_Tvpm1:
                bd_Tv_p[idx_Tvpm1[face_t], j] = (bd_Tv_p[idx_Tvpm1[face_t], j] + sign) % PRIME

    d_p_Tv = bd_Tv_p @ omega_basis_Tv.T % PRIME
    M_list = [[int(x) % PRIME for x in d_p_Tv[r]] for r in range(d_p_Tv.shape[0])]
    _, zycles_Tv = _gauss_nullbasis_modp(M_list, d_p_Tv.shape[0], d_p_Tv.shape[1], PRIME)
    if not zycles_Tv:
        return 0, bT, bTv
    zycles_Tv = np.array(zycles_Tv, dtype=np.int64)

    embedded_cycles = embedded_omega_Tv @ zycles_Tv.T % PRIME

    rank_B = _gauss_rank_np(im_dp1.copy() % PRIME, PRIME)
    combined = np.hstack([im_dp1, embedded_cycles]) % PRIME
    rank_combined = _gauss_rank_np(combined, PRIME)
    return rank_combined - rank_B, bT, bTv


def main():
    print("=" * 70)
    print("FAST SEARCH FOR beta_3=beta_4=1 TOURNAMENTS (n=8)")
    print("=" * 70)

    n = 8
    max_p = 7

    # Phase 1: Find beta_3=beta_4>0 tournaments quickly
    print(f"\n--- Phase 1: Finding beta_3*beta_4 > 0 tournaments ---")
    t0 = time.time()
    found_tours = []
    total = 0

    # Try many different seeds
    for seed in range(100):
        rng = np.random.RandomState(seed)
        for trial in range(200):
            A = random_tournament(n, rng)
            total += 1
            cc = full_chain_complex_modp(A, n, max_p)
            b3 = cc['bettis'].get(3, 0)
            b4 = cc['bettis'].get(4, 0)
            if b3 > 0 and b4 > 0:
                profile = tuple(cc['bettis'].get(p, 0) for p in range(max_p + 1))
                scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                found_tours.append((A.copy(), cc, profile, scores))
                print(f"  FOUND (seed={seed}, trial={trial}): "
                      f"bettis={profile}, scores={scores}")

        if len(found_tours) >= 5:
            break

        if (seed + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {total} checked, {len(found_tours)} found, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  Phase 1 done: {total} checked, {len(found_tours)} found, {elapsed:.1f}s")

    if not found_tours:
        print("  No beta_3=beta_4>0 found. Try more seeds or larger sample.")
        return

    # Phase 2: Compute rank(i_*^3) for all vertices of found tournaments
    print(f"\n--- Phase 2: rank(i_*^3) analysis ---")
    for idx, (A, cc, profile, scores) in enumerate(found_tours):
        print(f"\n  Tournament #{idx+1}: bettis={profile}, scores={scores}")
        for v in range(n):
            ri3, bT, bTv = compute_istar_p3(A, n, v, max_p)
            bTv_profile = tuple(bTv.get(p, 0) for p in range(max_p + 1))
            b3_sub = bTv.get(3, 0)
            b4_sub = bTv.get(4, 0)
            vtype = "BAD" if b3_sub >= 1 else "GOOD"
            print(f"    v={v}: b(T\\v)={bTv_profile}, rank(i_*^3)={ri3}, {vtype}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
