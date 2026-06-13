"""
beta3_2_les_analysis.py — LES analysis of beta_3=2 vs beta_3=1 tournaments at n=8

Key question: Why does beta_3=2 occur? What goes wrong with the proof architecture?

For beta_3=1: Claim I (rank(i_*)=1 for bad vertices) holds at n<=7 but FAILS at n=8
For beta_3=2: Multiple H_3 generators survive — why?

This script studies:
1. For each beta_3=2 tournament, compute the full LES data for each vertex
2. Compare rank(i_*) patterns between beta_3=1 and beta_3=2
3. Check if bad vertices satisfy quasi-iso at n=8
4. Look for the structural mechanism allowing 2 independent H_3 generators

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp,
    enumerate_all_allowed, boundary_faces,
    _build_constraint_matrix, _gauss_nullbasis_modp, _gauss_rank_np,
    RANK_PRIME
)
from fast_beta3 import fast_beta3_nullbasis

PRIME = RANK_PRIME


def compute_istar_rank(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T)) using mod-p exact arithmetic.
    Returns (rank_istar, b3_T, b3_Tv)."""
    # Get b3(T)
    res_T = full_chain_complex_modp(A, n, max_p=5)
    b3_T = res_T['bettis'].get(3, 0)
    if b3_T == 0:
        return 0, 0, 0

    # Get b3(T\v)
    remaining = [i for i in range(n) if i != v]
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
    res_Tv = full_chain_complex_modp(A_sub, n-1, max_p=5)
    b3_Tv = res_Tv['bettis'].get(3, 0)

    if b3_Tv == 0:
        return 0, b3_T, 0

    # Compute rank(i_*) by embedding Z_3(T\v) into Z_3(T) and checking rank mod B_3(T)
    ap_T = enumerate_all_allowed(A, n, max_p=5)
    paths_3T = ap_T.get(3, [])
    paths_4T = ap_T.get(4, [])
    if not paths_3T:
        return 0, b3_T, b3_Tv

    ap_Tv = enumerate_all_allowed(A_sub, n-1, max_p=4)
    paths_3Tv = ap_Tv.get(3, [])
    if not paths_3Tv:
        return 0, b3_T, b3_Tv

    # Get Omega_3(T\v) null basis
    P_Tv, na_rows, na_cols = _build_constraint_matrix(ap_Tv, 3, prime)
    if P_Tv is None:
        omega3_Tv_dim = len(paths_3Tv)
        omega3_Tv_basis = np.eye(len(paths_3Tv), dtype=np.int64)
    else:
        rk, nbasis = _gauss_nullbasis_modp(P_Tv, na_rows, na_cols, prime)
        omega3_Tv_dim = na_cols - rk
        if omega3_Tv_dim == 0:
            return 0, b3_T, b3_Tv
        omega3_Tv_basis = np.array(nbasis, dtype=np.int64)

    # ker(d3) in Omega_3(T\v) — Z_3(T\v)
    paths_2Tv = ap_Tv.get(2, [])
    idx2_Tv = {p: i for i, p in enumerate(paths_2Tv)}
    bd3_Tv = np.zeros((len(paths_2Tv), len(paths_3Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3Tv):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % prime

    d3_om = bd3_Tv @ omega3_Tv_basis.T % prime
    rk_d3 = _gauss_rank_np(d3_om.copy(), prime)
    ker_d3_dim = omega3_Tv_dim - rk_d3
    if ker_d3_dim == 0:
        return 0, b3_T, b3_Tv

    # Build Z_3(T\v) in A_3(T\v) coordinates
    d3_list = [[int(d3_om[i, j]) for j in range(d3_om.shape[1])] for i in range(d3_om.shape[0])]
    _, null_d3 = _gauss_nullbasis_modp(d3_list, d3_om.shape[0], d3_om.shape[1], prime)
    if not null_d3:
        return 0, b3_T, b3_Tv
    Z3_Tv_omega = np.array(null_d3, dtype=np.int64)
    Z3_Tv_Acoords = Z3_Tv_omega @ omega3_Tv_basis % prime

    # Embed T\v paths into T paths
    idx3_T = {p: i for i, p in enumerate(paths_3T)}
    embed = np.zeros((len(paths_3Tv), len(paths_3T)), dtype=np.int64)
    for j, path_tv in enumerate(paths_3Tv):
        embedded = tuple(remaining[k] for k in path_tv)
        if embedded in idx3_T:
            embed[j, idx3_T[embedded]] = 1

    Z3_embedded = Z3_Tv_Acoords @ embed % prime

    # Build im(d_4^T) to quotient by boundaries
    if not paths_4T:
        rk_Z = _gauss_rank_np(Z3_embedded.copy(), prime)
        return rk_Z, b3_T, b3_Tv

    idx3_T_map = {p: i for i, p in enumerate(paths_3T)}
    bd4_T = np.zeros((len(paths_3T), len(paths_4T)), dtype=np.int64)
    for j, path in enumerate(paths_4T):
        for sign, face in boundary_faces(path):
            if face in idx3_T_map:
                bd4_T[idx3_T_map[face], j] = (bd4_T[idx3_T_map[face], j] + sign) % prime

    P_4T, na4_rows, na4_cols = _build_constraint_matrix(ap_T, 4, prime)
    if P_4T is None:
        omega4_T_basis = np.eye(len(paths_4T), dtype=np.int64)
    else:
        rk4, nb4 = _gauss_nullbasis_modp(P_4T, na4_rows, na4_cols, prime)
        if na4_cols - rk4 == 0:
            rk_Z = _gauss_rank_np(Z3_embedded.copy(), prime)
            return rk_Z, b3_T, b3_Tv
        omega4_T_basis = np.array(nb4, dtype=np.int64)

    im_d4 = bd4_T @ omega4_T_basis.T % prime
    im_d4_rank = _gauss_rank_np(im_d4.copy(), prime)

    combined = np.concatenate([im_d4, Z3_embedded.T], axis=1) % prime
    combined_rank = _gauss_rank_np(combined.copy(), prime)

    return combined_rank - im_d4_rank, b3_T, b3_Tv


def analyze_les(A, n, label=""):
    """Full LES analysis for one tournament."""
    print(f"\n  {label}")

    res = full_chain_complex_modp(A, n, max_p=5)
    b3 = res['bettis'].get(3, 0)
    scores = sorted([int(sum(A[i])) for i in range(n)])
    print(f"    b3={b3}, scores={scores}")
    print(f"    omega_dims: {res['omega_dims']}")
    print(f"    ranks: {res['ranks']}")
    print(f"    bettis: {res['bettis']}")

    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
        res_v = full_chain_complex_modp(A_sub, n-1, max_p=5)
        b3_v = res_v['bettis'].get(3, 0)
        b4_v = res_v['bettis'].get(4, 0)
        out_deg = int(sum(A[v]))

        if b3_v > 0 and b3 > 0:
            rank_istar, _, _ = compute_istar_rank(A, n, v)
            h3_rel = b3 - rank_istar
            tag = ""
            if b3 == 1 and rank_istar != 1:
                tag = " <== CLAIM I VIOLATION"
            elif b3 == 2 and rank_istar == 0:
                tag = " <== i_* ZERO"
            print(f"    v={v}(out={out_deg}): b3_Tv={b3_v}, b4_Tv={b4_v}, "
                  f"rank(i_*)={rank_istar}, H3_rel={h3_rel}{tag}")
        elif b3_v == 0:
            h3_rel = b3  # rank(i_*) = 0 trivially
            print(f"    v={v}(out={out_deg}): b3_Tv=0, b4_Tv={b4_v}, "
                  f"H3_rel={h3_rel} [GOOD]")
        else:
            print(f"    v={v}(out={out_deg}): b3_Tv={b3_v}, b4_Tv={b4_v} [b3(T)=0]")


def main():
    print("=" * 70)
    print("BETA_3=2 LES ANALYSIS AT n=8")
    print("=" * 70)

    n = 8

    # Find beta_3=2 and beta_3=1 tournaments
    print("\n--- Finding beta_3=2 and beta_3=1 examples ---")
    rng = np.random.RandomState(12345)

    b3_2_tours = []
    b3_1_tours = []

    t0 = time.time()
    for trial in range(5000):
        A = random_tournament(n, rng)
        b3 = fast_beta3_nullbasis(A, n)

        if b3 == 2 and len(b3_2_tours) < 4:
            b3_2_tours.append((trial, A.copy()))
            print(f"  Found b3=2 at trial {trial}")
        elif b3 == 1 and len(b3_1_tours) < 2:
            b3_1_tours.append((trial, A.copy()))

        if len(b3_2_tours) >= 4 and len(b3_1_tours) >= 2:
            break

        if (trial + 1) % 1000 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/5000, {elapsed:.1f}s, b3=2: {len(b3_2_tours)}, b3=1: {len(b3_1_tours)}")

    print(f"\nFound {len(b3_2_tours)} b3=2 and {len(b3_1_tours)} b3=1 tournaments")

    # Analyze beta_3=2 tournaments
    print(f"\n{'='*60}")
    print("BETA_3=2 TOURNAMENTS")
    print(f"{'='*60}")
    for trial, A in b3_2_tours:
        analyze_les(A, n, label=f"b3=2, trial {trial}")

    # Analyze beta_3=1 tournaments for comparison
    print(f"\n{'='*60}")
    print("BETA_3=1 TOURNAMENTS (comparison)")
    print(f"{'='*60}")
    for trial, A in b3_1_tours:
        analyze_les(A, n, label=f"b3=1, trial {trial}")

    # Summary
    print(f"\n{'='*60}")
    print("KEY QUESTIONS:")
    print("  1. In b3=2 cases, do ALL bad vertices have rank(i_*)=0?")
    print("  2. In b3=2 cases, what is H3_rel for good vertices?")
    print("  3. What distinguishes b3=2 from b3=1 structurally?")
    print(f"{'='*60}")
    print("DONE.")


if __name__ == '__main__':
    main()
