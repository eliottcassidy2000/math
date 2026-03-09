"""
claims_n8_extended.py — Extended search for Claim I/II violations at n=8

From claims_test_coexistence.py Part 1:
  - Trial 498: Claims hold (all i_*=1 for bad, H_3^rel<=1)
  - Trial 851: Claims hold
  - Trial 994: CLAIM I VIOLATION at v=1: rank(i_*)=0 when b3(T)=1 & b3(T\v)=1
    BUT Claim II still holds: H_3^rel=1 <= 1

This script searches more broadly:
1. Tests Claim I and II on non-coexistence beta_3=1 tours at n=8
2. Searches for more coexistence cases
3. Tests good vertex existence (which is what actually proves beta_3<=1)

Author: kind-pasteur-S48 (2026-03-09)
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

PRIME = RANK_PRIME


def compute_rank_i_star(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T))."""
    ap_T = enumerate_all_allowed(A, n, max_p=5)
    paths_3T = ap_T.get(3, [])
    paths_4T = ap_T.get(4, [])
    if not paths_3T:
        return 0

    remaining = [i for i in range(n) if i != v]
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
    ap_Tv = enumerate_all_allowed(A_sub, n-1, max_p=4)
    paths_3Tv = ap_Tv.get(3, [])
    if not paths_3Tv:
        return 0

    P_Tv, na_rows, na_cols = _build_constraint_matrix(ap_Tv, 3, prime)
    if P_Tv is None:
        omega3_Tv_dim = len(paths_3Tv)
        omega3_Tv_basis = np.eye(len(paths_3Tv), dtype=np.int64)
    else:
        rk, nbasis = _gauss_nullbasis_modp(P_Tv, na_rows, na_cols, prime)
        omega3_Tv_dim = na_cols - rk
        if omega3_Tv_dim == 0:
            return 0
        omega3_Tv_basis = np.array(nbasis, dtype=np.int64)

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
        return 0

    d3_list = [[int(d3_om[i, j]) for j in range(d3_om.shape[1])] for i in range(d3_om.shape[0])]
    _, null_d3 = _gauss_nullbasis_modp(d3_list, d3_om.shape[0], d3_om.shape[1], prime)
    if not null_d3:
        return 0
    Z3_Tv_omega = np.array(null_d3, dtype=np.int64)
    Z3_Tv_Acoords = Z3_Tv_omega @ omega3_Tv_basis % prime

    idx3_T = {p: i for i, p in enumerate(paths_3T)}
    embed = np.zeros((len(paths_3Tv), len(paths_3T)), dtype=np.int64)
    for j, path_tv in enumerate(paths_3Tv):
        embedded = tuple(remaining[k] for k in path_tv)
        if embedded in idx3_T:
            embed[j, idx3_T[embedded]] = 1

    Z3_embedded = Z3_Tv_Acoords @ embed % prime

    if not paths_4T:
        rk_Z = _gauss_rank_np(Z3_embedded.copy(), prime)
        return rk_Z

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
            return rk_Z
        omega4_T_basis = np.array(nb4, dtype=np.int64)

    im_d4 = bd4_T @ omega4_T_basis.T % prime
    im_d4_rank = _gauss_rank_np(im_d4.copy(), prime)

    combined = np.concatenate([im_d4, Z3_embedded.T], axis=1) % prime
    combined_rank = _gauss_rank_np(combined.copy(), prime)

    return combined_rank - im_d4_rank


def main():
    print("=" * 70)
    print("EXTENDED CLAIMS TEST AT n=8 (mod-p exact)")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(12345)

    claim_I_violations = 0
    claim_II_violations = 0
    good_vertex_missing = 0
    total_b3_tours = 0
    total_coexist = 0
    total_checked = 0

    t0 = time.time()
    N_TRIALS = 5000

    for trial in range(N_TRIALS):
        A = random_tournament(n, rng)
        total_checked += 1

        # Quick b3 check
        res = full_chain_complex_modp(A, n, max_p=5)
        b3 = res['bettis'].get(3, 0)
        if b3 == 0:
            continue

        total_b3_tours += 1

        # Full Betti to check for coexistence
        res_full = full_chain_complex_modp(A, n, max_p=7)
        b4 = res_full['bettis'].get(4, 0)
        if b4 > 0:
            total_coexist += 1

        # Test good vertex existence and Claims I/II
        has_good = False
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            res_Tv = full_chain_complex_modp(A_sub, n-1, max_p=6)
            b3_Tv = res_Tv['bettis'].get(3, 0)

            if b3_Tv == 0:
                has_good = True
                # Good vertex: H_3^rel = b3(T) = 1
                if b3 > 1:
                    claim_II_violations += 1
            elif b3_Tv > 0:
                rank_istar = compute_rank_i_star(A, n, v)
                h3_rel = b3 - rank_istar
                if rank_istar != 1:
                    claim_I_violations += 1
                    score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                    coexist_tag = " [COEXIST]" if b4 > 0 else ""
                    print(f"  CLAIM I VIOLATION: trial={trial}, v={v}, "
                          f"b3={b3}, b3_Tv={b3_Tv}, rank(i_*)={rank_istar}, "
                          f"b4={b4}{coexist_tag}, score={score}")
                if h3_rel > 1:
                    claim_II_violations += 1
                    print(f"  CLAIM II VIOLATION: trial={trial}, v={v}, H_3^rel={h3_rel}")

        if not has_good:
            good_vertex_missing += 1
            score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  NO GOOD VERTEX: trial={trial}, b3={b3}, b4={b4}, score={score}")

        if (trial+1) % 500 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{N_TRIALS}, {elapsed:.1f}s, b3_tours={total_b3_tours}, "
                  f"coexist={total_coexist}, CI_viol={claim_I_violations}, "
                  f"CII_viol={claim_II_violations}, no_good={good_vertex_missing}", flush=True)

    t1 = time.time()
    print(f"\n{'='*70}")
    print(f"RESULTS ({total_checked} tournaments, {total_b3_tours} with b3>0, {t1-t0:.1f}s)")
    print(f"  Coexistence (b3>0 & b4>0): {total_coexist}")
    print(f"  Claim I violations (i_* not injective): {claim_I_violations}")
    print(f"  Claim II violations (H_3^rel > 1): {claim_II_violations}")
    print(f"  Missing good vertex: {good_vertex_missing}")
    print(f"\n  VERDICT:")
    print(f"    Claim I: {'HOLDS' if claim_I_violations == 0 else 'VIOLATED'}")
    print(f"    Claim II: {'HOLDS' if claim_II_violations == 0 else 'VIOLATED'}")
    print(f"    Good vertex always exists: {'YES' if good_vertex_missing == 0 else 'NO'}")
    if claim_II_violations == 0 and good_vertex_missing == 0:
        print(f"    => beta_3 <= 1 follows from good-vertex + Claim II at n=8")
    print("DONE.")


if __name__ == '__main__':
    main()
