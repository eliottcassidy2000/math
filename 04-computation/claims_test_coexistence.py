"""
claims_test_coexistence.py — Test Claims I & II on beta_3+beta_4 COEXISTENCE cases

HYP-394 (consecutive seesaw) is REFUTED at n=8 (kind-pasteur-S48).
But Claims I (i_*-injectivity) and II (H_3^rel ≤ 1) may still hold.

This script:
1. Regenerates the 3 confirmed coexistence tournaments from S45 (seed 11111)
2. Tests Claims I and II on each (T, v) pair using mod-p exact computation
3. Also searches for more coexistence cases and tests them

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


def compute_rank_i_star_modp(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T)) using mod-p exact arithmetic."""
    # Build T data
    ap_T = enumerate_all_allowed(A, n, max_p=5)
    paths_3T = ap_T.get(3, [])
    paths_4T = ap_T.get(4, [])
    if not paths_3T:
        return 0, 0, 0

    # Build T\v data
    remaining = [i for i in range(n) if i != v]
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
    n1 = n - 1
    ap_Tv = enumerate_all_allowed(A_sub, n1, max_p=4)
    paths_3Tv = ap_Tv.get(3, [])
    if not paths_3Tv:
        return 0, 0, 0

    # Get Omega_3(T\v) null basis
    P_Tv, na_rows, na_cols = _build_constraint_matrix(ap_Tv, 3, prime)
    if P_Tv is None:
        omega3_Tv_dim = len(paths_3Tv)
        omega3_Tv_basis = np.eye(len(paths_3Tv), dtype=np.int64)
    else:
        rk, nbasis = _gauss_nullbasis_modp(P_Tv, na_rows, na_cols, prime)
        omega3_Tv_dim = na_cols - rk
        if omega3_Tv_dim == 0:
            return 0, 0, 0
        omega3_Tv_basis = np.array(nbasis, dtype=np.int64)

    # Build boundary d3 for T\v in A coordinates
    paths_2Tv = ap_Tv.get(2, [])
    idx2_Tv = {p: i for i, p in enumerate(paths_2Tv)}
    bd3_Tv = np.zeros((len(paths_2Tv), len(paths_3Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3Tv):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % prime

    # ker(d3) in Omega_3(T\v)
    d3_om = bd3_Tv @ omega3_Tv_basis.T % prime
    rk_d3 = _gauss_rank_np(d3_om.copy(), prime)
    ker_d3_dim = omega3_Tv_dim - rk_d3
    if ker_d3_dim == 0:
        return 0, omega3_Tv_dim, rk_d3

    # Build d3 null space in omega coords
    d3_list = [[int(d3_om[i, j]) for j in range(d3_om.shape[1])] for i in range(d3_om.shape[0])]
    _, null_d3 = _gauss_nullbasis_modp(d3_list, d3_om.shape[0], d3_om.shape[1], prime)
    if not null_d3:
        return 0, omega3_Tv_dim, rk_d3
    Z3_Tv_omega = np.array(null_d3, dtype=np.int64)

    # Convert to A_3(T\v) coordinates
    Z3_Tv_Acoords = Z3_Tv_omega @ omega3_Tv_basis % prime

    # Embed T\v paths into T paths
    idx3_T = {p: i for i, p in enumerate(paths_3T)}
    embed = np.zeros((len(paths_3Tv), len(paths_3T)), dtype=np.int64)
    for j, path_tv in enumerate(paths_3Tv):
        embedded = tuple(remaining[k] for k in path_tv)
        if embedded in idx3_T:
            embed[j, idx3_T[embedded]] = 1

    Z3_embedded = Z3_Tv_Acoords @ embed % prime

    # Build im(d4_T) in A_3(T) coordinates
    if not paths_4T:
        im_d4_rank = 0
        rk_Z = _gauss_rank_np(Z3_embedded.copy(), prime)
        return rk_Z, omega3_Tv_dim, rk_d3

    idx3_T_map = {p: i for i, p in enumerate(paths_3T)}
    bd4_T = np.zeros((len(paths_3T), len(paths_4T)), dtype=np.int64)
    for j, path in enumerate(paths_4T):
        for sign, face in boundary_faces(path):
            if face in idx3_T_map:
                bd4_T[idx3_T_map[face], j] = (bd4_T[idx3_T_map[face], j] + sign) % prime

    # Get Omega_4(T) basis
    P_4T, na4_rows, na4_cols = _build_constraint_matrix(ap_T, 4, prime)
    if P_4T is None:
        omega4_T_basis = np.eye(len(paths_4T), dtype=np.int64)
    else:
        rk4, nb4 = _gauss_nullbasis_modp(P_4T, na4_rows, na4_cols, prime)
        if na4_cols - rk4 == 0:
            im_d4_rank = 0
            rk_Z = _gauss_rank_np(Z3_embedded.copy(), prime)
            return rk_Z, omega3_Tv_dim, rk_d3
        omega4_T_basis = np.array(nb4, dtype=np.int64)

    im_d4 = bd4_T @ omega4_T_basis.T % prime
    im_d4_rank = _gauss_rank_np(im_d4.copy(), prime)

    # Concatenate: [im_d4 | Z3_embedded^T]
    combined = np.concatenate([im_d4, Z3_embedded.T], axis=1) % prime
    combined_rank = _gauss_rank_np(combined.copy(), prime)

    rank_i_star = combined_rank - im_d4_rank
    return rank_i_star, omega3_Tv_dim, rk_d3


def test_claims_on_tournament(A, n, label=""):
    """Test Claims I and II on a tournament T."""
    res_T = full_chain_complex_modp(A, n, max_p=7)
    bettis_T = res_T['bettis']
    profile = tuple(bettis_T.get(p, 0) for p in range(n))
    b3_T = bettis_T.get(3, 0)
    b4_T = bettis_T.get(4, 0)
    score = tuple(sorted([int(sum(A[i])) for i in range(n)]))

    print(f"\n  {label}: bettis={profile}, score={score}")

    claim_I_ok = True
    claim_II_ok = True

    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
        res_Tv = full_chain_complex_modp(A_sub, n-1, max_p=6)
        b3_Tv = res_Tv['bettis'].get(3, 0)

        if b3_Tv > 0 and b3_T > 0:
            rank_istar, om_dim, rk_d3 = compute_rank_i_star_modp(A, n, v)
            h3_rel = b3_T - rank_istar
            tag = ""
            if rank_istar != 1:
                tag = " <== CLAIM I VIOLATION!"
                claim_I_ok = False
            elif h3_rel > 1:
                tag = " <== CLAIM II VIOLATION!"
                claim_II_ok = False
            else:
                tag = " OK"
            print(f"    v={v}: b3(T\\v)={b3_Tv}, b4(T)={b4_T}, "
                  f"rank(i_*)={rank_istar}, H_3^rel={h3_rel}{tag}")
        elif b3_Tv == 0 and b3_T > 0:
            h3_rel = b3_T  # rank(i_*) = 0
            tag = " OK" if h3_rel <= 1 else " <== CLAIM II VIOLATION!"
            if h3_rel > 1:
                claim_II_ok = False
            print(f"    v={v}: b3(T\\v)=0, b4(T)={b4_T}, H_3^rel={h3_rel}{tag}")

    return claim_I_ok, claim_II_ok


def main():
    print("=" * 70)
    print("CLAIMS I & II TEST ON BETA_3+BETA_4 COEXISTENCE (mod-p exact)")
    print("=" * 70)

    n = 8

    # Part 1: Test the 3 confirmed coexistence tournaments
    print("\n--- Part 1: 3 confirmed coexistence tournaments (seed 11111) ---")
    rng = np.random.RandomState(11111)
    target_trials = {498, 851, 994}

    coexist_As = []
    for trial in range(max(target_trials) + 1):
        A = random_tournament(n, rng)
        if trial in target_trials:
            coexist_As.append((trial, A.copy()))

    all_I_ok = True
    all_II_ok = True
    for trial, A in coexist_As:
        I_ok, II_ok = test_claims_on_tournament(A, n, label=f"Trial {trial}")
        all_I_ok = all_I_ok and I_ok
        all_II_ok = all_II_ok and II_ok

    print(f"\n  CLAIM I (i_*-injectivity): {'HOLDS' if all_I_ok else 'VIOLATED'}")
    print(f"  CLAIM II (H_3^rel <= 1):   {'HOLDS' if all_II_ok else 'VIOLATED'}")

    # Part 2: Search for more coexistence cases and test them
    print(f"\n--- Part 2: Search for more coexistence + test ---")
    rng2 = np.random.RandomState(54321)
    found = 0
    t0 = time.time()

    for trial in range(10000):
        A = random_tournament(n, rng2)
        res = full_chain_complex_modp(A, n, max_p=5)
        b3 = res['bettis'].get(3, 0)
        if b3 == 0:
            continue

        res_full = full_chain_complex_modp(A, n, max_p=7)
        b4 = res_full['bettis'].get(4, 0)
        if b4 > 0:
            found += 1
            I_ok, II_ok = test_claims_on_tournament(A, n, label=f"New #{found} (trial {trial})")
            all_I_ok = all_I_ok and I_ok
            all_II_ok = all_II_ok and II_ok

        if (trial+1) % 2000 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/10000, {elapsed:.1f}s, found {found} coexistence", flush=True)

    print(f"\n  Part 2: found {found} coexistence cases in 10000 trials")
    print(f"\n{'='*70}")
    print(f"FINAL VERDICT:")
    print(f"  CLAIM I (i_*-injectivity): {'HOLDS' if all_I_ok else 'VIOLATED'}")
    print(f"  CLAIM II (H_3^rel <= 1):   {'HOLDS' if all_II_ok else 'VIOLATED'}")
    print(f"{'='*70}")
    print("\nDONE.")


if __name__ == '__main__':
    main()
