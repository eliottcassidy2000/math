"""
beta34_coexistence_test.py — Find and analyze beta_3+beta_4 coexistence at n=8

From beta_exotic_patterns_n8.md: profile [1,0,0,1,1,0,0,0] occurs at ~0.15% rate.
Key question: Do Claims I and II of THM-110 still hold when beta_4 > 0?

Claim I: rank(i_*) = 1 when b3(T\v) = 1 and b3(T) >= 1
Claim II: dim H_3(T,T\v) <= 1 for all (T,v)

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp, compute_betti_hybrid,
    enumerate_all_allowed, boundary_faces,
    _build_constraint_matrix, _gauss_nullbasis_modp, _gauss_rank_np,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_rank_i_star(A, n, v, prime=PRIME):
    """Compute rank(i_*: H_3(T\\v) -> H_3(T)).

    Uses the formula: rank(i_*) = rank([im(d4_T) | Z3_Tv_embedded]) - rank(im(d4_T))
    where Z3_Tv = ker(d3) in Omega_3(T\\v) embedded in A_3(T).
    """
    # Build T data
    ap_T = enumerate_all_allowed(A, n, max_p=5)
    paths_3T = ap_T.get(3, [])
    paths_4T = ap_T.get(4, [])
    if not paths_3T:
        return 0, 0, 0

    # Build T\\v data
    remaining = [i for i in range(n) if i != v]
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
    n1 = n - 1
    ap_Tv = enumerate_all_allowed(A_sub, n1, max_p=4)
    paths_3Tv = ap_Tv.get(3, [])
    if not paths_3Tv:
        return 0, 0, 0

    # Get Omega_3(T\\v) null basis
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

    # Build boundary d3 for T\\v in Omega coordinates
    paths_2Tv = ap_Tv.get(2, [])
    idx2_Tv = {p: i for i, p in enumerate(paths_2Tv)}
    bd3_Tv = np.zeros((len(paths_2Tv), len(paths_3Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3Tv):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % prime

    # ker(d3) in Omega_3(T\\v)
    d3_om = bd3_Tv @ omega3_Tv_basis.T % prime
    rk_d3 = _gauss_rank_np(d3_om.copy(), prime)
    ker_d3_dim = omega3_Tv_dim - rk_d3
    if ker_d3_dim == 0:
        return 0, omega3_Tv_dim, rk_d3

    # Build d3 null space in omega coords
    # First get RREF of d3_om
    d3_list = [[int(d3_om[i, j]) for j in range(d3_om.shape[1])] for i in range(d3_om.shape[0])]
    _, null_d3 = _gauss_nullbasis_modp(d3_list, d3_om.shape[0], d3_om.shape[1], prime)
    if not null_d3:
        return 0, omega3_Tv_dim, rk_d3
    Z3_Tv_omega = np.array(null_d3, dtype=np.int64)  # (ker_dim x omega3_Tv_dim)

    # Convert to A_3(T\\v) coordinates
    Z3_Tv_Acoords = Z3_Tv_omega @ omega3_Tv_basis % prime  # (ker_dim x |A_3(T\\v)|)

    # Embed T\\v paths into T paths
    idx3_T = {p: i for i, p in enumerate(paths_3T)}
    embed = np.zeros((len(paths_3Tv), len(paths_3T)), dtype=np.int64)
    for j, path_tv in enumerate(paths_3Tv):
        embedded = tuple(remaining[k] for k in path_tv)
        if embedded in idx3_T:
            embed[j, idx3_T[embedded]] = 1

    Z3_embedded = Z3_Tv_Acoords @ embed % prime  # (ker_dim x |A_3(T)|)

    # Build im(d4_T) in A_3(T) coordinates
    if not paths_4T:
        im_d4_rank = 0
        # rank(i_*) = rank(Z3_embedded)
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

    im_d4 = bd4_T @ omega4_T_basis.T % prime  # (|A_3(T)| x omega4_dim)
    im_d4_rank = _gauss_rank_np(im_d4.copy(), prime)

    # Concatenate: [im_d4 | Z3_embedded^T]
    combined = np.concatenate([im_d4, Z3_embedded.T], axis=1) % prime
    combined_rank = _gauss_rank_np(combined.copy(), prime)

    rank_i_star = combined_rank - im_d4_rank
    return rank_i_star, omega3_Tv_dim, rk_d3


def main():
    print("=" * 70)
    print("BETA_3 + BETA_4 COEXISTENCE: LES CLAIMS TEST")
    print("=" * 70)

    # Search for coexistence at n=8
    print("\n--- Searching for beta_3+beta_4 coexistence at n=8 ---")
    n = 8
    rng = np.random.RandomState(999)  # Different seed from other scripts

    coexist_tours = []
    b3_only_tours = []
    total = 0

    t0 = time.time()
    for trial in range(20000):
        A = random_tournament(n, rng)
        total += 1

        # Fast beta_3 check
        b3 = compute_betti_hybrid(A, n, 3, max_p=5)
        if b3 == 0:
            continue

        # Full check for beta_4
        b4 = compute_betti_hybrid(A, n, 4, max_p=6)
        if b4 > 0:
            coexist_tours.append(A.copy())
            score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            print(f"  FOUND coexistence at trial={trial}: b3={b3}, b4={b4}, score={score}")
        else:
            if len(b3_only_tours) < 20:
                b3_only_tours.append(A.copy())

        if (trial+1) % 2000 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/20000, {elapsed:.1f}s, found {len(coexist_tours)} coexistence",
                  flush=True)

    t1 = time.time()
    print(f"\n  Searched {total} tournaments in {t1-t0:.1f}s")
    print(f"  Found {len(coexist_tours)} with beta_3+beta_4 coexistence")

    if not coexist_tours:
        print("  No coexistence found. Testing claims on beta_3-only tournaments...")
        test_tours = b3_only_tours[:10]
        label = "beta_3=1, beta_4=0"
    else:
        test_tours = coexist_tours[:5] + b3_only_tours[:5]
        label = "mixed"

    # Test Claims I and II
    print(f"\n--- Testing Claims I and II on {label} tournaments ---")

    for idx, A in enumerate(test_tours):
        res = full_chain_complex_modp(A, n, max_p=7)
        b = res['bettis']
        bv = tuple(b.get(p, 0) for p in range(8))
        score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        print(f"\n  Tour {idx}: bettis={bv}, score={score}")

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            data_Tv = full_chain_complex_modp(A_sub, n-1, max_p=6)
            b3_Tv = data_Tv['bettis'].get(3, 0)

            if b3_Tv > 0:
                rank_istar, _, _ = compute_rank_i_star(A, n, v)
                h3_rel = b.get(3, 0) - rank_istar
                print(f"    v={v}: b3(T\\v)={b3_Tv}, rank(i_*)={rank_istar}, H_3^rel={h3_rel}", end="")
                if rank_istar != 1:
                    print(" <== CLAIM I VIOLATION!")
                elif h3_rel > 1:
                    print(" <== CLAIM II VIOLATION!")
                else:
                    print(" OK")
            elif b3_Tv == 0:
                # For good vertices: H_3^rel should = beta_3(T)
                h3_rel = b.get(3, 0)  # since rank(i_*) = 0
                if h3_rel > 1:
                    print(f"    v={v}: b3(T\\v)=0, H_3^rel={h3_rel} <== CLAIM II VIOLATION!")

    print("\nDONE.")


if __name__ == '__main__':
    main()
