"""
verify_failure.py — Double-check n=8 i_*-injectivity failures

Find an n=8 failure case with safe arithmetic, then verify it using
two completely independent methods:
1. The ψ-containment method (from psi_image_analysis)
2. Direct embedding of H_3(T\v) generator into T and checking mod im(d_4^T)

If both agree, the failure is genuine (not a numerical artifact).

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    matmul_mod,
    RANK_PRIME
)

PRIME = RANK_PRIME


def verify_istar_two_methods(A, n, v):
    """Check rank(i_*) via two independent methods.

    Method 1: ψ-containment (relative cycle boundary check)
    Method 2: Direct embedding of H_3(T\v) gen into T, check mod im(d_4^T)
    """
    max_p = min(n - 1, 6)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1-1))

    paths_T = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    def get_omega(ap, deg):
        ps = ap.get(deg, [])
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    ob4_T = get_omega(ap_T, 4)
    ob3_T = get_omega(ap_T, 3)
    ob4_Tv = get_omega(ap_Tv, 4)
    ob3_Tv = get_omega(ap_Tv, 3)

    idx3_T = {p: i for i, p in enumerate(paths_T.get(3, []))}
    idx3_Tv = {p: i for i, p in enumerate(paths_Tv.get(3, []))}
    idx2_T = {p: i for i, p in enumerate(paths_T.get(2, []))}
    idx2_Tv = {p: i for i, p in enumerate(paths_Tv.get(2, []))}

    # ===== Build boundary maps SAFELY =====
    bd4_T = np.zeros((len(paths_T.get(3, [])), len(paths_T.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_T.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    bd3_T = np.zeros((len(paths_T.get(2, [])), len(paths_T.get(3, []))), dtype=np.int64)
    for j, path in enumerate(paths_T.get(3, [])):
        for sign, face in boundary_faces(path):
            if face in idx2_T:
                bd3_T[idx2_T[face], j] = (bd3_T[idx2_T[face], j] + sign) % PRIME

    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    bd3_Tv = np.zeros((len(paths_Tv.get(2, [])), len(paths_Tv.get(3, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(3, [])):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % PRIME

    # ===== METHOD 2: Direct embedding =====
    # Compute H_3(T\v) generator in A_3(T\v) coords
    # ker(d_3^{T\v}) in Omega_3 coords
    d3o_Tv = matmul_mod(bd3_Tv, ob3_Tv.T)
    if d3o_Tv.size > 0 and d3o_Tv.shape[1] > 0:
        _, kv_Tv = _gauss_nullbasis_modp(
            d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3_Tv = np.array(kv_Tv, dtype=np.int64) if kv_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3_Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3_Tv_A = matmul_mod(ker_d3_Tv, ob3_Tv) if ker_d3_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)

    # im(d_4^{T\v}) in A_3(T\v) coords
    im_d4_Tv_omega = matmul_mod(bd4_Tv, ob4_Tv.T) if ob4_Tv.shape[0] > 0 else np.zeros((len(paths_Tv.get(3, [])), 0), dtype=np.int64)
    im_d4_Tv = im_d4_Tv_omega.T

    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0
    b3_Tv = ker_d3_Tv_A.shape[0] - rk_d4_Tv
    if b3_Tv != 1:
        return None

    # H_3(T\v) representative: pick a cycle not in im(d_4^{T\v})
    # Use: last row of ker_d3_Tv_A (should be independent of im d_4)
    # Actually just take any row that extends the image
    h3_gen_Tv = None
    for row_i in range(ker_d3_Tv_A.shape[0]):
        if im_d4_Tv.shape[0] > 0:
            test = np.vstack([im_d4_Tv, ker_d3_Tv_A[row_i:row_i+1]]) % PRIME
            rk_test = int(_gauss_rank_np(test.copy(), PRIME))
            if rk_test > rk_d4_Tv:
                h3_gen_Tv = ker_d3_Tv_A[row_i]
                break
        else:
            if np.any(ker_d3_Tv_A[row_i] % PRIME != 0):
                h3_gen_Tv = ker_d3_Tv_A[row_i]
                break

    if h3_gen_Tv is None:
        return None

    # Embed into A_3(T) coords
    embed_3 = np.zeros((len(paths_Tv.get(3, [])), len(paths_T.get(3, []))), dtype=np.int64)
    for jv, pv in enumerate(paths_Tv.get(3, [])):
        mapped = tuple(remaining[x] for x in pv)
        if mapped in idx3_T:
            embed_3[jv, idx3_T[mapped]] = 1

    h3_gen_T = (h3_gen_Tv @ embed_3) % PRIME  # in A_3(T) coords

    # im(d_4^T) in A_3(T) coords
    im_d4_T_omega = matmul_mod(bd4_T, ob4_T.T)
    im_d4_T = im_d4_T_omega.T
    rk_d4_T = int(_gauss_rank_np(im_d4_T.copy(), PRIME)) if im_d4_T.shape[0] > 0 else 0

    # ker(d_3^T) in Omega_3 coords
    d3o_T = matmul_mod(bd3_T, ob3_T.T)
    if d3o_T.size > 0 and d3o_T.shape[1] > 0:
        _, kv_T = _gauss_nullbasis_modp(
            d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME)
        ker_d3_T = np.array(kv_T, dtype=np.int64) if kv_T else np.zeros((0, ob3_T.shape[0]), dtype=np.int64)
    else:
        ker_d3_T = np.eye(ob3_T.shape[0], dtype=np.int64)

    ker_d3_T_A = matmul_mod(ker_d3_T, ob3_T) if ker_d3_T.shape[0] > 0 else np.zeros((0, len(paths_T.get(3, []))), dtype=np.int64)
    b3_T = ker_d3_T_A.shape[0] - rk_d4_T
    if b3_T != 1:
        return None

    # Check: is h3_gen_T a cycle of T? (it should be since it's embedded from a T\v cycle)
    check_cycle = matmul_mod(bd3_T, h3_gen_T.reshape(-1, 1)).flatten() % PRIME
    is_cycle = np.all(check_cycle == 0)

    # Check: is h3_gen_T in im(d_4^T)?
    if im_d4_T.shape[0] > 0:
        test2 = np.vstack([im_d4_T, h3_gen_T.reshape(1, -1)]) % PRIME
        rk_test2 = int(_gauss_rank_np(test2.copy(), PRIME))
        in_im_d4_T = (rk_test2 == rk_d4_T)
    else:
        in_im_d4_T = np.all(h3_gen_T % PRIME == 0)

    # METHOD 2 result: rank(i_*) = 1 iff h3_gen_T is NOT in im(d_4^T)
    method2_rank = 0 if in_im_d4_T else 1

    # ===== METHOD 1: ψ-containment =====
    tv4 = [i for i, p in enumerate(paths_T.get(4, [])) if v in p]
    old3 = [i for i, p in enumerate(paths_T.get(3, [])) if v not in p]
    tv3 = [i for i, p in enumerate(paths_T.get(3, [])) if v in p]

    if not tv4:
        return None

    d4_tvtv = matmul_mod(bd4_T[tv3, :], ob4_T.T)
    if d4_tvtv.size > 0 and d4_tvtv.shape[1] > 0:
        _, kv = _gauss_nullbasis_modp(
            d4_tvtv.astype(np.int32), d4_tvtv.shape[0], d4_tvtv.shape[1], PRIME)
        ker_tvtv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
    else:
        ker_tvtv = np.eye(ob4_T.shape[0], dtype=np.int64)

    if ker_tvtv.shape[0] > 0:
        ker_A4 = matmul_mod(ker_tvtv, ob4_T)
        psi_full = matmul_mod(bd4_T[old3, :], ker_A4.T).T
    else:
        psi_full = np.zeros((0, len(old3)), dtype=np.int64)

    # Map psi_full to T\v coords
    psi_Tv = np.zeros((psi_full.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
    for local_j, global_i in enumerate(old3):
        path_T = paths_T[3][global_i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            j_Tv = idx3_Tv[path_Tv]
            psi_Tv[:, j_Tv] = (psi_Tv[:, j_Tv] + psi_full[:, local_j]) % PRIME

    # Check containment
    if psi_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        method1_rank = 0 if rk_combined > rk_d4_Tv else 1
    elif psi_Tv.shape[0] > 0:
        rk_psi = int(_gauss_rank_np(psi_Tv.copy(), PRIME))
        method1_rank = 0 if rk_psi > 0 else 1
    else:
        method1_rank = 1

    return {
        'method1_rank': method1_rank,
        'method2_rank': method2_rank,
        'is_cycle': is_cycle,
        'b3_T': b3_T,
        'b3_Tv': b3_Tv,
        'agree': method1_rank == method2_rank,
    }


def main():
    print("VERIFY n=8 FAILURES WITH TWO METHODS")
    print("="*70)

    rng = np.random.RandomState(42)
    failures_verified = 0
    total = 0
    agree_count = 0
    t0 = time.time()

    for trial in range(80000):
        if failures_verified >= 10 and total >= 500:
            break
        A = random_tournament(8, rng)
        cc = full_chain_complex_modp(A, 8, 7)
        if cc['bettis'].get(3, 0) != 1:
            continue

        for v in range(8):
            r = verify_istar_two_methods(A, 8, v)
            if r is None:
                continue

            total += 1
            if r['agree']:
                agree_count += 1

            if r['method1_rank'] == 0 or r['method2_rank'] == 0:
                failures_verified += 1
                print(f"\n  FAILURE #{failures_verified}: method1={r['method1_rank']}, method2={r['method2_rank']}, "
                      f"agree={r['agree']}, is_cycle={r['is_cycle']}, b3_T={r['b3_T']}, b3_Tv={r['b3_Tv']}")

            if total % 100 == 0:
                print(f"  ... {total} checked, {failures_verified} failures, {agree_count}/{total} agree")

            if total >= 500 and failures_verified >= 10:
                break

    elapsed = time.time() - t0
    print(f"\n  TOTAL: {total} (T,v) pairs, {failures_verified} failures, {elapsed:.1f}s")
    print(f"  Methods agree: {agree_count}/{total}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
