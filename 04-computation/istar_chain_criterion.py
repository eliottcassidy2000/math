"""
istar_chain_criterion.py — Chain-level criterion for rank(i_*)

From the block decomposition of d_4:
  d_4^{old→old}: Ω_4^old → A_3^old
  d_4^{tv→old}: Ω_4^tv → A_3^old   (v-deletion face map ψ)
  d_4^{tv→tv}:  Ω_4^tv → A_3^tv
  d_4^{old→tv} = 0                   (proved S58)

The embedded H_3(T\v) gen z_0 ∈ im(d_4^T) iff:
  z_0 = d_4^{old→old}(w_old) + ψ(w_tv)  with d_4^{tv→tv}(w_tv) = 0

Since z_0 ∉ im(d_4^{T\v}) = im(d_4^{old→old}|_{Ω_4^{T\v}}):
  rank(i_*) = 0  iff  ψ(ker(d_4^{tv→tv})) reaches the H_3(T\v) direction

Concretely: map ψ(ker(d_4^{tv→tv})) → H_3(T\v) = ker(d_3^{T\v})/im(d_4^{T\v}).
If this map is nonzero, rank(i_*) = 0.

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, deletion_adj,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def chain_criterion(A, n, v):
    """Test the chain-level criterion for rank(i_*)."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_4_Tv = ap_Tv.get(4, [])
    paths_2_Tv = ap_Tv.get(2, [])

    if not paths_3_T or not paths_4_T:
        return None

    tv3 = [i for i, p in enumerate(paths_3_T) if v in p]
    old3 = [i for i, p in enumerate(paths_3_T) if v not in p]
    tv4 = [i for i, p in enumerate(paths_4_T) if v in p]
    old4 = [i for i, p in enumerate(paths_4_T) if v not in p]

    if not tv3 or not old3 or not tv4:
        return None

    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob4_T = get_omega(ap_T, 4)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)

    # Build d_4^T in A-path coords
    idx3_T = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    # Block: d_4^{tv→tv} = bd4_T[tv3, tv4]
    d4_tv_tv = bd4_T[np.ix_(tv3, tv4)] % PRIME
    # Block: ψ = d_4^{tv→old} = bd4_T[old3, tv4]
    psi = bd4_T[np.ix_(old3, tv4)] % PRIME

    # Omega_4^tv: project ob4_T to tv4 and old4 indices
    ob4_tv = ob4_T[:, tv4] % PRIME if tv4 else np.zeros((ob4_T.shape[0], 0), dtype=np.int64)
    ob4_old = ob4_T[:, old4] % PRIME if old4 else np.zeros((ob4_T.shape[0], 0), dtype=np.int64)

    # d_4^{tv→tv} on Omega_4: d4_tv_tv @ ob4_tv^T
    # But we need tv-only Omega_4 vectors too... actually no.
    # ker(d_4^{tv→tv}) is on ALL Omega_4 vectors (not just tv-only).
    # d_4^{tv→tv}(w) picks the tv component of d_4(w). For an old w, d_4^{tv→tv}(w) = 0 (old→tv = 0).
    # So ker(d_4^{tv→tv}) = Omega_4^old + ker(d_4^{tv→tv}|_{Omega_4^tv}).

    # Actually: we work with Omega_4 vectors. d_4^{tv→tv}(w) = d4_tv_tv @ w[tv4].
    # In Omega coords: d_4^{tv→tv} on Omega basis = d4_tv_tv @ ob4_tv^T.
    d4_tvtv_omega = d4_tv_tv @ ob4_tv.T % PRIME  # shape: (len(tv3), dim_omega4)

    # ker(d_4^{tv→tv} on Omega) = null space of d4_tvtv_omega
    if d4_tvtv_omega.size > 0:
        _, kv = _gauss_nullbasis_modp(
            d4_tvtv_omega.astype(np.int32), d4_tvtv_omega.shape[0], d4_tvtv_omega.shape[1], PRIME)
        ker_d4_tvtv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
    else:
        ker_d4_tvtv = np.eye(ob4_T.shape[0], dtype=np.int64)

    dim_ker_d4_tvtv = ker_d4_tvtv.shape[0]

    # ψ(ker(d4^{tv→tv})): map these Omega vectors through ψ
    # First, get the A_4 coords of these kernel vectors
    ker_A4 = ker_d4_tvtv @ ob4_T % PRIME  # shape: (dim_ker, len(paths_4_T))
    # Then apply ψ: old3 components of d_4(ker_vectors)
    # ψ(w) = bd4_T[old3, :] @ w^T for each w
    psi_of_ker = (bd4_T[old3, :] @ ker_A4.T % PRIME).T % PRIME  # shape: (dim_ker, len(old3))

    # Map to T\v coords
    idx3_Tv = {p: i for i, p in enumerate(paths_3_Tv)}
    old3_T_to_Tv = {}
    for i in old3:
        path_T = paths_3_T[i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_T_to_Tv[i] = idx3_Tv[path_Tv]

    psi_ker_Tv = np.zeros((psi_of_ker.shape[0], len(paths_3_Tv)), dtype=np.int64)
    for local_j, global_i in enumerate(old3):
        if global_i in old3_T_to_Tv:
            j_Tv = old3_T_to_Tv[global_i]
            psi_ker_Tv[:, j_Tv] = (psi_ker_Tv[:, j_Tv] + psi_of_ker[:, local_j]) % PRIME

    rk_psi_ker = int(_gauss_rank_np(psi_ker_Tv.copy(), PRIME)) if psi_ker_Tv.shape[0] > 0 else 0

    # Now compute H_3(T\v): ker(d_3^{T\v})/im(d_4^{T\v})
    idx2_Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2_Tv:
                bd3_Tv[idx2_Tv[face], j] = (bd3_Tv[idx2_Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3_Tv = np.array(kv_Tv, dtype=np.int64) if kv_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3_Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3_Tv_A = ker_d3_Tv @ ob3_Tv % PRIME

    bd4_Tv = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0
    beta3_Tv = ker_d3_Tv_A.shape[0] - rk_d4_Tv

    if beta3_Tv != 1:
        return None  # only BAD vertices

    # Does ψ(ker(d4^{tv→tv})) reach H_3(T\v)?
    # Check: ψ_ker_Tv combined with im(d_4^{T\v}) — does rank increase?
    if psi_ker_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_ker_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        psi_reaches_h3 = (rk_combined > rk_d4_Tv)
        extra = rk_combined - rk_d4_Tv
    elif psi_ker_Tv.shape[0] > 0:
        rk_combined = int(_gauss_rank_np(psi_ker_Tv.copy(), PRIME))
        psi_reaches_h3 = (rk_combined > 0)
        extra = rk_combined
    else:
        psi_reaches_h3 = False
        extra = 0

    # Also: does ψ(ker(d4^{tv→tv})) lie in ker(d_3^{T\v})?
    # It SHOULD, because d_3 ∘ d_4 = 0 in the full complex
    # ψ(w) = old part of d_4(w), and d_3(d_4(w)) = 0, so d_3(ψ(w) + tv_part) = 0
    # d_3^{old→old}(ψ(w)) + d_3 applied to tv part = 0
    # But this doesn't immediately put ψ(w) in ker(d_3^{T\v})...
    # Actually ψ(w) is in A_3^old ≅ A_3(T\v). Is d_3^{T\v}(ψ(w)) = 0?
    psi_in_ker = False
    if psi_ker_Tv.shape[0] > 0:
        residual = bd3_Tv @ psi_ker_Tv.T % PRIME
        rk_res = int(_gauss_rank_np(residual.copy(), PRIME))
        psi_in_ker = (rk_res == 0)

    # Also compute rank(i_*) directly for comparison
    # embed ker(d_3^{T\v}) gen into A_3(T)
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3_T:
            embed[jv, idx3_T[mapped]] = 1

    embedded_ker = ker_d3_Tv_A @ embed % PRIME
    # Full im(d_4^T)
    im_d4_T = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4_T = int(_gauss_rank_np(im_d4_T.copy(), PRIME)) if im_d4_T.shape[0] > 0 else 0

    if im_d4_T.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined_full = np.vstack([im_d4_T, embedded_ker]) % PRIME
        rk_comb_full = int(_gauss_rank_np(combined_full.copy(), PRIME))
        rank_istar = rk_comb_full - rk_d4_T
    else:
        rank_istar = 0

    return {
        'dim_ker_d4_tvtv': dim_ker_d4_tvtv,
        'rk_psi_ker': rk_psi_ker,
        'psi_reaches_h3': psi_reaches_h3,
        'extra_directions': extra,
        'psi_in_ker_d3Tv': psi_in_ker,
        'rank_istar': rank_istar,
        'criterion_matches': psi_reaches_h3 == (rank_istar == 0),
        'beta3_Tv': beta3_Tv,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"CHAIN-LEVEL CRITERION FOR rank(i_*) AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(12345)  # use seed that catches failures at n=8
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 800  # larger sample at n=8 to catch failures

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = chain_criterion(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} BAD (T,v) pairs, {elapsed:.1f}s")

        # ψ reaches H_3(T\v)?
        reaches = sum(1 for r in results if r['psi_reaches_h3'])
        print(f"  ψ(ker(d4^tvtv)) reaches H_3(T\\v): {reaches}/{len(results)}")

        # rank(i_*) = 0
        rk0 = sum(1 for r in results if r['rank_istar'] == 0)
        print(f"  rank(i_*) = 0: {rk0}/{len(results)}")

        # Criterion matches?
        matches = sum(1 for r in results if r['criterion_matches'])
        print(f"  Criterion matches: {matches}/{len(results)}")

        # ψ(ker) in ker(d_3^{T\v})?
        in_ker = sum(1 for r in results if r['psi_in_ker_d3Tv'])
        print(f"  ψ(ker) ⊂ ker(d_3^Tv): {in_ker}/{len(results)}")

        # Dimensions
        for key in ['dim_ker_d4_tvtv', 'rk_psi_ker', 'extra_directions']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        if rk0 > 0:
            print(f"\n  --- rank(i_*)=0 cases ---")
            fail_cases = [r for r in results if r['rank_istar'] == 0]
            for r in fail_cases[:5]:
                print(f"    psi_reaches={r['psi_reaches_h3']}, extra={r['extra_directions']}, "
                      f"dim_ker_tvtv={r['dim_ker_d4_tvtv']}, rk_psi={r['rk_psi_ker']}, "
                      f"psi_in_ker={r['psi_in_ker_d3Tv']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
