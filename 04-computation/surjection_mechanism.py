"""
surjection_mechanism.py — Why does ψ(ker(d_4^{tv→tv})) = im(d_4^{T\v})?

At n=7 (b3=b3_Tv=1), we have the identity:
  ψ(ker(d_4^{tv→tv})) = im(d_4^{T\v})

where ψ = d_4^T|_{old3-rows} applied to ker(d_4^{tv→tv}) vectors.

This means: for every w_old ∈ Ω_4(T\v), d_4^{T\v}(w_old) can be written as
the old3-projection of d_4^T(w) for some w ∈ Ω_4(T) with d_4^{tv→tv}(w) = 0.

Key question: is there a STRUCTURAL chain-map reason this holds?

Hypothesis A: the map d_4^{old→old} already surjects onto im(d_4^{T\v}),
and ker(d_4^{tv→tv}) includes enough "old-only" vectors.

Hypothesis B: there's a functorial lifting — every T\v boundary lifts
to a relative cycle boundary in T.

Let's test both.

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


def test_surjection(A, n, v):
    """Dissect why ψ(ker(d_4^{tv→tv})) = im(d_4^{T\\v})."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    remaining_inv = {remaining[i]: i for i in range(n1)}
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    def classify(ps, v):
        tv = [i for i, p in enumerate(ps) if v in p]
        old = [i for i, p in enumerate(ps) if v not in p]
        return tv, old

    tv4, old4 = classify(paths.get(4, []), v)
    tv3, old3 = classify(paths.get(3, []), v)

    if not tv4 or not old3:
        return None

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
    ob4_Tv = get_omega(ap_Tv, 4)
    ob3_Tv = get_omega(ap_Tv, 3)

    idx3_T = {p: i for i, p in enumerate(paths.get(3, []))}
    idx4_T = {p: i for i, p in enumerate(paths.get(4, []))}
    idx3_Tv = {p: i for i, p in enumerate(paths_Tv.get(3, []))}

    bd4_T = np.zeros((len(paths.get(3, [])), len(paths.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_T:
                bd4_T[idx3_T[face], j] = (bd4_T[idx3_T[face], j] + sign) % PRIME

    # Block decomposition
    d4_oo = bd4_T[np.ix_(old3, old4)] % PRIME  # old→old
    d4_to = bd4_T[np.ix_(old3, tv4)] % PRIME   # tv→old (= ψ)
    d4_tt = bd4_T[np.ix_(tv3, tv4)] % PRIME    # tv→tv

    # Omega_4(T) in A_4(T) coords
    # d_4^{tv→tv} through Omega: d4_tt @ ob4_T[:, tv4].T
    d4tt_omega = d4_tt @ ob4_T[:, tv4].T % PRIME

    # ker(d_4^{tv→tv})
    if d4tt_omega.size > 0 and d4tt_omega.shape[1] > 0:
        _, kv = _gauss_nullbasis_modp(
            d4tt_omega.astype(np.int32), d4tt_omega.shape[0], d4tt_omega.shape[1], PRIME)
        ker_tvtv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob4_T.shape[0]), dtype=np.int64)
    else:
        ker_tvtv = np.eye(ob4_T.shape[0], dtype=np.int64) if ob4_T.shape[0] > 0 else np.zeros((0, 0), dtype=np.int64)

    ker_A4 = ker_tvtv @ ob4_T % PRIME if ker_tvtv.shape[0] > 0 else np.zeros((0, len(paths.get(4, []))), dtype=np.int64)

    # ψ(ker): full boundary restricted to old3
    psi_ker = (bd4_T[old3, :] @ ker_A4.T % PRIME).T if ker_A4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)

    # TEST A: d_4^{old→old} alone — does it already surject onto im(d_4^{T\v})?
    # d_4^{old→old} in Omega coords: d4_oo @ ob4_T[:, old4].T
    # But wait: ob4_T has mixed support. Need the old4-part of each Omega basis vector.
    # Some Omega vectors may have BOTH tv4 and old4 components.
    # The "old-only" Omega vectors are those supported entirely on old4.
    #
    # Actually: ker_A4 has BOTH old4 and tv4 components.
    # psi_ker = d4_oo @ ker_A4[:, old4].T + d4_to @ ker_A4[:, tv4].T
    # = (old→old contribution) + (tv→old contribution)

    # Contribution from old→old block only
    oo_contrib = (d4_oo @ ker_A4[:, old4].T % PRIME).T if ker_A4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)
    # Contribution from tv→old block only
    to_contrib = (d4_to @ ker_A4[:, tv4].T % PRIME).T if ker_A4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)

    # Map everything to T\v coords for comparison
    old3_to_Tv = {}
    for local_j, global_i in enumerate(old3):
        path_T = paths[3][global_i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_to_Tv[local_j] = idx3_Tv[path_Tv]

    def to_Tv_coords(vecs):
        result = np.zeros((vecs.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
        for local_j in range(len(old3)):
            if local_j in old3_to_Tv:
                j_Tv = old3_to_Tv[local_j]
                result[:, j_Tv] = (result[:, j_Tv] + vecs[:, local_j]) % PRIME
        return result

    psi_ker_Tv = to_Tv_coords(psi_ker) if psi_ker.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    oo_Tv = to_Tv_coords(oo_contrib) if oo_contrib.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    to_Tv = to_Tv_coords(to_contrib) if to_contrib.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)

    rk_psi_ker = int(_gauss_rank_np(psi_ker_Tv.copy(), PRIME)) if psi_ker_Tv.shape[0] > 0 else 0
    rk_oo = int(_gauss_rank_np(oo_Tv.copy(), PRIME)) if oo_Tv.shape[0] > 0 else 0
    rk_to = int(_gauss_rank_np(to_Tv.copy(), PRIME)) if to_Tv.shape[0] > 0 else 0

    # im(d_4^{T\v}) in A_3(T\v) coords
    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    im_d4_Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    rk_im_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

    # Test A: does old→old alone surject onto im(d_4^{T\v})?
    if oo_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combo_A = np.vstack([im_d4_Tv, oo_Tv]) % PRIME
        rk_A = int(_gauss_rank_np(combo_A.copy(), PRIME))
        oo_surjects = (rk_A == rk_im_d4_Tv)
    elif oo_Tv.shape[0] > 0:
        oo_surjects = (rk_oo == 0 and rk_im_d4_Tv == 0)
    else:
        oo_surjects = (rk_im_d4_Tv == 0)

    # Test B: does tv→old alone surject onto im(d_4^{T\v})?
    if to_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combo_B = np.vstack([im_d4_Tv, to_Tv]) % PRIME
        rk_B = int(_gauss_rank_np(combo_B.copy(), PRIME))
        to_surjects = (rk_B == rk_im_d4_Tv)
    elif to_Tv.shape[0] > 0:
        to_surjects = (rk_to == 0 and rk_im_d4_Tv == 0)
    else:
        to_surjects = (rk_im_d4_Tv == 0)

    # How much does old→old contribute vs tv→old?
    # Proportion: rk_oo / rk_psi_ker
    return {
        'rk_psi_ker': rk_psi_ker,
        'rk_oo': rk_oo,
        'rk_to': rk_to,
        'rk_im_d4_Tv': rk_im_d4_Tv,
        'oo_surjects': oo_surjects,
        'to_surjects': to_surjects,
        'dim_ker_tvtv': ker_tvtv.shape[0],
        'dim_ob4': ob4_T.shape[0],
        'n_old4_in_ker': sum(1 for i in range(ker_A4.shape[0]) if np.all(ker_A4[i, tv4] % PRIME == 0)) if ker_A4.shape[0] > 0 else 0,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"SURJECTION MECHANISM AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 300

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, min(n - 1, 6))
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                n1 = n - 1
                rem = [i for i in range(n) if i != v_cand]
                A_sub = [[A[rem[i]][rem[j]] for j in range(n1)] for i in range(n1)]
                cc_Tv = full_chain_complex_modp(A_sub, n1, min(n1-1, 5))
                if cc_Tv['bettis'].get(3, 0) != 1:
                    continue

                r = test_surjection(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
                if len(results) >= target:
                    break

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        oo_surj = sum(1 for r in results if r['oo_surjects'])
        to_surj = sum(1 for r in results if r['to_surjects'])
        print(f"\n  old→old alone surjects onto im(d_4^Tv): {oo_surj}/{len(results)}")
        print(f"  tv→old alone surjects onto im(d_4^Tv): {to_surj}/{len(results)}")

        for key in ['rk_psi_ker', 'rk_oo', 'rk_to', 'rk_im_d4_Tv', 'dim_ker_tvtv', 'n_old4_in_ker']:
            vals = [r[key] for r in results]
            if vals:
                print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

        # How often rk_oo + rk_to > rk_psi_ker (indicating overlap)?
        overlap = sum(1 for r in results if r['rk_oo'] + r['rk_to'] > r['rk_psi_ker'])
        print(f"\n  Overlap (rk_oo + rk_to > rk_psi): {overlap}/{len(results)}")

        # Distribution of n_old4_in_ker (pure old-only kernel vectors)
        old_only_dist = Counter(r['n_old4_in_ker'] for r in results)
        print(f"  Pure old-only ker vectors: {dict(sorted(old_only_dist.items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
