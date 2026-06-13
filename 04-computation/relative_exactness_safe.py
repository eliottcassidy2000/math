"""
relative_exactness_safe.py — H_4^rel with safe modular arithmetic

Re-verifies that H_4^rel = 0 universally at n=7 AND n=8, using
matmul_mod to avoid int64 overflow (MISTAKE-019).

Key question: does rk(d_5^rel) = dim(ker d_4^rel) at BOTH n=7 and n=8?
(Kind-pasteur-S50 says YES with safe arithmetic.)

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
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


def safe_relative_exactness(A, n, v):
    """Compute H_4^rel using safe arithmetic throughout."""
    max_p = min(n - 1, 6)
    ap_T = enumerate_all_allowed(A, n, max_p)

    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1-1))

    cc_Tv = full_chain_complex_modp(A_sub, n1, min(n1-1, 5))
    if cc_Tv['bettis'].get(3, 0) != 1:
        return None

    paths = {p: ap_T.get(p, []) for p in range(2, max_p + 1)}
    paths_Tv = {p: ap_Tv.get(p, []) for p in range(2, min(max_p, n1-1) + 1)}

    def classify(ps, v):
        tv = [i for i, p in enumerate(ps) if v in p]
        old = [i for i, p in enumerate(ps) if v not in p]
        return tv, old

    tv = {}; old = {}
    for p in range(2, max_p + 1):
        tv[p], old[p] = classify(paths.get(p, []), v)

    if not tv.get(4):
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

    ob_T = {}
    for p in range(2, max_p + 1):
        ob_T[p] = get_omega(ap_T, p)

    ob_Tv = {}
    for p in range(2, min(max_p, n1-1) + 1):
        ob_Tv[p] = get_omega(ap_Tv, p)

    # Build boundary maps
    idx = {}
    for p in range(2, max_p + 1):
        idx[p] = {path: i for i, path in enumerate(paths[p])}

    bd = {}
    for p in range(3, max_p + 1):
        src = paths.get(p, [])
        tgt = paths.get(p-1, [])
        if not src or not tgt:
            continue
        mat = np.zeros((len(tgt), len(src)), dtype=np.int64)
        for j, path in enumerate(src):
            for sign, face in boundary_faces(path):
                if face in idx[p-1]:
                    mat[idx[p-1][face], j] = (mat[idx[p-1][face], j] + sign) % PRIME
        bd[p] = mat

    # d_4^{tv→tv} through Omega: SAFE
    if 4 in bd and ob_T[4].shape[0] > 0:
        d4_tv_omega = matmul_mod(bd[4][tv[3], :], ob_T[4].T)
    else:
        d4_tv_omega = np.zeros((len(tv.get(3, [])), 0), dtype=np.int64)

    # ker(d_4^{tv→tv})
    if d4_tv_omega.size > 0 and d4_tv_omega.shape[1] > 0:
        rk_d4_tv = int(_gauss_rank_np(d4_tv_omega.copy(), PRIME))
        _, kv = _gauss_nullbasis_modp(
            d4_tv_omega.astype(np.int32), d4_tv_omega.shape[0], d4_tv_omega.shape[1], PRIME)
        dim_ker_d4_rel = len(kv) if kv else 0
        ker_tvtv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob_T[4].shape[0]), dtype=np.int64)
    else:
        rk_d4_tv = 0
        dim_ker_d4_rel = ob_T[4].shape[0]
        ker_tvtv = np.eye(ob_T[4].shape[0], dtype=np.int64) if ob_T[4].shape[0] > 0 else np.zeros((0, 0), dtype=np.int64)

    # ψ(ker): SAFE — full boundary of relative cycles projected to old3
    if ker_tvtv.shape[0] > 0:
        ker_A4 = matmul_mod(ker_tvtv, ob_T[4])  # SAFE
        psi_full = matmul_mod(bd[4][old[3], :], ker_A4.T).T  # SAFE
    else:
        ker_A4 = np.zeros((0, len(paths.get(4, []))), dtype=np.int64)
        psi_full = np.zeros((0, len(old[3])), dtype=np.int64)

    # Map to T\v coords
    idx3_Tv = {path: i for i, path in enumerate(paths_Tv.get(3, []))}
    remaining_inv = {remaining[i]: i for i in range(n1)}
    old3_to_Tv = {}
    for local_j, global_i in enumerate(old[3]):
        path_T = paths[3][global_i]
        path_Tv = tuple(remaining_inv[x] for x in path_T)
        if path_Tv in idx3_Tv:
            old3_to_Tv[local_j] = idx3_Tv[path_Tv]

    psi_Tv = np.zeros((psi_full.shape[0], len(paths_Tv.get(3, []))), dtype=np.int64)
    for local_j in range(len(old[3])):
        if local_j in old3_to_Tv:
            j_Tv = old3_to_Tv[local_j]
            psi_Tv[:, j_Tv] = (psi_Tv[:, j_Tv] + psi_full[:, local_j]) % PRIME

    rk_psi = int(_gauss_rank_np(psi_Tv.copy(), PRIME)) if psi_Tv.shape[0] > 0 else 0

    # im(d_4^{T\v}) in A_3(T\v) coords: SAFE
    bd4_Tv = np.zeros((len(paths_Tv.get(3, [])), len(paths_Tv.get(4, []))), dtype=np.int64)
    for j, path in enumerate(paths_Tv.get(4, [])):
        for sign, face in boundary_faces(path):
            if face in idx3_Tv:
                bd4_Tv[idx3_Tv[face], j] = (bd4_Tv[idx3_Tv[face], j] + sign) % PRIME

    ob4_Tv = ob_Tv.get(4, np.zeros((0, 0), dtype=np.int64))
    if ob4_Tv.shape[0] > 0:
        im_d4_Tv = matmul_mod(bd4_Tv, ob4_Tv.T).T  # SAFE
    else:
        im_d4_Tv = np.zeros((0, len(paths_Tv.get(3, []))), dtype=np.int64)
    rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

    # Check containment
    if psi_Tv.shape[0] > 0 and im_d4_Tv.shape[0] > 0:
        combined = np.vstack([im_d4_Tv, psi_Tv]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        psi_in_im = (rk_combined == rk_d4_Tv)
        extra = rk_combined - rk_d4_Tv
    elif psi_Tv.shape[0] > 0:
        psi_in_im = (rk_psi == 0)
        extra = rk_psi
    else:
        psi_in_im = True
        extra = 0

    # Compute d_5^rel rank (for relative exactness): SAFE
    if 5 in bd and ob_T.get(5, np.zeros((0,0))).shape[0] > 0:
        bd5_omega = matmul_mod(bd[5], ob_T[5].T)  # SAFE
    else:
        bd5_omega = np.zeros((len(paths.get(4, [])), 0), dtype=np.int64)

    im_d5 = bd5_omega.T  # in A_4(T) coords

    # Embed Omega_4(T\v) into A_4(T) coords
    embed_4 = np.zeros((len(paths_Tv.get(4, [])), len(paths.get(4, []))), dtype=np.int64)
    idx4_T = {path: i for i, path in enumerate(paths.get(4, []))}
    for jv, pv in enumerate(paths_Tv.get(4, [])):
        mapped = tuple(remaining[x] for x in pv)
        if mapped in idx4_T:
            embed_4[jv, idx4_T[mapped]] = 1

    if ob4_Tv.shape[0] > 0:
        i_Omega4_Tv = matmul_mod(ob4_Tv, embed_4)  # SAFE
    else:
        i_Omega4_Tv = np.zeros((0, len(paths.get(4, []))), dtype=np.int64)

    dim_Omega4_Tv = i_Omega4_Tv.shape[0]

    if im_d5.shape[0] > 0 and i_Omega4_Tv.shape[0] > 0:
        combined_5 = np.vstack([im_d5, i_Omega4_Tv]) % PRIME
        rk_combined_5 = int(_gauss_rank_np(combined_5.copy(), PRIME))
        rk_d5_rel = rk_combined_5 - dim_Omega4_Tv
    elif im_d5.shape[0] > 0:
        rk_d5_rel = int(_gauss_rank_np(im_d5.copy(), PRIME))
    else:
        rk_d5_rel = 0

    dim_C4_rel = ob_T[4].shape[0] - dim_Omega4_Tv
    dim_ker_d4_rel_check = dim_C4_rel - rk_d4_tv  # cross-check
    h4_rel = dim_ker_d4_rel - rk_d5_rel

    return {
        'h4_rel': h4_rel,
        'dim_ker_d4_rel': dim_ker_d4_rel,
        'rk_d5_rel': rk_d5_rel,
        'rk_psi': rk_psi,
        'rk_d4_Tv': rk_d4_Tv,
        'psi_in_im': psi_in_im,
        'extra': extra,
        'dim_C4_rel': dim_C4_rel,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"SAFE RELATIVE EXACTNESS AT DEGREE 4, n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 200 if n <= 7 else 500

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, min(n - 1, 6))
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = safe_relative_exactness(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
                if len(results) >= target:
                    break

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        h4_dist = Counter(r['h4_rel'] for r in results)
        print(f"\n  H_4^rel dimension: {dict(sorted(h4_dist.items()))}")

        in_im = sum(1 for r in results if r['psi_in_im'])
        print(f"  ψ(ker) ⊂ im(d_4^Tv): {in_im}/{len(results)}")

        extra_dist = Counter(r['extra'] for r in results)
        print(f"  Extra directions in H_3(T\\v): {dict(sorted(extra_dist.items()))}")

        codims = [r['rk_d4_Tv'] - r['rk_psi'] for r in results if r['psi_in_im']]
        codim_dist = Counter(codims)
        print(f"  codim(ψ(ker), im d_4^Tv): {dict(sorted(codim_dist.items()))}")

        slacks = [r['dim_ker_d4_rel'] - r['rk_d5_rel'] for r in results]
        slack_dist = Counter(slacks)
        print(f"  Slack (ker d4^rel - im d5^rel): {dict(sorted(slack_dist.items()))}")

        for key in ['dim_ker_d4_rel', 'rk_d5_rel', 'rk_psi', 'rk_d4_Tv', 'dim_C4_rel']:
            vals = [r[key] for r in results]
            print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
