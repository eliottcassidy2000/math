"""
d4_preimage_omega.py — Corrected analysis in OMEGA coordinates

The previous script worked in A-path coords, which is wrong.
The correct chain complex is Omega_*, not A_*.

For BAD vertices (b3(T)=b3(T\\v)=1), we know:
- z = H_3(T\\v) generator, a nonzero class in H_3(T) (rank(i_*)=1)
- z ∉ im(d_4) in OMEGA coordinates

This script computes:
1. z_omega: the embedded H_3(T\\v) generator in Omega_3(T) coordinates
2. im(d_4)_omega: im(d_4|_{Omega_4(T)}) projected to Omega_3(T)
3. Split im(d_4)_omega into "old" (from Omega_4(T\\v)) and "new" contributions
4. Whether the new contribution contains z_omega

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
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    full_chain_complex_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME


def analyze_d4_omega(A, n, v, verbose=False):
    """Analyze d_4 preimage in Omega coordinates."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    # Omega bases
    def get_omega_basis(ap, deg, prime):
        paths = ap.get(deg, [])
        if not paths:
            return None
        P, nr, nc = _build_constraint_matrix(ap, deg, prime)
        if P is not None:
            r, nb = _gauss_nullbasis_modp(P, nr, nc, prime)
            return np.array(nb, dtype=np.int64) if nb else None
        else:
            return np.eye(len(paths), dtype=np.int64)

    ob_T3 = get_omega_basis(ap_T, 3, PRIME)  # (dim_omega3_T x |A_3(T)|)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)
    ob_Tv3 = get_omega_basis(ap_Tv, 3, PRIME)
    ob_Tv4 = get_omega_basis(ap_Tv, 4, PRIME)

    if ob_T3 is None or ob_Tv3 is None:
        return None

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_2 = ap_Tv.get(2, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])

    # Step 1: Compute H_3(T\v) generator z in Omega_3(T\v) coords
    idx_Tv2 = {tuple(q): i for i, q in enumerate(paths_Tv_2)}
    bd_Tv3 = np.zeros((len(paths_Tv_2), len(paths_Tv_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv2:
                bd_Tv3[idx_Tv2[ft], j] = (bd_Tv3[idx_Tv2[ft], j] + sign) % PRIME

    d3_Tv_omega = bd_Tv3 @ ob_Tv3.T % PRIME  # (|A_2(T\v)| x dim_omega3_Tv)
    M_list = [[int(x) % PRIME for x in d3_Tv_omega[r]] for r in range(d3_Tv_omega.shape[0])]
    _, ker_d3Tv = _gauss_nullbasis_modp(M_list, d3_Tv_omega.shape[0], d3_Tv_omega.shape[1], PRIME)
    ker_d3Tv = np.array(ker_d3Tv, dtype=np.int64)

    # z in A_3(T\v) coords
    z_Tv_Apath = (ob_Tv3.T @ ker_d3Tv[0].reshape(-1, 1)).flatten() % PRIME

    # Step 2: Embed z into A_3(T) coords
    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}
    z_T_Apath = np.zeros(len(paths_T_3), dtype=np.int64)
    for j, path_tv in enumerate(paths_Tv_3):
        if z_Tv_Apath[j] % PRIME != 0:
            embedded = embed_path(path_tv)
            if embedded in idx_T3:
                z_T_Apath[idx_T3[embedded]] = z_Tv_Apath[j] % PRIME

    # Step 3: Express z in Omega_3(T) coords
    # z_T_Apath is in A_3(T) coords. We need to find coefficients c such that
    # ob_T3.T @ c = z_T_Apath (mod PRIME).
    # This is a least-squares / consistency check.
    # Use Gauss: [ob_T3.T | z] and check rank increase
    z_col = z_T_Apath.reshape(-1, 1) % PRIME
    rank_ob = _gauss_rank_np(ob_T3.T.copy() % PRIME, PRIME)
    combined_z = np.hstack([ob_T3.T, z_col]) % PRIME
    rank_z_combined = _gauss_rank_np(combined_z, PRIME)
    z_in_omega = (rank_z_combined == rank_ob)

    if not z_in_omega and verbose:
        print(f"    WARNING: z not in Omega_3(T)! rank increase={rank_z_combined - rank_ob}")

    # Step 4: im(d_4|_{Omega_4(T)}) in A_3(T) coordinates
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    if ob_T4 is not None and ob_T4.shape[0] > 0:
        im_d4_Apath = bd_T4 @ ob_T4.T % PRIME  # (|A_3| x dim_omega4_T)
    else:
        im_d4_Apath = np.zeros((len(paths_T_3), 1), dtype=np.int64)

    # Step 5: Check z vs im(d_4) in A_3(T) coords
    rank_im = _gauss_rank_np(im_d4_Apath.copy() % PRIME, PRIME)
    combined_im = np.hstack([im_d4_Apath, z_col]) % PRIME
    rank_im_combined = _gauss_rank_np(combined_im, PRIME)
    z_in_im_d4 = (rank_im_combined == rank_im)

    # Step 6: Split im(d_4) into "old" (from Omega_4(T\v)) and "new"
    if ob_Tv4 is not None and ob_Tv4.shape[0] > 0:
        idx_T4 = {tuple(q): i for i, q in enumerate(paths_T_4)}
        embed_4 = np.zeros((len(paths_T_4), len(paths_Tv_4)), dtype=np.int64)
        for j, ptv in enumerate(paths_Tv_4):
            emb = embed_path(ptv)
            if emb in idx_T4:
                embed_4[idx_T4[emb], j] = 1
        embedded_omega_Tv4 = embed_4 @ ob_Tv4.T % PRIME
        im_old = bd_T4 @ embedded_omega_Tv4 % PRIME

        rank_old = _gauss_rank_np(im_old.copy() % PRIME, PRIME)
        combined_old = np.hstack([im_old, z_col]) % PRIME
        rank_old_combined = _gauss_rank_np(combined_old, PRIME)
        z_in_old = (rank_old_combined == rank_old)
    else:
        rank_old = 0
        z_in_old = False

    result = {
        'z_in_omega': z_in_omega,
        'z_in_im_d4': z_in_im_d4,
        'z_in_im_old': z_in_old,
        'rank_im_d4': rank_im,
        'rank_im_old': rank_old,
        'rank_new_incr': rank_im - rank_old,
        'dim_omega3': ob_T3.shape[0],
        'dim_omega4': ob_T4.shape[0] if ob_T4 is not None else 0,
        'ker_d3': cc_T['kers'].get(3, 0),
        'beta3_T': cc_T['bettis'].get(3, 0),
    }

    if verbose:
        print(f"    z in Omega_3(T): {z_in_omega}")
        print(f"    z in im(d_4|_Omega): {z_in_im_d4}")
        print(f"    z in im(d_4|_old): {z_in_old}")
        print(f"    rank(d_4|_Omega) = {rank_im}, rank(d_4|_old) = {rank_old}, "
              f"new increment = {rank_im - rank_old}")
        print(f"    ker(d_3) = {cc_T['kers'].get(3, 0)}, "
              f"Omega_3 = {ob_T3.shape[0]}, Omega_4 = {ob_T4.shape[0] if ob_T4 is not None else 0}")

    return result


def main():
    print("=" * 70)
    print("d_4 PREIMAGE IN OMEGA COORDINATES")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found = 0
    target = 60
    t0 = time.time()
    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            r = analyze_d4_omega(A, n, v, verbose=(found <= 2))
            if r is not None:
                results.append(r)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tours, {len(results)} bad verts, {elapsed:.1f}s")

    # Summary
    print(f"\n  KEY RESULTS:")
    z_in_omega = sum(1 for r in results if r['z_in_omega'])
    z_in_im = sum(1 for r in results if r['z_in_im_d4'])
    z_in_old = sum(1 for r in results if r['z_in_im_old'])
    print(f"    z in Omega_3(T): {z_in_omega}/{len(results)}")
    print(f"    z in im(d_4|_Omega): {z_in_im}/{len(results)}")
    print(f"    z in im(d_4|_old): {z_in_old}/{len(results)}")

    print(f"\n  Rank distributions:")
    rank_im_dist = Counter(r['rank_im_d4'] for r in results)
    rank_old_dist = Counter(r['rank_im_old'] for r in results)
    incr_dist = Counter(r['rank_new_incr'] for r in results)
    print(f"    rank(d_4|_Omega): {dict(sorted(rank_im_dist.items()))}")
    print(f"    rank(d_4|_old): {dict(sorted(rank_old_dist.items()))}")
    print(f"    new increment: {dict(sorted(incr_dist.items()))}")

    # Check: does rank(d_4) = ker(d_3) - 1 for beta_3=1?
    ker_minus_rank = Counter(r['ker_d3'] - r['rank_im_d4'] for r in results)
    print(f"\n  ker(d_3) - rank(d_4) = {dict(ker_minus_rank)}")
    print(f"  (Should be 1 for beta_3=1)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
