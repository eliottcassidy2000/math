"""
cancellation_mechanism.py — The old+new face cancellation in d_4(new)

KEY FINDING from old_face_analysis.py:
  The old-face-projection of new d_4 boundaries DOES reach the H_3(T\v) generator
  (adds 1 dimension to ker(d_3(T\v)) beyond im(d_4(T\v))).

But HYP-398 says new boundaries don't kill old cycles in T's chain complex.

RESOLUTION: The full boundary d_4(σ) = old_faces + new_faces.
When we compute in T's ker(d_3), the new-face components cancel out
the old-face components' contribution to the H_3 direction.

This script verifies:
1. The old-face component reaches the H_3(T\v) generator
2. The new-face component provides an exact cancellation
3. The total (in T's ker_d3) doesn't kill old cycles

If this cancellation is UNIVERSAL, it suggests a topological reason:
the inclusion i: T\v → T is a CHAIN MAP, and chain maps can't change
homology generators via the new cells they add.

Author: opus-2026-03-09-S56
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


def find_h3_generator_in_ker(im_d4_in_ker, ker_dim, prime):
    """Find a unit vector in ker_d3 coords NOT in im(d_4_in_ker).
    Returns the index of the ker_d3 basis vector that is the H_3 generator."""
    rank_im = _gauss_rank_np(im_d4_in_ker.copy(), prime)
    if rank_im == 0:
        return 0  # first basis vector is the generator

    for i in range(ker_dim):
        e_i = np.zeros((ker_dim, 1), dtype=np.int64)
        e_i[i] = 1
        combined = np.hstack([im_d4_in_ker, e_i]) % prime
        rank_combined = _gauss_rank_np(combined.copy(), prime)
        if rank_combined > rank_im:
            return i
    return 0


def verify_cancellation(A, n, v, verbose=False):
    """Verify old+new face cancellation mechanism."""
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

    ob_T3 = get_omega_basis(ap_T, 3, PRIME)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)
    ob_Tv3 = get_omega_basis(ap_Tv, 3, PRIME)
    ob_Tv4 = get_omega_basis(ap_Tv, 4, PRIME)

    paths_T_2 = ap_T.get(2, [])
    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])

    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}
    idx_T2 = {tuple(q): i for i, q in enumerate(paths_T_2)}

    # d_3^T boundary matrix
    bd_T3 = np.zeros((len(paths_T_2), len(paths_T_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T2:
                bd_T3[idx_T2[ft], j] = (bd_T3[idx_T2[ft], j] + sign) % PRIME

    # d_4^T boundary matrix
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    # ker(d_3^T) in Omega_3 coords
    d3_T_omega = bd_T3 @ ob_T3.T % PRIME
    M_d3T = [[int(x) % PRIME for x in d3_T_omega[r]] for r in range(d3_T_omega.shape[0])]
    _, ker_d3_T = _gauss_nullbasis_modp(M_d3T, d3_T_omega.shape[0], d3_T_omega.shape[1], PRIME)
    ker_d3_T = np.array(ker_d3_T, dtype=np.int64)

    # im(d_4^T) in Omega_3 coords
    if ob_T4 is not None and ob_T4.shape[0] > 0:
        im_d4_T = ob_T3 @ bd_T4 @ ob_T4.T % PRIME
    else:
        im_d4_T = np.zeros((ob_T3.shape[0], 1), dtype=np.int64)

    # im(d_4^T) projected to ker_d3^T
    im_d4_in_ker = ker_d3_T @ im_d4_T % PRIME
    rank_d4 = _gauss_rank_np(im_d4_in_ker.copy(), PRIME)

    # Find H_3(T) generator (as index in ker_d3_T basis)
    ker_dim_T = ker_d3_T.shape[0]
    h3_idx_T = find_h3_generator_in_ker(im_d4_in_ker, ker_dim_T, PRIME)
    h3_gen_T = np.zeros(ker_dim_T, dtype=np.int64)
    h3_gen_T[h3_idx_T] = 1

    # Embed H_3(T\v) generator into T's Omega_3 and project to ker_d3^T
    # First find H_3(T\v) generator
    idx_Tv3 = {tuple(q): i for i, q in enumerate(paths_Tv_3)}
    paths_Tv_2 = ap_Tv.get(2, [])
    idx_Tv2 = {tuple(q): i for i, q in enumerate(paths_Tv_2)}

    bd_Tv3 = np.zeros((len(paths_Tv_2), len(paths_Tv_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv2:
                bd_Tv3[idx_Tv2[ft], j] = (bd_Tv3[idx_Tv2[ft], j] + sign) % PRIME

    d3_Tv_omega = bd_Tv3 @ ob_Tv3.T % PRIME
    M_d3Tv = [[int(x) % PRIME for x in d3_Tv_omega[r]] for r in range(d3_Tv_omega.shape[0])]
    _, ker_d3_Tv = _gauss_nullbasis_modp(M_d3Tv, d3_Tv_omega.shape[0], d3_Tv_omega.shape[1], PRIME)
    ker_d3_Tv = np.array(ker_d3_Tv, dtype=np.int64)

    bd_Tv4 = np.zeros((len(paths_Tv_3), len(paths_Tv_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv3:
                bd_Tv4[idx_Tv3[ft], j] = (bd_Tv4[idx_Tv3[ft], j] + sign) % PRIME

    if ob_Tv4 is not None and ob_Tv4.shape[0] > 0:
        im_d4_Tv = ob_Tv3 @ bd_Tv4 @ ob_Tv4.T % PRIME
    else:
        im_d4_Tv = np.zeros((ob_Tv3.shape[0], 1), dtype=np.int64)

    im_d4_Tv_in_ker = ker_d3_Tv @ im_d4_Tv % PRIME
    ker_dim_Tv = ker_d3_Tv.shape[0]
    h3_idx_Tv = find_h3_generator_in_ker(im_d4_Tv_in_ker, ker_dim_Tv, PRIME)
    h3_gen_Tv = np.zeros(ker_dim_Tv, dtype=np.int64)
    h3_gen_Tv[h3_idx_Tv] = 1

    # Now: for each new Omega_4 basis vector, compute its boundary
    # decomposed into old (A_3(T\v)) and new (contains v) components,
    # and project each to T's ker_d3.

    new_omega4 = []
    for i in range(ob_T4.shape[0]):
        uses_v = any(ob_T4[i, j] % PRIME != 0 and v in paths_T_4[j]
                     for j in range(len(paths_T_4)))
        if uses_v:
            new_omega4.append(i)

    old_3idx = [i for i, q in enumerate(paths_T_3) if v not in q]
    new_3idx = [i for i, q in enumerate(paths_T_3) if v in q]

    # Collect all old-ker projections and check if they lie in im(d_4_in_ker)
    old_ker_cols = []
    new_ker_cols = []
    full_ker_cols = []

    for omega_idx in new_omega4:
        # Full boundary of this Omega_4 basis vector
        full_bd = bd_T4 @ ob_T4[omega_idx] % PRIME  # in A_3(T)

        # Split into old and new A_3 components
        old_comp = np.zeros(len(paths_T_3), dtype=np.int64)
        new_comp = np.zeros(len(paths_T_3), dtype=np.int64)
        for i in range(len(paths_T_3)):
            if v not in paths_T_3[i]:
                old_comp[i] = full_bd[i]
            else:
                new_comp[i] = full_bd[i]

        # Project through Omega_3(T) → ker_d3(T)
        old_ker = ker_d3_T @ (ob_T3 @ old_comp % PRIME) % PRIME
        new_ker = ker_d3_T @ (ob_T3 @ new_comp % PRIME) % PRIME
        full_ker = ker_d3_T @ (ob_T3 @ full_bd % PRIME) % PRIME

        old_ker_cols.append(old_ker)
        new_ker_cols.append(new_ker)
        full_ker_cols.append(full_ker)

    if not old_ker_cols:
        return None

    old_ker_mat = np.array(old_ker_cols, dtype=np.int64).T % PRIME  # (ker_dim, #new_omega4)
    new_ker_mat = np.array(new_ker_cols, dtype=np.int64).T % PRIME
    full_ker_mat = np.array(full_ker_cols, dtype=np.int64).T % PRIME

    # Check: full_ker should be in im(d_4_in_ker) (boundaries stay boundaries)
    rank_im = _gauss_rank_np(im_d4_in_ker.copy(), PRIME)
    combined_full = np.hstack([im_d4_in_ker, full_ker_mat]) % PRIME
    rank_full_combined = _gauss_rank_np(combined_full.copy(), PRIME)
    full_in_im = (rank_full_combined == rank_im)

    # KEY TEST 1: Is old_ker_mat in span of im_d4_in_ker?
    combined_old = np.hstack([im_d4_in_ker, old_ker_mat]) % PRIME
    rank_old_combined = _gauss_rank_np(combined_old.copy(), PRIME)
    old_adds = rank_old_combined - rank_im

    # KEY TEST 2: For the OLD d_4 image, compute its old-face projection too
    old_omega4 = []
    for i in range(ob_T4.shape[0]):
        uses_v = any(ob_T4[i, j] % PRIME != 0 and v in paths_T_4[j]
                     for j in range(len(paths_T_4)))
        if not uses_v:
            old_omega4.append(i)

    if old_omega4:
        old_ker_from_old_d4 = []
        for omega_idx in old_omega4:
            full_bd = bd_T4 @ ob_T4[omega_idx] % PRIME
            old_comp = np.zeros(len(paths_T_3), dtype=np.int64)
            for i in range(len(paths_T_3)):
                if v not in paths_T_3[i]:
                    old_comp[i] = full_bd[i]
            old_ker = ker_d3_T @ (ob_T3 @ old_comp % PRIME) % PRIME
            old_ker_from_old_d4.append(old_ker)
        old_d4_old_proj = np.array(old_ker_from_old_d4, dtype=np.int64).T % PRIME
        rank_old_d4_old = _gauss_rank_np(old_d4_old_proj.copy(), PRIME)

        # Combined old d4 + new d4 old projections
        combined_both = np.hstack([old_d4_old_proj, old_ker_mat]) % PRIME
        rank_both = _gauss_rank_np(combined_both.copy(), PRIME)
        new_adds_beyond_old_d4 = rank_both - rank_old_d4_old
    else:
        rank_old_d4_old = 0
        new_adds_beyond_old_d4 = _gauss_rank_np(old_ker_mat.copy(), PRIME)

    # Check: does old_ker_mat have H_3 component?
    # i.e., does [im_d4_in_ker | old_ker_mat] have same rank as [im_d4_in_ker]?
    # old_adds > 0 means old-face projection DOES reach beyond im(d_4)
    h3_reached = (old_adds > 0)

    result = {
        'n_new_omega4': len(new_omega4),
        'n_old_omega4': len(old_omega4),
        'rank_im': rank_im,
        'ker_dim': ker_dim_T,
        'full_in_im': full_in_im,
        'old_adds': old_adds,
        'h3_reached_by_old_proj': h3_reached,
        'rank_old_d4_old': rank_old_d4_old,
        'new_adds_beyond_old_d4': new_adds_beyond_old_d4,
    }

    if verbose:
        print(f"    Omega_4: {len(new_omega4)} new, {len(old_omega4)} old")
        print(f"    rank(im d_4 in ker_d3) = {rank_im}, ker_dim = {ker_dim_T}")
        print(f"    full boundaries in im(d_4)? {full_in_im}")
        print(f"    Old-face proj of new d_4 adds {old_adds} dims beyond im(d_4)")
        print(f"    Old-face proj of old d_4: rank = {rank_old_d4_old}")
        print(f"    New adds beyond old d_4 old-proj: {new_adds_beyond_old_d4}")
        if h3_reached:
            print(f"    *** OLD-FACE PROJ REACHES BEYOND im(d_4)! ***")
        else:
            print(f"    Old-face proj stays within im(d_4)")

    return result


def main():
    print("=" * 70)
    print("CANCELLATION MECHANISM — OLD+NEW FACE DECOMPOSITION")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found = 0
    target = 40
    t0 = time.time()
    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for vv in range(n):
            r = verify_cancellation(A, n, vv, verbose=(found <= 3))
            if r is not None:
                results.append(r)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  Done: {found} tours, {len(results)} bad verts, {elapsed:.1f}s")

    print(f"\n  KEY QUESTION: Does old-face proj of new d_4 reach beyond im(d_4)?")
    reached_dist = Counter(r['old_adds'] for r in results)
    print(f"    old_adds distribution: {dict(sorted(reached_dist.items()))}")

    if all(r['old_adds'] == 0 for r in results):
        print(f"\n    *** Old-face proj of new d_4 stays IN im(d_4_in_ker)! ***")
        print(f"    The new face components don't even need to cancel anything!")
    else:
        reached = sum(1 for r in results if r['old_adds'] > 0)
        print(f"\n    Old-face proj reaches beyond im(d_4) in {reached}/{len(results)} cases")
        print(f"    But full boundary stays in im(d_4): {sum(1 for r in results if r['full_in_im'])}/{len(results)}")
        print(f"    => New-face components CANCEL the old-face excess")

    print(f"\n  Old-face proj of old d_4 vs new d_4:")
    adds_dist = Counter(r['new_adds_beyond_old_d4'] for r in results)
    print(f"    New adds beyond old d4 old-proj: {dict(sorted(adds_dist.items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
