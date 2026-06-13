"""
d4_preimage_analysis.py — Why can't the H_3(T\\v) generator become a boundary in T?

For BAD vertices (b3(T)=b3(T\\v)=1):
- z = H_3(T\\v) generator, embedded in Omega_3(T)
- z is a cycle: d_3(z) = 0
- z is NOT a boundary: z ∉ im(d_4^T)

Question: What is the obstruction to z being in im(d_4^T)?

Idea: The "new" d_4 content (from 5-paths through v) maps to a subspace
of Omega_3(T) that AVOIDS z. Why?

This script analyzes:
1. The "new" d_4 image: d_4(Omega_4^new) where Omega_4^new uses vertex v
2. The "old" d_4 image: d_4(Omega_4^old) = im(d_4^{T\\v}) embedded
3. The relationship between z and im(d_4^T) = im(d_4^old) + im(d_4^new)

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


def analyze_d4_preimage(A, n, v, verbose=False):
    """Analyze the structure of d_4 preimage for a specific tournament and vertex.

    Returns dict with:
    - is_bad: whether b3(T\v) = 1
    - rank_old: rank of d_4 from "old" (non-v) 5-paths
    - rank_new: rank of d_4 from "new" (v-containing) 5-paths
    - rank_total: total rank of d_4
    - z_in_old: whether z is in im(d_4^old)
    - z_in_new: whether z is in im(d_4^new)
    - z_in_total: whether z is in im(d_4^T)
    """
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    if b3_T != 1 or b3_Tv != 1:
        return {'is_bad': False, 'skip': True}

    # Build A-path data
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])
    paths_Tv_2 = ap_Tv.get(2, [])

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    # Omega_3(T) basis
    P_T3, nr3, nc3 = _build_constraint_matrix(ap_T, 3, PRIME)
    if P_T3 is not None:
        r3, nb3 = _gauss_nullbasis_modp(P_T3, nr3, nc3, PRIME)
        omega_basis_T3 = np.array(nb3, dtype=np.int64)
    else:
        omega_basis_T3 = np.eye(len(paths_T_3), dtype=np.int64)

    # Omega_4(T) basis
    P_T4, nr4, nc4 = _build_constraint_matrix(ap_T, 4, PRIME)
    if P_T4 is not None:
        r4, nb4 = _gauss_nullbasis_modp(P_T4, nr4, nc4, PRIME)
        omega_basis_T4 = np.array(nb4, dtype=np.int64) if nb4 else None
    else:
        omega_basis_T4 = np.eye(len(paths_T_4), dtype=np.int64) if paths_T_4 else None

    # Omega_3(T\v) basis
    P_Tv3, nrv3, ncv3 = _build_constraint_matrix(ap_Tv, 3, PRIME)
    if P_Tv3 is not None:
        rv3, nbv3 = _gauss_nullbasis_modp(P_Tv3, nrv3, ncv3, PRIME)
        omega_basis_Tv3 = np.array(nbv3, dtype=np.int64)
    else:
        omega_basis_Tv3 = np.eye(len(paths_Tv_3), dtype=np.int64)

    # Compute H_3(T\v) generator: cycle in ker(d_3^{T\v})
    idx_Tv2 = {tuple(q): i for i, q in enumerate(paths_Tv_2)}
    bd_Tv3 = np.zeros((len(paths_Tv_2), len(paths_Tv_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv2:
                bd_Tv3[idx_Tv2[ft], j] = (bd_Tv3[idx_Tv2[ft], j] + sign) % PRIME

    d3_Tv_omega = bd_Tv3 @ omega_basis_Tv3.T % PRIME
    M_list = [[int(x) % PRIME for x in d3_Tv_omega[r]] for r in range(d3_Tv_omega.shape[0])]
    _, ker_d3Tv = _gauss_nullbasis_modp(M_list, d3_Tv_omega.shape[0], d3_Tv_omega.shape[1], PRIME)
    ker_d3Tv = np.array(ker_d3Tv, dtype=np.int64)  # rows = cycle basis in Omega_3(T\v) coords

    # H_3(T\v) generator: first cycle basis vector (representative)
    z_Tv_omega = ker_d3Tv[0]  # in Omega_3(T\v) coords
    z_Tv_Apath = (omega_basis_Tv3.T @ z_Tv_omega.reshape(-1, 1)).flatten() % PRIME  # in A_3(T\v) coords

    # Embed z in A_3(T) coordinates
    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}
    z_T_Apath = np.zeros(len(paths_T_3), dtype=np.int64)
    for j, path_tv in enumerate(paths_Tv_3):
        if z_Tv_Apath[j] % PRIME != 0:
            embedded = embed_path(path_tv)
            if embedded in idx_T3:
                z_T_Apath[idx_T3[embedded]] = z_Tv_Apath[j] % PRIME

    # Boundary matrix d_4: A_4(T) -> A_3(T)
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    # Classify 5-paths as old (not using v) or new (using v)
    old_indices = []
    new_indices = []
    for j, path in enumerate(paths_T_4):
        if v in path:
            new_indices.append(j)
        else:
            old_indices.append(j)

    # im(d_4) from old and new paths (in A_3(T) coords)
    if omega_basis_T4 is not None:
        im_d4_total_Apath = bd_T4 @ omega_basis_T4.T % PRIME

        # Split Omega_4 basis vectors into old/new components
        # An Omega_4 basis vector is a linear combination of A_4 paths.
        # "Old" component: coefficients on non-v paths only
        # "New" component: coefficients on v-paths only

        # Actually, we should split at the Omega level, not path level.
        # But that's hard. Instead: compute d_4 image from:
        # - Omega_4(T\v) embedded (the "old" image)
        # - The complement (the "new" contribution)

        # Omega_4(T\v) embedded in A_4(T)
        if paths_Tv_4:
            P_Tv4, nrv4, ncv4 = _build_constraint_matrix(ap_Tv, 4, PRIME)
            if P_Tv4 is not None:
                rv4, nbv4 = _gauss_nullbasis_modp(P_Tv4, nrv4, ncv4, PRIME)
                omega_basis_Tv4 = np.array(nbv4, dtype=np.int64) if nbv4 else None
            else:
                omega_basis_Tv4 = np.eye(len(paths_Tv_4), dtype=np.int64)

            if omega_basis_Tv4 is not None:
                # Embed Omega_4(T\v) into A_4(T)
                idx_T4 = {tuple(q): i for i, q in enumerate(paths_T_4)}
                embed_4 = np.zeros((len(paths_T_4), len(paths_Tv_4)), dtype=np.int64)
                for j, ptv in enumerate(paths_Tv_4):
                    emb = embed_path(ptv)
                    if emb in idx_T4:
                        embed_4[idx_T4[emb], j] = 1

                embedded_omega_Tv4 = embed_4 @ omega_basis_Tv4.T % PRIME
                im_old = bd_T4 @ embedded_omega_Tv4 % PRIME  # "old" d_4 image in A_3(T)
            else:
                im_old = np.zeros((len(paths_T_3), 1), dtype=np.int64)
        else:
            im_old = np.zeros((len(paths_T_3), 1), dtype=np.int64)

        # Check: is z in im(d_4^T)?
        z_col = z_T_Apath.reshape(-1, 1) % PRIME
        rank_im = _gauss_rank_np(im_d4_total_Apath.copy() % PRIME, PRIME)
        combined = np.hstack([im_d4_total_Apath, z_col]) % PRIME
        rank_combined = _gauss_rank_np(combined, PRIME)
        z_in_total = (rank_combined == rank_im)

        # Is z in im(d_4^old)?
        rank_old = _gauss_rank_np(im_old.copy() % PRIME, PRIME)
        combined_old = np.hstack([im_old, z_col]) % PRIME
        rank_combined_old = _gauss_rank_np(combined_old, PRIME)
        z_in_old = (rank_combined_old == rank_old)

        # im(d_4^new) = complement of old in total
        # Compute: rank of [im_old | im_new_supplement | z]
        # where im_new_supplement = total minus old
        rank_total = _gauss_rank_np(im_d4_total_Apath.copy() % PRIME, PRIME)

        result = {
            'is_bad': True,
            'skip': False,
            'dim_omega3_T': omega_basis_T3.shape[0],
            'dim_omega4_T': omega_basis_T4.shape[0] if omega_basis_T4 is not None else 0,
            'rank_d4_total': rank_total,
            'rank_d4_old': rank_old,
            'rank_d4_new_increment': rank_total - rank_old,
            'z_in_total': z_in_total,
            'z_in_old': z_in_old,
            'n_old_5paths': len(old_indices),
            'n_new_5paths': len(new_indices),
        }

        if verbose:
            print(f"    Omega dims: Omega_3={result['dim_omega3_T']}, Omega_4={result['dim_omega4_T']}")
            print(f"    rank(d_4^total)={rank_total}, rank(d_4^old)={rank_old}, "
                  f"new increment={rank_total - rank_old}")
            print(f"    z in im(d_4^total)? {z_in_total}")
            print(f"    z in im(d_4^old)? {z_in_old}")
            print(f"    A-paths: {len(old_indices)} old, {len(new_indices)} new 5-paths")

        return result
    else:
        return {'is_bad': True, 'skip': True, 'no_omega4': True}


def main():
    print("=" * 70)
    print("d_4 PREIMAGE ANALYSIS — WHY z ∉ im(d_4^T)?")
    print("=" * 70)

    n = 7
    max_p = n - 1
    rng = np.random.RandomState(42)

    print(f"\n--- n={n}: beta_3=1 tournaments, BAD vertices ---")

    found = 0
    target = 50
    t0 = time.time()

    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            result = analyze_d4_preimage(A, n, v, verbose=(found <= 2))
            if not result.get('skip', True):
                results.append(result)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments, {len(results)} bad vertex records, {elapsed:.1f}s")

    # Analyze
    print(f"\n  KEY QUESTION: Does z ever fall in im(d_4^old)?")
    z_in_old_count = sum(1 for r in results if r.get('z_in_old', False))
    z_in_total_count = sum(1 for r in results if r.get('z_in_total', False))
    print(f"    z in im(d_4^old): {z_in_old_count}/{len(results)}")
    print(f"    z in im(d_4^total): {z_in_total_count}/{len(results)}")

    print(f"\n  Rank analysis:")
    rank_old_dist = Counter(r.get('rank_d4_old', 0) for r in results)
    rank_total_dist = Counter(r.get('rank_d4_total', 0) for r in results)
    incr_dist = Counter(r.get('rank_d4_new_increment', 0) for r in results)
    print(f"    rank(d_4^old) distribution: {dict(rank_old_dist)}")
    print(f"    rank(d_4^total) distribution: {dict(rank_total_dist)}")
    print(f"    new increment distribution: {dict(incr_dist)}")

    # Key insight: if rank(d_4^old) = rank(d_4^total), then new paths add NOTHING
    # to the d_4 image. This would trivially imply z is not in im(d_4) since
    # z is not in im(d_4^old) (by the inductive hypothesis that b3(T\v)=1).
    no_new_contribution = sum(1 for r in results if r.get('rank_d4_new_increment', 0) == 0)
    print(f"\n  New paths add NOTHING to d_4 image: {no_new_contribution}/{len(results)}")

    # What about the size of the new contribution?
    print(f"\n  Omega dimensions:")
    omega3_dist = Counter(r.get('dim_omega3_T', 0) for r in results)
    omega4_dist = Counter(r.get('dim_omega4_T', 0) for r in results)
    print(f"    Omega_3(T): {dict(omega3_dist)}")
    print(f"    Omega_4(T): {dict(omega4_dist)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
