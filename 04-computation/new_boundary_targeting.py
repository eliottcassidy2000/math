"""
new_boundary_targeting.py — Do new d_4 boundaries ONLY kill new d_3 cycles?

The claim: im(d_4^new) ∩ span(old_ker_d3) = im(d_4^old)
i.e., the new boundaries don't add anything to the killing of old cycles.

This would explain WHY the old H_3 generator survives.

If true: im(d_4^T) = im(d_4^old) + im(d_4^new)
and im(d_4^new) is entirely in the "new cycle" subspace.

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


def check_targeting(A, n, v, verbose=False):
    """Check: does im(d_4^new) target only new ker_d3 cycles?"""
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
    paths_Tv_2 = ap_Tv.get(2, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])

    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}

    # ker(d_3^{T\v}) in A_3(T) coords (embedded)
    idx_Tv2 = {tuple(q): i for i, q in enumerate(paths_Tv_2)}
    bd_Tv3 = np.zeros((len(paths_Tv_2), len(paths_Tv_3)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv2:
                bd_Tv3[idx_Tv2[ft], j] = (bd_Tv3[idx_Tv2[ft], j] + sign) % PRIME

    d3_Tv_omega = bd_Tv3 @ ob_Tv3.T % PRIME
    M_dTv3 = [[int(x) % PRIME for x in d3_Tv_omega[r]] for r in range(d3_Tv_omega.shape[0])]
    _, ker_d3_Tv = _gauss_nullbasis_modp(M_dTv3, d3_Tv_omega.shape[0], d3_Tv_omega.shape[1], PRIME)
    ker_d3_Tv = np.array(ker_d3_Tv, dtype=np.int64)

    # Embed ker_d3(T\v) into A_3(T)
    embed_3 = np.zeros((len(paths_T_3), len(paths_Tv_3)), dtype=np.int64)
    for j, ptv in enumerate(paths_Tv_3):
        emb = embed_path(ptv)
        if emb in idx_T3:
            embed_3[idx_T3[emb], j] = 1
    embedded_Tv3 = embed_3 @ ob_Tv3.T % PRIME

    # Embedded old cycles in A_3(T): embedded_Tv3 @ ker_d3_Tv.T
    old_cycles_Apath = embedded_Tv3 @ ker_d3_Tv.T % PRIME  # (|A_3(T)| x ker_dim_Tv)

    # im(d_4^old) in A_3(T)
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    if ob_Tv4 is not None and ob_Tv4.shape[0] > 0:
        idx_T4 = {tuple(q): i for i, q in enumerate(paths_T_4)}
        embed_4 = np.zeros((len(paths_T_4), len(paths_Tv_4)), dtype=np.int64)
        for j, ptv in enumerate(paths_Tv_4):
            emb = embed_path(ptv)
            if emb in idx_T4:
                embed_4[idx_T4[emb], j] = 1
        embedded_Tv4 = embed_4 @ ob_Tv4.T % PRIME
        im_old = bd_T4 @ embedded_Tv4 % PRIME
    else:
        im_old = np.zeros((len(paths_T_3), 1), dtype=np.int64)

    # im(d_4^T) total
    if ob_T4 is not None and ob_T4.shape[0] > 0:
        im_total = bd_T4 @ ob_T4.T % PRIME
    else:
        im_total = np.zeros((len(paths_T_3), 1), dtype=np.int64)

    rank_old = _gauss_rank_np(im_old.copy() % PRIME, PRIME)
    rank_total = _gauss_rank_np(im_total.copy() % PRIME, PRIME)

    # KEY TEST: rank of [im_total | old_cycles] vs rank of [im_old | old_cycles]
    # If im(d_4^new) doesn't reach old cycles beyond im(d_4^old), then:
    # rank([im_total | old_cycles]) = rank_total + (ker_dim_Tv - rank_old)
    # = same as rank([im_old | old_cycles]) + (rank_total - rank_old)

    combined_old_cycles = np.hstack([im_old, old_cycles_Apath]) % PRIME
    rank_old_and_cycles = _gauss_rank_np(combined_old_cycles, PRIME)

    combined_total_cycles = np.hstack([im_total, old_cycles_Apath]) % PRIME
    rank_total_and_cycles = _gauss_rank_np(combined_total_cycles, PRIME)

    # How many old cycles does im_old kill?
    old_killed_by_old = old_cycles_Apath.shape[1] - (rank_old_and_cycles - rank_old)

    # How many old cycles does im_total kill?
    old_killed_by_total = old_cycles_Apath.shape[1] - (rank_total_and_cycles - rank_total)

    # Do new boundaries kill MORE old cycles than old boundaries?
    new_kills_extra_old = old_killed_by_total - old_killed_by_old

    result = {
        'ker_dim_Tv': ker_d3_Tv.shape[0],
        'rank_old': rank_old,
        'rank_total': rank_total,
        'old_killed_by_old': old_killed_by_old,
        'old_killed_by_total': old_killed_by_total,
        'new_kills_extra_old': new_kills_extra_old,
    }

    if verbose:
        print(f"    ker_d3(T\\v) = {ker_d3_Tv.shape[0]}")
        print(f"    rank(d_4^old) = {rank_old}, rank(d_4^total) = {rank_total}")
        print(f"    Old cycles killed by im(d_4^old): {old_killed_by_old}")
        print(f"    Old cycles killed by im(d_4^total): {old_killed_by_total}")
        print(f"    New boundaries kill {new_kills_extra_old} EXTRA old cycles")

    return result


def main():
    print("=" * 70)
    print("NEW BOUNDARY TARGETING ANALYSIS")
    print("Do new d_4 boundaries kill old ker_d3 cycles?")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found = 0
    target = 50
    t0 = time.time()
    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            r = check_targeting(A, n, v, verbose=(found <= 2))
            if r is not None:
                results.append(r)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tours, {len(results)} bad verts, {elapsed:.1f}s")

    print(f"\n  KEY RESULT: Do new boundaries kill extra old cycles?")
    extra_kills = Counter(r['new_kills_extra_old'] for r in results)
    print(f"    Extra old cycles killed by new boundaries: {dict(sorted(extra_kills.items()))}")

    if all(r['new_kills_extra_old'] == 0 for r in results):
        print(f"\n    *** CONFIRMED: New boundaries NEVER kill old cycles ***")
        print(f"    This explains i_*-injectivity:")
        print(f"    - Old H_3 generator is an old cycle")
        print(f"    - im(d_4^old) kills ker_dim_Tv - 1 old cycles (leaving the generator)")
        print(f"    - New boundaries add nothing to old cycle killing")
        print(f"    - Therefore the old generator survives in H_3(T)")
    else:
        zero_count = sum(1 for r in results if r['new_kills_extra_old'] == 0)
        print(f"\n    New boundaries DO sometimes kill old cycles: "
              f"{len(results) - zero_count}/{len(results)} cases")


if __name__ == '__main__':
    main()
    print("\nDONE.")
