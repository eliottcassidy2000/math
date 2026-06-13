"""
istar_failure_n8.py — Find and analyze i_*-injectivity failures at n=8

At n=7: rank(i_*^3) = 1 for ALL BAD vertices (HYP-380 confirmed).
At n=8: rank(i_*^3) = 0 for ~0.26% of BAD vertices (HYP-380 REFUTED).

When rank(i_*)=0, the H_3(T\v) generator IS killed in T.
Question: HOW is it killed? Via old boundaries or new boundaries?
Does HYP-398 (new boundaries target only new cycles) still hold in these cases?

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


def analyze_failure(A, n, v, verbose=False):
    """Detailed analysis when rank(i_*)=0 at n=8."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)
    b4_T = cc_T['bettis'].get(4, 0)
    b4_Tv = cc_Tv['bettis'].get(4, 0)

    if b3_T != 1 or b3_Tv != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    ob_T3 = get_omega_basis(ap_T, 3, PRIME)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)
    ob_Tv3 = get_omega_basis(ap_Tv, 3, PRIME)
    ob_Tv4 = get_omega_basis(ap_Tv, 4, PRIME)

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])

    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}

    # ker(d_3^{T\v}) and im(d_4^{T\v})
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

    idx_Tv3 = {tuple(q): i for i, q in enumerate(paths_Tv_3)}
    bd_Tv4 = np.zeros((len(paths_Tv_3), len(paths_Tv_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_Tv_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_Tv3:
                bd_Tv4[idx_Tv3[ft], j] = (bd_Tv4[idx_Tv3[ft], j] + sign) % PRIME

    if ob_Tv4 is not None and ob_Tv4.shape[0] > 0:
        im_d4_Tv = bd_Tv4 @ ob_Tv4.T % PRIME
    else:
        im_d4_Tv = np.zeros((len(paths_Tv_3), 1), dtype=np.int64)

    # Embed into T
    embed_3 = np.zeros((len(paths_T_3), len(paths_Tv_3)), dtype=np.int64)
    for j, ptv in enumerate(paths_Tv_3):
        emb = embed_path(ptv)
        if emb in idx_T3:
            embed_3[idx_T3[emb], j] = 1
    embedded_Tv3 = embed_3 @ ob_Tv3.T % PRIME
    old_cycles_Apath = embedded_Tv3 @ ker_d3_Tv.T % PRIME

    # d_4^T boundary
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    # im(d_4^old) and im(d_4^T)
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

    if ob_T4 is not None and ob_T4.shape[0] > 0:
        im_total = bd_T4 @ ob_T4.T % PRIME
    else:
        im_total = np.zeros((len(paths_T_3), 1), dtype=np.int64)

    rank_old = _gauss_rank_np(im_old.copy() % PRIME, PRIME)
    rank_total = _gauss_rank_np(im_total.copy() % PRIME, PRIME)

    # Targeting test
    combined_old_cycles = np.hstack([im_old, old_cycles_Apath]) % PRIME
    rank_old_and_cycles = _gauss_rank_np(combined_old_cycles, PRIME)
    combined_total_cycles = np.hstack([im_total, old_cycles_Apath]) % PRIME
    rank_total_and_cycles = _gauss_rank_np(combined_total_cycles, PRIME)

    old_killed_by_old = old_cycles_Apath.shape[1] - (rank_old_and_cycles - rank_old)
    old_killed_by_total = old_cycles_Apath.shape[1] - (rank_total_and_cycles - rank_total)
    new_kills_extra_old = old_killed_by_total - old_killed_by_old

    # rank(i_*): number of old cycles that survive in H_3(T)
    rank_istar = rank_total_and_cycles - rank_total

    return {
        'b3_T': b3_T, 'b3_Tv': b3_Tv, 'b4_T': b4_T, 'b4_Tv': b4_Tv,
        'ker_dim_Tv': ker_d3_Tv.shape[0],
        'rank_old': rank_old, 'rank_total': rank_total,
        'old_killed_by_old': old_killed_by_old,
        'old_killed_by_total': old_killed_by_total,
        'new_kills_extra_old': new_kills_extra_old,
        'rank_istar': rank_istar,
    }


def main():
    print("=" * 70)
    print("i_*-INJECTIVITY FAILURE ANALYSIS AT n=8")
    print("=" * 70)

    n = 8
    # Use seed 12345 which kind-pasteur used to find failures
    rng = np.random.RandomState(12345)

    found = 0
    target = 2000  # need ~2000 to find ~5 failures
    t0 = time.time()
    results_injective = []
    results_failure = []
    checked = 0

    while checked < target:
        A = random_tournament(n, rng)
        checked += 1
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            r = analyze_failure(A, n, v)
            if r is not None:
                if r['rank_istar'] == 0:
                    results_failure.append(r)
                    scores = sorted(sum(A[i]) for i in range(n))
                    print(f"  FAILURE #{len(results_failure)}: trial={checked}, v={v}, "
                          f"b4={r['b4_T']}, extra_kills={r['new_kills_extra_old']}, "
                          f"scores={tuple(scores)}")
                else:
                    results_injective.append(r)

        if checked % 500 == 0:
            elapsed = time.time() - t0
            print(f"  [{checked}/{target}, b3=1: {found}, failures: {len(results_failure)}, "
                  f"injective: {len(results_injective)}, {elapsed:.1f}s]")

    elapsed = time.time() - t0
    print(f"\n  Total: {checked} checked, {found} b3=1, {elapsed:.1f}s")

    print(f"\n{'='*70}")
    print(f"RESULTS:")
    print(f"  Injective (rank_i*=1): {len(results_injective)}")
    print(f"  Failure (rank_i*=0): {len(results_failure)}")

    if results_failure:
        print(f"\n  FAILURE ANALYSIS:")
        extra_dist = Counter(r['new_kills_extra_old'] for r in results_failure)
        print(f"    extra_kills distribution: {dict(sorted(extra_dist.items()))}")

        b4_dist = Counter(r['b4_T'] for r in results_failure)
        print(f"    b4(T) distribution: {dict(sorted(b4_dist.items()))}")

        b4v_dist = Counter(r['b4_Tv'] for r in results_failure)
        print(f"    b4(T\\v) distribution: {dict(sorted(b4v_dist.items()))}")

        for r in results_failure:
            print(f"    b4_T={r['b4_T']}, b4_Tv={r['b4_Tv']}, "
                  f"rank_old={r['rank_old']}, rank_total={r['rank_total']}, "
                  f"killed_by_old={r['old_killed_by_old']}, "
                  f"killed_by_total={r['old_killed_by_total']}, "
                  f"extra={r['new_kills_extra_old']}")

    if results_injective:
        print(f"\n  INJECTIVE ANALYSIS:")
        extra_dist = Counter(r['new_kills_extra_old'] for r in results_injective)
        print(f"    extra_kills distribution: {dict(sorted(extra_dist.items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
