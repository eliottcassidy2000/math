"""
targeting_n8.py — Test HYP-398 (new→new targeting) at n=8

At n=7, new d_4 boundaries NEVER kill old cycles (34/34 BAD vertices).
Question: Does this still hold at n=8?
- For beta_3=1 tournaments (should still hold since i_*-injectivity holds for MOST)
- For beta_3=2 tournaments (likely FAILS — this is where the mechanism breaks)

Also test: when beta_3(T\v)=beta_3(T) for BAD vertices, does the mechanism hold?
When beta_3(T\v) < beta_3(T) for GOOD vertices, what changes?

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


def check_targeting_general(A, n, v):
    """General targeting check — works for any beta_3 values."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    if b3_Tv == 0:
        return None  # Only care about BAD vertices (b3_Tv > 0)

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

    # ker(d_3^{T\v}) in Omega coords
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

    ker_dim_Tv = ker_d3_Tv.shape[0]

    # Embed ker_d3(T\v) into A_3(T) space
    embed_3 = np.zeros((len(paths_T_3), len(paths_Tv_3)), dtype=np.int64)
    for j, ptv in enumerate(paths_Tv_3):
        emb = embed_path(ptv)
        if emb in idx_T3:
            embed_3[idx_T3[emb], j] = 1
    embedded_Tv3 = embed_3 @ ob_Tv3.T % PRIME

    old_cycles_Apath = embedded_Tv3 @ ker_d3_Tv.T % PRIME

    # d_4 boundary matrix in A_3(T)
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    # im(d_4^old) from T\v 5-paths embedded
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

    combined_old_cycles = np.hstack([im_old, old_cycles_Apath]) % PRIME
    rank_old_and_cycles = _gauss_rank_np(combined_old_cycles, PRIME)

    combined_total_cycles = np.hstack([im_total, old_cycles_Apath]) % PRIME
    rank_total_and_cycles = _gauss_rank_np(combined_total_cycles, PRIME)

    old_killed_by_old = old_cycles_Apath.shape[1] - (rank_old_and_cycles - rank_old)
    old_killed_by_total = old_cycles_Apath.shape[1] - (rank_total_and_cycles - rank_total)
    new_kills_extra_old = old_killed_by_total - old_killed_by_old

    # Also compute rank(i_*^3)
    # i_* maps ker_d3(T\v) into ker_d3(T) / im(d_4^T)
    # rank(i_*) = rank of [im_total | old_cycles] - rank_total
    rank_istar = rank_total_and_cycles - rank_total

    return {
        'b3_T': b3_T,
        'b3_Tv': b3_Tv,
        'ker_dim_Tv': ker_dim_Tv,
        'rank_old': rank_old,
        'rank_total': rank_total,
        'old_killed_by_old': old_killed_by_old,
        'old_killed_by_total': old_killed_by_total,
        'new_kills_extra_old': new_kills_extra_old,
        'rank_istar': rank_istar,
    }


def main():
    print("=" * 70)
    print("NEW BOUNDARY TARGETING AT n=8")
    print("Testing HYP-398 for beta_3=1 and beta_3=2 tournaments")
    print("=" * 70)

    n = 8
    target_b3_1 = 30   # beta_3=1 tournaments
    target_b3_2 = 5    # beta_3=2 tournaments (rare ~0.08%)
    rng = np.random.RandomState(42)

    t0 = time.time()
    found_b3_1 = 0
    found_b3_2 = 0
    checked = 0
    results_b3_1 = []
    results_b3_2 = []

    while found_b3_1 < target_b3_1 or found_b3_2 < target_b3_2:
        A = random_tournament(n, rng)
        checked += 1
        cc = full_chain_complex_modp(A, n, n - 1)
        b3 = cc['bettis'].get(3, 0)

        if b3 == 0:
            continue

        if b3 == 1 and found_b3_1 >= target_b3_1:
            continue
        if b3 == 2 and found_b3_2 >= target_b3_2:
            continue

        scores = sorted(sum(A[i]) for i in range(n))
        bettis = tuple(cc['bettis'].get(p, 0) for p in range(n))

        if b3 == 1:
            found_b3_1 += 1
            label = f"b3=1 #{found_b3_1}"
        elif b3 == 2:
            found_b3_2 += 1
            label = f"b3=2 #{found_b3_2}"
        else:
            label = f"b3={b3}"

        print(f"\n  {label}: bettis={bettis}, scores={tuple(scores)}")

        for v in range(n):
            r = check_targeting_general(A, n, v)
            if r is not None:
                if b3 == 1:
                    results_b3_1.append(r)
                elif b3 == 2:
                    results_b3_2.append(r)
                extra = r['new_kills_extra_old']
                tag = "OK" if extra == 0 else f"EXTRA={extra}"
                print(f"    v={v}: b3_Tv={r['b3_Tv']}, rank_i*={r['rank_istar']}, "
                      f"extra_kills={extra} [{tag}]")

        if checked % 500 == 0:
            elapsed = time.time() - t0
            print(f"  [{checked} checked, b3=1: {found_b3_1}/{target_b3_1}, "
                  f"b3=2: {found_b3_2}/{target_b3_2}, {elapsed:.1f}s]")

    elapsed = time.time() - t0
    print(f"\n  Total: {checked} checked, {elapsed:.1f}s")

    # Summary for beta_3=1
    print(f"\n{'='*70}")
    print(f"RESULTS FOR beta_3=1 ({len(results_b3_1)} BAD vertices):")
    if results_b3_1:
        extra_dist = Counter(r['new_kills_extra_old'] for r in results_b3_1)
        print(f"  Extra old kills: {dict(sorted(extra_dist.items()))}")
        istar_dist = Counter(r['rank_istar'] for r in results_b3_1)
        print(f"  rank(i_*): {dict(sorted(istar_dist.items()))}")
        if all(r['new_kills_extra_old'] == 0 for r in results_b3_1):
            print(f"  *** HYP-398 CONFIRMED at n=8 for beta_3=1 ***")
        else:
            fail = sum(1 for r in results_b3_1 if r['new_kills_extra_old'] != 0)
            print(f"  HYP-398 FAILS: {fail}/{len(results_b3_1)} cases")

    # Summary for beta_3=2
    print(f"\nRESULTS FOR beta_3=2 ({len(results_b3_2)} BAD vertices):")
    if results_b3_2:
        extra_dist = Counter(r['new_kills_extra_old'] for r in results_b3_2)
        print(f"  Extra old kills: {dict(sorted(extra_dist.items()))}")
        istar_dist = Counter(r['rank_istar'] for r in results_b3_2)
        print(f"  rank(i_*): {dict(sorted(istar_dist.items()))}")

        # Detailed: for b3=2 tours, what are b3_Tv values?
        b3tv_dist = Counter(r['b3_Tv'] for r in results_b3_2)
        print(f"  b3(T\\v) values: {dict(sorted(b3tv_dist.items()))}")

        if all(r['new_kills_extra_old'] == 0 for r in results_b3_2):
            print(f"  *** HYP-398 HOLDS even for beta_3=2! ***")
        else:
            fail = sum(1 for r in results_b3_2 if r['new_kills_extra_old'] != 0)
            print(f"  HYP-398 FAILS for beta_3=2: {fail}/{len(results_b3_2)}")
            # Show which ones fail
            for r in results_b3_2:
                if r['new_kills_extra_old'] != 0:
                    print(f"    b3_T={r['b3_T']}, b3_Tv={r['b3_Tv']}, "
                          f"rank_i*={r['rank_istar']}, extra={r['new_kills_extra_old']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
