"""
failure_structure_n8.py — Structural analysis of rank(i_*)=0 BAD vertices at n=8

Recreates the exact RNG sequence from istar_failure_n8.py (seed 12345)
to find the 5 failure tournaments and analyze them.

Key question: what structural property makes i_* fail?

Author: opus-2026-03-09-S57
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


def compute_istar_rank_detailed(A, n, v):
    """Compute rank(i_*: H_3(T\\v) → H_3(T)) and structural details."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)
    b4_T = cc_T['bettis'].get(4, 0)
    b4_Tv = cc_Tv['bettis'].get(4, 0)

    if b3_T < 1 or b3_Tv < 1:
        return None  # Need both nonzero for meaningful i_*

    # Compute i_* rank via embedding ker(d_3^{T\v}) into H_3(T)
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_4_Tv = ap_Tv.get(4, [])

    # Ω_3(T) and Ω_3(T\v) bases
    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob3_T = get_omega(ap_T, 3)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_T = get_omega(ap_T, 4)
    ob4_Tv = get_omega(ap_Tv, 4)

    # d_3(T\v) to get ker(d_3^{T\v}) in Ω_3 coords
    paths_2_Tv = ap_Tv.get(2, [])
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, ker_vecs = _gauss_nullbasis_modp(
            d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3Tv = np.array(ker_vecs, dtype=np.int64) if ker_vecs else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64) if ob3_Tv.shape[0] > 0 else np.zeros((0, 0), dtype=np.int64)

    # ker in A_3(T\v) coords
    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Embed into A_3(T)
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # im(d_4^T) in A_3
    idx3T_d = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T_d:
                bd4_T[idx3T_d[face], j] = (bd4_T[idx3T_d[face], j] + sign) % PRIME

    if ob4_T.shape[0] > 0:
        im_d4 = (bd4_T @ ob4_T.T % PRIME).T
    else:
        im_d4 = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    if im_d4.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4

    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4

    return {
        'b3_T': b3_T, 'b3_Tv': b3_Tv, 'b4_T': b4_T, 'b4_Tv': b4_Tv,
        'rank_istar': rank_istar,
        'rk_d4': rk_d4,
        'dim_omega3_T': ob3_T.shape[0],
        'dim_omega3_Tv': ob3_Tv.shape[0],
        'dim_omega4_T': ob4_T.shape[0],
        'dim_omega4_Tv': ob4_Tv.shape[0] if ob4_Tv.size > 0 else 0,
        'dim_ker_d3Tv': ker_d3Tv.shape[0],
        'n_paths3_T': len(paths_3_T),
        'n_paths4_T': len(paths_4_T),
        'omega_dims_T': {p: cc_T.get('omega_dims', {}).get(p, 0) for p in range(n)},
        'omega_dims_Tv': {p: cc_Tv.get('omega_dims', {}).get(p, 0) for p in range(n - 1)},
        'bettis_T': {p: cc_T['bettis'].get(p, 0) for p in range(n)},
        'bettis_Tv': {p: cc_Tv['bettis'].get(p, 0) for p in range(n - 1)},
        'scores': sorted([int(sum(A[i])) for i in range(n)]),
    }


def main():
    print("=" * 70)
    print("STRUCTURAL ANALYSIS OF rank(i_*)=0 FAILURES AT n=8")
    print("Matching istar_failure_n8.py exactly (seed 12345)")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(12345)

    checked = 0
    found_b3_1 = 0
    failures = []
    successes = []  # first few for comparison
    t0 = time.time()

    while checked < 2000:
        A = random_tournament(n, rng)
        checked += 1

        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found_b3_1 += 1

        for v in range(n):
            r = compute_istar_rank_detailed(A, n, v)
            if r is None:
                continue

            if r['rank_istar'] == 0:
                failures.append((checked, v, A.copy(), r))
                print(f"  FAILURE #{len(failures)}: checked={checked}, v={v}, "
                      f"scores={r['scores']}, b3_T={r['b3_T']}, b3_Tv={r['b3_Tv']}")
            elif len(successes) < 10:
                successes.append((checked, v, A.copy(), r))

    elapsed = time.time() - t0
    print(f"\n  Total: {checked} checked, {found_b3_1} b3=1, "
          f"{len(failures)} failures, {len(successes)} success samples, {elapsed:.1f}s")

    # Detailed analysis of failures
    print(f"\n{'='*60}")
    print("FAILURE DETAILS")
    print(f"{'='*60}")

    for i, (trial, v, A, r) in enumerate(failures):
        print(f"\n  Failure #{i+1}: trial={trial}, v={v}")
        print(f"    scores = {r['scores']}")
        print(f"    bettis(T) = {[r['bettis_T'].get(p, 0) for p in range(n)]}")
        print(f"    bettis(T\\v) = {[r['bettis_Tv'].get(p, 0) for p in range(n-1)]}")
        print(f"    Ω(T) = {[r['omega_dims_T'].get(p, 0) for p in range(n)]}")
        print(f"    Ω(T\\v) = {[r['omega_dims_Tv'].get(p, 0) for p in range(n-1)]}")
        print(f"    Ω_3(T)={r['dim_omega3_T']}, Ω_3(T\\v)={r['dim_omega3_Tv']}, "
              f"Ω_4(T)={r['dim_omega4_T']}, Ω_4(T\\v)={r['dim_omega4_Tv']}")
        print(f"    rank(i_*)={r['rank_istar']}, rk(d_4)={r['rk_d4']}, "
              f"dim_ker_d3Tv={r['dim_ker_d3Tv']}")

        out_deg = int(sum(A[v]))
        print(f"    out_deg(v)={out_deg}")

        # Check ALL vertices of this tournament
        print(f"    All BAD vertices in this tournament:")
        for u in range(n):
            r2 = compute_istar_rank_detailed(A, n, u)
            if r2 is not None:
                out_u = int(sum(A[u]))
                tag = "FAIL" if r2['rank_istar'] == 0 else "ok"
                print(f"      v={u} (out={out_u}): b3_Tv={r2['b3_Tv']}, rank(i_*)={r2['rank_istar']} [{tag}]")

    # Detailed analysis of successes for comparison
    print(f"\n{'='*60}")
    print("SUCCESS EXAMPLES (for comparison)")
    print(f"{'='*60}")

    for i, (trial, v, A, r) in enumerate(successes[:5]):
        print(f"\n  Success #{i+1}: trial={trial}, v={v}")
        print(f"    scores = {r['scores']}")
        print(f"    bettis(T) = {[r['bettis_T'].get(p, 0) for p in range(n)]}")
        print(f"    Ω_3(T)={r['dim_omega3_T']}, Ω_3(T\\v)={r['dim_omega3_Tv']}, "
              f"Ω_4(T)={r['dim_omega4_T']}")
        print(f"    rank(i_*)={r['rank_istar']}, rk(d_4)={r['rk_d4']}")

    # Statistical comparison
    print(f"\n{'='*60}")
    print("STATISTICAL COMPARISON")
    print(f"{'='*60}")

    if failures and successes:
        fail_omega3 = [r['dim_omega3_T'] for _, _, _, r in failures]
        succ_omega3 = [r['dim_omega3_T'] for _, _, _, r in successes]
        fail_omega4 = [r['dim_omega4_T'] for _, _, _, r in failures]
        succ_omega4 = [r['dim_omega4_T'] for _, _, _, r in successes]
        fail_ratio = [r['dim_omega4_T'] / r['dim_omega3_T'] if r['dim_omega3_T'] > 0 else 0
                      for _, _, _, r in failures]
        succ_ratio = [r['dim_omega4_T'] / r['dim_omega3_T'] if r['dim_omega3_T'] > 0 else 0
                      for _, _, _, r in successes]

        print(f"  Ω_3(T) — fail: {fail_omega3}, mean={np.mean(fail_omega3):.1f}")
        print(f"  Ω_3(T) — succ: mean={np.mean(succ_omega3):.1f} (n={len(succ_omega3)})")
        print(f"  Ω_4(T) — fail: {fail_omega4}, mean={np.mean(fail_omega4):.1f}")
        print(f"  Ω_4(T) — succ: mean={np.mean(succ_omega4):.1f}")
        print(f"  Ω_4/Ω_3 — fail: {[f'{x:.2f}' for x in fail_ratio]}")
        print(f"  Ω_4/Ω_3 — succ: mean={np.mean(succ_ratio):.2f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
