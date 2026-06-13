"""
h3_projection_vs_istar.py — Relate H_3(T) generator old-projection to rank(i_*)

The key relationship:
  codim(im(d_4)_old, ker(d_3)_old) = 1 (HYP-408)
means that ker(d_3)_old / im(d_4)_old is 1-dimensional.

The question is: does the embedded H_3(T\\v) gen span this 1-dim quotient?
- YES => rank(i_*) = 1
- NO => rank(i_*) = 0

But we can also ask: does the H_3(T) gen's old-projection span this quotient?
By definition it does (the generator IS the quotient).

So the real question: is the embedded H_3(T\\v) gen's old-projection
PROPORTIONAL to the H_3(T) gen's old-projection (mod im(d_4)_old)?

If yes: rank(i_*) = 1.
If no (emb gen is in im(d_4)_old while H_3 gen is not): rank(i_*) = 0.

This is exactly what rank(i_*) measures. But we can look at it from a
different angle: compute the actual projection coefficient.

Author: opus-2026-03-09-S58
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


def projection_coefficient(A, n, v):
    """Compute the projection coefficient of emb H_3(T\\v) gen onto H_3(T) gen."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
    paths_4_Tv = ap_Tv.get(4, [])
    paths_2_T = ap_T.get(2, [])
    paths_2_Tv = ap_Tv.get(2, [])

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
    ob4_T = get_omega(ap_T, 4)
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_Tv = get_omega(ap_Tv, 4)

    # ker(d_3^T) and im(d_4^T) in A_3(T)
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    idx3T = {p: i for i, p in enumerate(paths_3_T)}

    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2T:
                bd3_T[idx2T[face], j] = (bd3_T[idx2T[face], j] + sign) % PRIME

    d3o_T = bd3_T @ ob3_T.T % PRIME
    _, kv_T = _gauss_nullbasis_modp(d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME) if d3o_T.size > 0 else (0, [])
    ker_d3T = (np.array(kv_T, dtype=np.int64) @ ob3_T % PRIME) if kv_T else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    im_d4T = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4T = int(_gauss_rank_np(im_d4T.copy(), PRIME)) if im_d4T.shape[0] > 0 else 0

    # H_3(T) generator
    h3T_gen = None
    for i in range(ker_d3T.shape[0]):
        if im_d4T.shape[0] > 0:
            test = np.vstack([im_d4T, ker_d3T[i:i+1]]) % PRIME
            if int(_gauss_rank_np(test.copy(), PRIME)) > rk_d4T:
                h3T_gen = ker_d3T[i]
                break
        else:
            h3T_gen = ker_d3T[i]
            break

    if h3T_gen is None:
        return None

    # ker(d_3^{T\\v}) and im(d_4^{T\\v})
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}

    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME) if d3o_Tv.size > 0 else (0, [])
    ker_d3Tv = (np.array(kv_Tv, dtype=np.int64) @ ob3_Tv % PRIME) if kv_Tv else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    bd4_Tv = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3Tv:
                bd4_Tv[idx3Tv[face], j] = (bd4_Tv[idx3Tv[face], j] + sign) % PRIME

    im_d4Tv = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
    rk_d4Tv = int(_gauss_rank_np(im_d4Tv.copy(), PRIME)) if im_d4Tv.shape[0] > 0 else 0

    # H_3(T\\v) generator
    h3Tv_gen = None
    for i in range(ker_d3Tv.shape[0]):
        if im_d4Tv.shape[0] > 0:
            test = np.vstack([im_d4Tv, ker_d3Tv[i:i+1]]) % PRIME
            if int(_gauss_rank_np(test.copy(), PRIME)) > rk_d4Tv:
                h3Tv_gen = ker_d3Tv[i]
                break
        else:
            h3Tv_gen = ker_d3Tv[i]
            break

    if h3Tv_gen is None:
        return None

    # Embed H_3(T\\v) gen into A_3(T)
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    emb_gen = h3Tv_gen @ embed % PRIME  # 1D vector in A_3(T)

    # rank(i_*)
    if im_d4T.shape[0] > 0:
        combined = np.vstack([im_d4T, emb_gen.reshape(1, -1)]) % PRIME
        rk_comb = int(_gauss_rank_np(combined.copy(), PRIME))
        rank_istar = rk_comb - rk_d4T
    else:
        rank_istar = 1 if np.any(emb_gen % PRIME != 0) else 0

    # Support analysis
    h3T_supp = set(i for i in range(len(h3T_gen)) if h3T_gen[i] % PRIME != 0)
    emb_supp = set(i for i in range(len(emb_gen)) if emb_gen[i] % PRIME != 0)

    # Vertices in H_3(T) support paths
    h3T_vertices = set()
    for i in h3T_supp:
        h3T_vertices.update(paths_3_T[i])

    # Per-vertex: old support size
    old_support_by_v = {}
    for u in range(n):
        old_paths = [i for i in h3T_supp if u not in paths_3_T[i]]
        old_support_by_v[u] = len(old_paths)

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'h3T_support_size': len(h3T_supp),
        'emb_support_size': len(emb_supp),
        'h3T_old_support_for_v': old_support_by_v[v],
        'h3T_min_old_support': min(old_support_by_v.values()),
        'h3T_through_v_frac': 1 - old_support_by_v[v] / max(len(h3T_supp), 1),
        'n_h3T_vertices': len(h3T_vertices),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"PROJECTION COEFFICIENT ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(12345 if n == 8 else 42)
        results = []
        target_checked = 2000 if n == 8 else 500
        checked = 0
        t0 = time.time()

        for trial in range(50000):
            if checked >= target_checked:
                break
            A = random_tournament(n, rng)
            checked += 1

            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = projection_coefficient(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

            if len(results) % 200 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)} BAD vertices, {elapsed:.1f}s")

        elapsed = time.time() - t0
        fails = [r for r in results if r['rank_istar'] == 0]
        succs = [r for r in results if r['rank_istar'] == 1]
        print(f"  {len(results)} BAD vertices ({len(fails)} fail), {elapsed:.1f}s")

        # H_3(T) old support for v (the deleted vertex)
        print(f"\n  H_3(T) old support for deleted vertex v:")
        if fails:
            fail_old = [r['h3T_old_support_for_v'] for r in fails]
            print(f"    FAILURES: {fail_old[:20]}")
        if succs:
            succ_old = [r['h3T_old_support_for_v'] for r in succs]
            print(f"    SUCCESSES: min={min(succ_old)}, max={max(succ_old)}, mean={np.mean(succ_old):.1f}")

        # Through-v fraction
        print(f"\n  Fraction of H_3(T) support through v:")
        if fails:
            fail_frac = [r['h3T_through_v_frac'] for r in fails]
            print(f"    FAILURES: {[f'{x:.3f}' for x in fail_frac[:10]]}")
        if succs:
            succ_frac = [r['h3T_through_v_frac'] for r in succs]
            print(f"    SUCCESSES: min={min(succ_frac):.3f}, max={max(succ_frac):.3f}, "
                  f"mean={np.mean(succ_frac):.3f}")

        # Min old support across ALL vertices
        print(f"\n  Min old support across all vertices:")
        min_olds = [r['h3T_min_old_support'] for r in results]
        print(f"    min={min(min_olds)}, max={max(min_olds)}, mean={np.mean(min_olds):.1f}")
        print(f"    distribution: {dict(sorted(Counter(min_olds).items()))}")

        # Correlation: does low old_support_for_v predict failure?
        if n == 8:
            by_old_support = {}
            for r in results:
                os = r['h3T_old_support_for_v']
                if os not in by_old_support:
                    by_old_support[os] = {'fail': 0, 'succ': 0}
                if r['rank_istar'] == 0:
                    by_old_support[os]['fail'] += 1
                else:
                    by_old_support[os]['succ'] += 1

            print(f"\n  Failure rate by H_3(T) old support for v:")
            for os in sorted(by_old_support.keys()):
                d = by_old_support[os]
                total = d['fail'] + d['succ']
                rate = 100 * d['fail'] / total if total > 0 else 0
                if total >= 5:
                    print(f"    old_support={os}: {d['fail']}/{total} ({rate:.2f}%)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
