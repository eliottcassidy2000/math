"""
failure_predictor_n8.py — What predicts rank(i_*)=0 at n=8?

We know:
  - Overall rate: ~0.85% of BAD vertices
  - Out-degree correlation: symmetric, peaks at |d - 3.5| = 1.5
  - Failure vertices have H_3(T) gen concentrated on through-v paths
  - Score sequences: concentrated on a few types

This script collects a large sample and tests multiple predictors:
1. Out-degree of v
2. Score sequence of T
3. Omega_3 / Omega_4 ratio
4. Number of 4-paths through v
5. b_6(T) (higher Betti)

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
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


def quick_istar(A, n, v, cc_T=None):
    """Compute rank(i_*) quickly."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    if cc_T is None:
        cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    paths_3_T = ap_T.get(3, [])
    paths_4_T = ap_T.get(4, [])
    paths_3_Tv = ap_Tv.get(3, [])
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

    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    idx3T = {p: i for i, p in enumerate(paths_3_T)}

    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3Tv = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    im_d4 = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    if im_d4.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4
    rk_comb = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_comb - rk_d4

    out_deg = int(sum(A[v]))
    n_4_through_v = sum(1 for p in paths_4_T if v in p)

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'n_4_through_v': n_4_through_v,
        'n_3_through_v': sum(1 for p in paths_3_T if v in p),
        'omega3_T': ob3_T.shape[0],
        'omega4_T': ob4_T.shape[0],
        'omega3_Tv': ob3_Tv.shape[0],
        'dim_ker_d3Tv': ker_d3Tv.shape[0],
    }


def main():
    n = 8
    print(f"{'='*70}")
    print(f"FAILURE PREDICTOR ANALYSIS AT n={n}")
    print(f"{'='*70}")

    rng = np.random.RandomState(777)
    results = []
    tournament_data = {}  # keyed by tournament index
    checked = 0
    t_idx = 0
    t0 = time.time()

    while checked < 3000:
        A = random_tournament(n, rng)
        checked += 1

        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue

        scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
        b6 = cc['bettis'].get(6, 0)
        omega_dims = {p: cc.get('omega_dims', {}).get(p, 0) for p in range(n)}

        for v_cand in range(n):
            r = quick_istar(A, n, v_cand, cc_T=cc)
            if r is None:
                continue
            r['scores'] = scores
            r['b6'] = b6
            r['omega_dims'] = omega_dims
            r['t_idx'] = t_idx
            results.append(r)

        t_idx += 1

        if len(results) % 500 == 0 and len(results) > 0:
            elapsed = time.time() - t0
            n_fail = sum(1 for r in results if r['rank_istar'] == 0)
            print(f"  Progress: {len(results)} BAD, {n_fail} fail, {elapsed:.1f}s")

    elapsed = time.time() - t0
    fails = [r for r in results if r['rank_istar'] == 0]
    succs = [r for r in results if r['rank_istar'] == 1]
    print(f"\n  TOTAL: {len(results)} BAD vertices, {len(fails)} failures "
          f"({100*len(fails)/max(len(results),1):.2f}%), {elapsed:.1f}s")

    # 1. Out-degree
    print(f"\n{'='*50}")
    print("PREDICTOR 1: Out-degree")
    print(f"{'='*50}")
    by_od = defaultdict(lambda: [0, 0])
    for r in results:
        by_od[r['out_deg']][r['rank_istar']] += 1
    for od in sorted(by_od.keys()):
        total = sum(by_od[od])
        fail = by_od[od][0]
        print(f"  out_deg={od}: {fail}/{total} fail ({100*fail/total:.2f}%)")

    # 2. Score sequence
    print(f"\n{'='*50}")
    print("PREDICTOR 2: Score sequence")
    print(f"{'='*50}")
    by_score = defaultdict(lambda: [0, 0])
    for r in results:
        by_score[r['scores']][r['rank_istar']] += 1
    for scores in sorted(by_score.keys(), key=lambda s: -sum(by_score[s])):
        total = sum(by_score[scores])
        if total < 10:
            continue
        fail = by_score[scores][0]
        rate = 100 * fail / total if total > 0 else 0
        print(f"  {list(scores)}: {fail}/{total} ({rate:.1f}%)")

    # 3. b_6(T)
    print(f"\n{'='*50}")
    print("PREDICTOR 3: b_6(T)")
    print(f"{'='*50}")
    by_b6 = defaultdict(lambda: [0, 0])
    for r in results:
        by_b6[r['b6']][r['rank_istar']] += 1
    for b6 in sorted(by_b6.keys()):
        total = sum(by_b6[b6])
        if total < 5:
            continue
        fail = by_b6[b6][0]
        rate = 100 * fail / total if total > 0 else 0
        print(f"  b6={b6}: {fail}/{total} ({rate:.1f}%)")

    # 4. Omega ratio
    print(f"\n{'='*50}")
    print("PREDICTOR 4: Omega_4/Omega_3 ratio")
    print(f"{'='*50}")
    if fails:
        fail_ratio = [r['omega4_T'] / max(r['omega3_T'], 1) for r in fails]
        succ_ratio = [r['omega4_T'] / max(r['omega3_T'], 1) for r in succs[:200]]
        print(f"  Failures: min={min(fail_ratio):.3f}, max={max(fail_ratio):.3f}, mean={np.mean(fail_ratio):.3f}")
        print(f"  Successes: min={min(succ_ratio):.3f}, max={max(succ_ratio):.3f}, mean={np.mean(succ_ratio):.3f}")

    # 5. n_4_through_v ratio
    print(f"\n{'='*50}")
    print("PREDICTOR 5: 4-paths through v / total Omega_4")
    print(f"{'='*50}")
    if fails:
        fail_4ratio = [r['n_4_through_v'] / max(r['omega4_T'], 1) for r in fails]
        succ_4ratio = [r['n_4_through_v'] / max(r['omega4_T'], 1) for r in succs[:200]]
        print(f"  Failures: min={min(fail_4ratio):.3f}, max={max(fail_4ratio):.3f}, mean={np.mean(fail_4ratio):.3f}")
        print(f"  Successes: min={min(succ_4ratio):.3f}, max={max(succ_4ratio):.3f}, mean={np.mean(succ_4ratio):.3f}")

    # 6. Combined: score + outdeg
    print(f"\n{'='*50}")
    print("PREDICTOR 6: Score × out-degree interaction")
    print(f"{'='*50}")
    by_combo = defaultdict(lambda: [0, 0])
    for r in results:
        key = (r['scores'], r['out_deg'])
        by_combo[key][r['rank_istar']] += 1
    for key in sorted(by_combo.keys(), key=lambda k: -by_combo[k][0]):
        total = sum(by_combo[key])
        fail = by_combo[key][0]
        if fail == 0:
            continue
        scores, od = key
        rate = 100 * fail / total
        print(f"  scores={list(scores)}, out_deg={od}: {fail}/{total} ({rate:.1f}%)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
