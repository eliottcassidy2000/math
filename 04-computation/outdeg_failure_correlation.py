"""
outdeg_failure_correlation.py — Test whether low out-degree correlates with i_* failure

At n=8, rank(i_*)=0 failures occur ~0.65% of BAD vertices.
From the 5 known failures: 3/5 have out_deg=2, 1 has out_deg=3, 1 has out_deg=4.
Question: Is low out-degree a significant predictor of failure?

Also tests: does the embedded H_3(T\v) generator align with H_3(T) generator
in successes but not in failures? (HYP-399 direction mismatch theory)

Author: opus-2026-03-09-S57
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


def compute_istar_and_direction(A, n, v):
    """Compute rank(i_*) and direction alignment info."""
    max_p = min(n - 1, 6)
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    b3_T = cc_T['bettis'].get(3, 0)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    if b3_T != 1 or b3_Tv != 1:
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
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_T = get_omega(ap_T, 4)

    # d_3(T\v): Omega_3(T\v) -> A_2(T\v)
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

    # H_3(T\v) generator in A_3(T\v) coords (should be 1-dim since b3=1)
    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Embed into A_3(T)
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # im(d_4^T) in A_3(T)
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

    # DIRECTION ANALYSIS: How many of the embedded ker vectors are "new" vs in im(d_4)?
    # For each row of embedded_ker, check if it's in span(im_d4)
    n_embedded = embedded_ker.shape[0]
    n_in_imd4 = 0
    if n_embedded > 0 and im_d4.shape[0] > 0:
        for row_idx in range(n_embedded):
            row = embedded_ker[row_idx:row_idx+1]
            test = np.vstack([im_d4, row]) % PRIME
            rk_test = int(_gauss_rank_np(test.copy(), PRIME))
            if rk_test == rk_d4:
                n_in_imd4 += 1

    out_deg = int(sum(A[v]))
    in_deg = n - 1 - out_deg

    # Score of v relative to others
    all_scores = sorted([int(sum(A[i])) for i in range(n)])
    v_score = out_deg

    # Count how many 3-paths go through v
    n_paths_through_v = sum(1 for p in paths_3_T if v in p)
    n_paths_total = len(paths_3_T)

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'in_deg': in_deg,
        'v_score': v_score,
        'scores': all_scores,
        'rk_d4': rk_d4,
        'dim_ker_d3Tv': ker_d3Tv.shape[0],
        'n_embedded': n_embedded,
        'n_in_imd4': n_in_imd4,
        'n_paths_through_v': n_paths_through_v,
        'n_paths_total': n_paths_total,
        'omega3_T': ob3_T.shape[0],
        'omega4_T': ob4_T.shape[0],
    }


def main():
    print("=" * 70)
    print("OUT-DEGREE vs i_* FAILURE CORRELATION AT n=8")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(999)  # different seed for independent sample

    # Collect data by out-degree
    by_outdeg = defaultdict(lambda: {'fail': 0, 'success': 0, 'details': []})
    total_checked = 0
    total_b3_1 = 0
    total_bad = 0
    total_fail = 0
    t0 = time.time()

    target_bad = 2000  # collect 2000 BAD vertices

    for trial in range(50000):
        if total_bad >= target_bad:
            break
        A = random_tournament(n, rng)
        total_checked += 1

        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        total_b3_1 += 1

        for v in range(n):
            r = compute_istar_and_direction(A, n, v)
            if r is None:
                continue

            total_bad += 1
            od = r['out_deg']
            if r['rank_istar'] == 0:
                by_outdeg[od]['fail'] += 1
                total_fail += 1
                by_outdeg[od]['details'].append(r)
            else:
                by_outdeg[od]['success'] += 1

        if total_bad % 200 == 0 and total_bad > 0:
            elapsed = time.time() - t0
            print(f"  Progress: {total_bad} BAD vertices, {total_fail} failures, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  TOTAL: {total_checked} tournaments, {total_b3_1} with b3=1, "
          f"{total_bad} BAD vertices, {total_fail} failures ({100*total_fail/max(total_bad,1):.2f}%), {elapsed:.1f}s")

    # Out-degree breakdown
    print(f"\n{'='*60}")
    print("FAILURE RATE BY OUT-DEGREE")
    print(f"{'='*60}")
    for od in sorted(by_outdeg.keys()):
        d = by_outdeg[od]
        total_od = d['fail'] + d['success']
        rate = 100 * d['fail'] / total_od if total_od > 0 else 0
        print(f"  out_deg={od}: {d['fail']}/{total_od} fail ({rate:.2f}%)")

    # In-degree analysis
    print(f"\n{'='*60}")
    print("FAILURE RATE BY IN-DEGREE")
    print(f"{'='*60}")
    by_indeg = defaultdict(lambda: {'fail': 0, 'success': 0})
    for od, d in by_outdeg.items():
        ind = n - 1 - od
        by_indeg[ind]['fail'] += d['fail']
        by_indeg[ind]['success'] += d['success']
    for ind in sorted(by_indeg.keys()):
        d = by_indeg[ind]
        total_ind = d['fail'] + d['success']
        rate = 100 * d['fail'] / total_ind if total_ind > 0 else 0
        print(f"  in_deg={ind}: {d['fail']}/{total_ind} fail ({rate:.2f}%)")

    # Direction analysis for failures
    print(f"\n{'='*60}")
    print("DIRECTION ANALYSIS (failure cases)")
    print(f"{'='*60}")
    for od in sorted(by_outdeg.keys()):
        for r in by_outdeg[od]['details'][:3]:
            print(f"  out_deg={od}: n_embedded_ker={r['n_embedded']}, "
                  f"n_in_im_d4={r['n_in_imd4']}, "
                  f"dim_ker_d3Tv={r['dim_ker_d3Tv']}, rk_d4={r['rk_d4']}")

    # Score sequence analysis
    print(f"\n{'='*60}")
    print("SCORE SEQUENCE ANALYSIS")
    print(f"{'='*60}")
    fail_scores = Counter()
    for od, d in by_outdeg.items():
        for r in d['details']:
            fail_scores[tuple(r['scores'])] += 1
    print("  Failure score sequences:")
    for scores, count in fail_scores.most_common(10):
        print(f"    {list(scores)}: {count} failures")

    # Paths through v analysis
    print(f"\n{'='*60}")
    print("PATHS THROUGH v ANALYSIS")
    print(f"{'='*60}")
    fail_pct = []
    succ_pct = []
    for od, d in by_outdeg.items():
        for r in d['details']:
            fail_pct.append(r['n_paths_through_v'] / max(r['n_paths_total'], 1))
    # Need success details too — collect separately
    # (We didn't store success details to save memory, so just report failures)
    if fail_pct:
        print(f"  Failures: paths_through_v / total_paths_3:")
        print(f"    min={min(fail_pct):.3f}, max={max(fail_pct):.3f}, mean={np.mean(fail_pct):.3f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
