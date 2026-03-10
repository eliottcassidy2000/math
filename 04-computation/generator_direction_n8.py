"""
generator_direction_n8.py — Compare H_3(T\v) and H_3(T) generators directly

Theory: When b3(T)=1 and b3(T\v)=1, ker(d_3^T) has dim = rk(d_4)+1.
The H_3(T) generator spans the unique direction in ker(d_3) not in im(d_4).

rank(i_*)=1: embedded H_3(T\v) generator has nonzero component along this direction
rank(i_*)=0: embedded H_3(T\v) generator lies entirely in im(d_4^T)

This script computes:
1. The actual H_3(T) generator direction
2. The embedded H_3(T\v) generator direction
3. Their "alignment" — the projection coefficient
4. Whether the embedded generator decomposes as im(d_4) + alpha * h3_generator

Author: opus-2026-03-09-S57
"""
import sys
import time
import numpy as np
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


def analyze_generator_direction(A, n, v, verbose=False):
    """Full analysis of generator alignment between T\v and T."""
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
    ob3_Tv = get_omega(ap_Tv, 3)
    ob4_T = get_omega(ap_T, 4)

    # === ker(d_3^T) in Omega_3(T) coords ===
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2T:
                bd3_T[idx2T[face], j] = (bd3_T[idx2T[face], j] + sign) % PRIME

    d3o_T = bd3_T @ ob3_T.T % PRIME
    if d3o_T.size > 0:
        _, ker_vecs_T = _gauss_nullbasis_modp(
            d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME)
        ker_d3T = np.array(ker_vecs_T, dtype=np.int64) if ker_vecs_T else np.zeros((0, ob3_T.shape[0]), dtype=np.int64)
    else:
        ker_d3T = np.eye(ob3_T.shape[0], dtype=np.int64)

    # === im(d_4^T) in Omega_3(T) coords ===
    # d_4: A_4(T) -> A_3(T)
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    # im(d_4) in A_3 coords, then project to Omega_3 coords
    if ob4_T.shape[0] > 0:
        im_d4_A = (bd4_T @ ob4_T.T % PRIME).T  # rows = images
    else:
        im_d4_A = np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # We need im(d_4) in Omega_3(T) coords
    # For that, we'd need to invert ob3_T... but ob3_T maps Omega -> A
    # Instead, work entirely in A_3(T) coords

    # ker(d_3^T) in A_3 coords
    ker_d3T_A = ker_d3T @ ob3_T % PRIME if ker_d3T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # im(d_4^T) in A_3 coords — intersect with ker(d_3) to get im(d_4) ∩ ker(d_3)
    # Actually d_3 ∘ d_4 = 0 so im(d_4) ⊂ ker(d_3) always
    rk_d4 = int(_gauss_rank_np(im_d4_A.copy(), PRIME)) if im_d4_A.shape[0] > 0 else 0

    dim_ker_d3T = ker_d3T_A.shape[0]

    if verbose:
        print(f"    dim ker(d_3^T) = {dim_ker_d3T}, rk(d_4^T) = {rk_d4}")
        print(f"    b3(T) = dim_ker - rk_d4 = {dim_ker_d3T - rk_d4}")

    # === The H_3(T) generator ===
    # It's a vector in ker(d_3^T) not in im(d_4^T)
    # Find it by extending im(d_4) with ker(d_3) rows and finding which ker row is new
    if im_d4_A.shape[0] > 0 and ker_d3T_A.shape[0] > 0:
        # Find the ker(d_3) vector that's NOT in span(im_d4)
        h3_gen_idx = -1
        for i in range(ker_d3T_A.shape[0]):
            test = np.vstack([im_d4_A, ker_d3T_A[i:i+1]]) % PRIME
            rk_test = int(_gauss_rank_np(test.copy(), PRIME))
            if rk_test > rk_d4:
                h3_gen_idx = i
                break
        if h3_gen_idx >= 0:
            h3_gen = ker_d3T_A[h3_gen_idx]
        else:
            return None  # shouldn't happen if b3=1
    elif ker_d3T_A.shape[0] > 0:
        h3_gen = ker_d3T_A[0]
        h3_gen_idx = 0
    else:
        return None

    # === ker(d_3^{T\v}) ===
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    if d3o_Tv.size > 0:
        _, ker_vecs_Tv = _gauss_nullbasis_modp(
            d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
        ker_d3Tv = np.array(ker_vecs_Tv, dtype=np.int64) if ker_vecs_Tv else np.zeros((0, ob3_Tv.shape[0]), dtype=np.int64)
    else:
        ker_d3Tv = np.eye(ob3_Tv.shape[0], dtype=np.int64)

    ker_d3Tv_A = ker_d3Tv @ ob3_Tv % PRIME if ker_d3Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    # Embed into A_3(T)
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # === rank(i_*) ===
    if im_d4_A.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4_A, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4_A
    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4

    # === SUPPORT ANALYSIS ===
    # How many nonzero coords does h3_gen have?
    # How many nonzero coords does embedded H_3(T\v) gen have?
    # What's their overlap?
    h3_gen_supp = set(i for i in range(len(h3_gen)) if h3_gen[i] % PRIME != 0)

    # The H_3(T\v) generator in ker(d_3^{T\v}) — take the one that maps to the
    # im(d_4^{T\v}) complement
    # Since b3(T\v)=1 and we have ker_d3Tv, the "generator" spans ker/im
    # For our purposes, any nonzero element of ker_d3Tv_A works
    # Take the first embedded ker vector
    if embedded_ker.shape[0] > 0:
        emb_gen = embedded_ker[0]
        emb_gen_supp = set(i for i in range(len(emb_gen)) if emb_gen[i] % PRIME != 0)
    else:
        return None

    # Overlap
    overlap = h3_gen_supp & emb_gen_supp
    union = h3_gen_supp | emb_gen_supp

    # "Through-v" paths: paths in A_3(T) that contain vertex v
    through_v_indices = set()
    for i, p in enumerate(paths_3_T):
        if v in p:
            through_v_indices.add(i)

    # Support of h3_gen on through-v vs not-through-v
    h3_through_v = h3_gen_supp & through_v_indices
    h3_not_through_v = h3_gen_supp - through_v_indices

    # Support of embedded gen on through-v (should be 0 — T\v has no v paths)
    emb_through_v = emb_gen_supp & through_v_indices

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'rk_d4': rk_d4,
        'dim_ker_d3T': dim_ker_d3T,
        'h3_gen_support_size': len(h3_gen_supp),
        'emb_gen_support_size': len(emb_gen_supp),
        'support_overlap': len(overlap),
        'support_union': len(union),
        'jaccard': len(overlap) / len(union) if union else 0,
        'h3_through_v': len(h3_through_v),
        'h3_not_through_v': len(h3_not_through_v),
        'emb_through_v': len(emb_through_v),  # should be 0
        'n_paths_3': len(paths_3_T),
        'n_through_v': len(through_v_indices),
    }


def main():
    print("=" * 70)
    print("GENERATOR DIRECTION ANALYSIS: H_3(T) vs embedded H_3(T\\v)")
    print("=" * 70)

    n = 8
    rng = np.random.RandomState(12345)  # same seed as failure_structure to get known failures

    failures = []
    successes = []
    checked = 0
    t0 = time.time()

    while checked < 2000:
        A = random_tournament(n, rng)
        checked += 1

        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue

        for v in range(n):
            r = analyze_generator_direction(A, n, v, verbose=False)
            if r is None:
                continue

            if r['rank_istar'] == 0:
                failures.append((checked, v, r))
            elif len(successes) < 50:
                successes.append((checked, v, r))

    elapsed = time.time() - t0
    print(f"  {checked} checked, {len(failures)} failures, {len(successes)} success samples, {elapsed:.1f}s")

    # Failure analysis
    print(f"\n{'='*60}")
    print("FAILURES — generator direction")
    print(f"{'='*60}")
    for i, (trial, v, r) in enumerate(failures):
        print(f"\n  Failure #{i+1}: trial={trial}, v={v}, out_deg={r['out_deg']}")
        print(f"    dim_ker(d3)={r['dim_ker_d3T']}, rk(d4)={r['rk_d4']}")
        print(f"    H3(T) gen support: {r['h3_gen_support_size']} paths")
        print(f"    Emb H3(Tv) gen support: {r['emb_gen_support_size']} paths")
        print(f"    Overlap: {r['support_overlap']}, Union: {r['support_union']}, Jaccard: {r['jaccard']:.3f}")
        print(f"    H3(T) gen on through-v paths: {r['h3_through_v']}")
        print(f"    H3(T) gen on non-through-v: {r['h3_not_through_v']}")
        print(f"    Emb gen on through-v (should=0): {r['emb_through_v']}")
        print(f"    Total 3-paths: {r['n_paths_3']}, through-v: {r['n_through_v']}")

    # Success summary
    print(f"\n{'='*60}")
    print("SUCCESSES — generator direction (summary)")
    print(f"{'='*60}")
    if successes:
        jaccards = [r['jaccard'] for _, _, r in successes]
        h3_supports = [r['h3_gen_support_size'] for _, _, r in successes]
        emb_supports = [r['emb_gen_support_size'] for _, _, r in successes]
        overlaps = [r['support_overlap'] for _, _, r in successes]
        h3_thru = [r['h3_through_v'] for _, _, r in successes]

        print(f"  Jaccard similarity: min={min(jaccards):.3f}, max={max(jaccards):.3f}, mean={np.mean(jaccards):.3f}")
        print(f"  H3(T) gen support: min={min(h3_supports)}, max={max(h3_supports)}, mean={np.mean(h3_supports):.1f}")
        print(f"  Emb gen support: min={min(emb_supports)}, max={max(emb_supports)}, mean={np.mean(emb_supports):.1f}")
        print(f"  Overlap: min={min(overlaps)}, max={max(overlaps)}, mean={np.mean(overlaps):.1f}")
        print(f"  H3(T) gen through-v: min={min(h3_thru)}, max={max(h3_thru)}, mean={np.mean(h3_thru):.1f}")

    # Comparison
    print(f"\n{'='*60}")
    print("COMPARISON: Failures vs Successes")
    print(f"{'='*60}")
    if failures and successes:
        fail_j = [r['jaccard'] for _, _, r in failures]
        succ_j = [r['jaccard'] for _, _, r in successes]
        fail_h3thru = [r['h3_through_v'] for _, _, r in failures]
        succ_h3thru = [r['h3_through_v'] for _, _, r in successes]
        fail_h3supp = [r['h3_gen_support_size'] for _, _, r in failures]
        succ_h3supp = [r['h3_gen_support_size'] for _, _, r in successes]

        print(f"  Jaccard — fail: {[f'{x:.3f}' for x in fail_j]}, succ mean: {np.mean(succ_j):.3f}")
        print(f"  H3(T) through-v — fail: {fail_h3thru}, succ mean: {np.mean(succ_h3thru):.1f}")
        print(f"  H3(T) support — fail: {fail_h3supp}, succ mean: {np.mean(succ_h3supp):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
