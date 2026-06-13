"""
codim1_universality.py — Test the codim=1 universality more thoroughly

KEY FINDING FROM codim_analysis.py:
  When we project ker(d_3^T) and im(d_4^T) onto the "old" coordinate
  subspace (paths not through v), the codimension is ALWAYS 1.

  This is equivalent to saying:
    The old-projection of H_3(T) is 1-dimensional.

  Why? Because:
    ker(d_3^T)_old_proj / im(d_4^T)_old_proj ≅ (something 1-dim)

  And rank(i_*)=1 iff the embedded H_3(T\v) gen maps to this 1-dim quotient
  nontrivially.

NEW QUESTION: Is the old-projection of ker(d_3^T) = ker(d_3^{T\v})?
That would explain everything — ker(d_3^{T\v}) projects to same thing as
ker(d_3^T)_old, and the H_3(T\v) gen already spans the codim-1 complement
in ker(d_3^{T\v})/im(d_4^{T\v}), so its projection should span the
codim-1 complement in the old-projected version too.

But wait — at n=8 failures, this breaks. So the question is more subtle.

This script tests:
1. dim(ker(d_3^T)_old_proj) vs dim(ker(d_3^{T\v}))
2. Whether ker(d_3^{T\v}) embeds into ker(d_3^T)_old_proj
3. The codim of im(d_4^{T\v}) in ker(d_3^{T\v}) vs codim of im(d_4^T)_old in ker(d_3^T)_old

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


def codim1_deep_analysis(A, n, v):
    """Deep analysis of old-projection universality."""
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

    # Build boundary matrices
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}
    idx3Tv = {p: i for i, p in enumerate(paths_3_Tv)}

    bd3_T = np.zeros((len(paths_2_T), len(paths_3_T)), dtype=np.int64)
    for j, path in enumerate(paths_3_T):
        for sign, face in boundary_faces(path):
            if face in idx2T:
                bd3_T[idx2T[face], j] = (bd3_T[idx2T[face], j] + sign) % PRIME

    bd4_T = np.zeros((len(paths_3_T), len(paths_4_T)), dtype=np.int64)
    for j, path in enumerate(paths_4_T):
        for sign, face in boundary_faces(path):
            if face in idx3T:
                bd4_T[idx3T[face], j] = (bd4_T[idx3T[face], j] + sign) % PRIME

    bd3_Tv = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_3_Tv):
        for sign, face in boundary_faces(path):
            if face in idx2Tv:
                bd3_Tv[idx2Tv[face], j] = (bd3_Tv[idx2Tv[face], j] + sign) % PRIME

    bd4_Tv = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64)
    for j, path in enumerate(paths_4_Tv):
        for sign, face in boundary_faces(path):
            if face in idx3Tv:
                bd4_Tv[idx3Tv[face], j] = (bd4_Tv[idx3Tv[face], j] + sign) % PRIME

    # ker(d_3^T) and im(d_4^T) in A_3(T) coords
    d3o_T = bd3_T @ ob3_T.T % PRIME
    _, kv_T = _gauss_nullbasis_modp(d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME) if d3o_T.size > 0 else (0, [])
    ker_d3T_A = (np.array(kv_T, dtype=np.int64) @ ob3_T % PRIME) if kv_T else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    im_d4T_A = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4T = int(_gauss_rank_np(im_d4T_A.copy(), PRIME)) if im_d4T_A.shape[0] > 0 else 0

    # ker(d_3^{T\v}) and im(d_4^{T\v}) in A_3(T\v) coords
    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME) if d3o_Tv.size > 0 else (0, [])
    ker_d3Tv_A = (np.array(kv_Tv, dtype=np.int64) @ ob3_Tv % PRIME) if kv_Tv else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    im_d4Tv_A = (bd4_Tv @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
    rk_d4Tv = int(_gauss_rank_np(im_d4Tv_A.copy(), PRIME)) if im_d4Tv_A.shape[0] > 0 else 0

    # Old 3-path indices
    old_3_idx = [i for i, p in enumerate(paths_3_T) if v not in p]

    # Embedding map: A_3(T\v) -> A_3(T)
    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    # Embed ker(d_3^{T\v}) into A_3(T)
    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # === TEST 1: dim(ker(d_3^T)_old_proj) vs dim(ker(d_3^{T\v})) ===
    ker_T_old_proj = ker_d3T_A[:, old_3_idx] % PRIME if ker_d3T_A.shape[0] > 0 and old_3_idx else np.zeros((0, 0), dtype=np.int64)
    rk_ker_T_old = int(_gauss_rank_np(ker_T_old_proj.copy(), PRIME)) if ker_T_old_proj.size > 0 else 0

    dim_ker_Tv = ker_d3Tv_A.shape[0]  # = dim ker(d_3^{T\v})

    # === TEST 2: Does embedded ker(d_3^{T\v}) span the same old-projection? ===
    emb_old_proj = embedded_ker[:, old_3_idx] % PRIME if embedded_ker.shape[0] > 0 and old_3_idx else np.zeros((0, 0), dtype=np.int64)
    rk_emb_old = int(_gauss_rank_np(emb_old_proj.copy(), PRIME)) if emb_old_proj.size > 0 else 0

    # Does span(ker_T_old_proj) = span(emb_old_proj)?
    if rk_ker_T_old > 0 and rk_emb_old > 0:
        test_span = np.vstack([ker_T_old_proj, emb_old_proj]) % PRIME
        rk_union = int(_gauss_rank_np(test_span.copy(), PRIME))
        same_span = (rk_union == rk_ker_T_old == rk_emb_old)
    else:
        rk_union = max(rk_ker_T_old, rk_emb_old)
        same_span = False

    # === TEST 3: im(d_4^T)_old vs im(d_4^{T\v}) ===
    im_T_old_proj = im_d4T_A[:, old_3_idx] % PRIME if im_d4T_A.shape[0] > 0 and old_3_idx else np.zeros((0, 0), dtype=np.int64)
    rk_im_T_old = int(_gauss_rank_np(im_T_old_proj.copy(), PRIME)) if im_T_old_proj.size > 0 else 0

    # Embed im(d_4^{T\v}) into A_3(T)
    emb_im_d4Tv = im_d4Tv_A @ embed % PRIME if im_d4Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    emb_im_d4Tv_old = emb_im_d4Tv[:, old_3_idx] % PRIME if emb_im_d4Tv.shape[0] > 0 and old_3_idx else np.zeros((0, 0), dtype=np.int64)
    rk_emb_im_Tv = int(_gauss_rank_np(emb_im_d4Tv_old.copy(), PRIME)) if emb_im_d4Tv_old.size > 0 else 0

    # Does im(d_4^{T\v}) ⊂ im(d_4^T)_old? (in old coords)
    if rk_im_T_old > 0 and rk_emb_im_Tv > 0:
        test_im = np.vstack([im_T_old_proj, emb_im_d4Tv_old]) % PRIME
        rk_im_union = int(_gauss_rank_np(test_im.copy(), PRIME))
        im_Tv_in_im_T = (rk_im_union == rk_im_T_old)
    else:
        rk_im_union = max(rk_im_T_old, rk_emb_im_Tv)
        im_Tv_in_im_T = (rk_emb_im_Tv == 0)

    # rank(i_*)
    if im_d4T_A.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4T_A, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4T_A
    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4T

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'rk_ker_T_old': rk_ker_T_old,
        'dim_ker_Tv': dim_ker_Tv,
        'rk_emb_old': rk_emb_old,
        'same_span': same_span,
        'rk_union': rk_union,
        'rk_im_T_old': rk_im_T_old,
        'rk_emb_im_Tv': rk_emb_im_Tv,
        'im_Tv_in_im_T': im_Tv_in_im_T,
        'rk_im_union': rk_im_union,
        'rk_d4T': rk_d4T,
        'rk_d4Tv': rk_d4Tv,
        'b3_Tv': cc_Tv['bettis'].get(3, 0),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"CODIM-1 UNIVERSALITY DEEP ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        checked = 0
        target = 200 if n == 7 else 120
        t0 = time.time()

        for trial in range(30000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            checked += 1
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue
            for v_cand in range(n):
                r = codim1_deep_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
            if len(results) % 50 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)} BAD vertices, {elapsed:.1f}s")

        elapsed = time.time() - t0
        fails = [r for r in results if r['rank_istar'] == 0]
        succs = [r for r in results if r['rank_istar'] == 1]
        print(f"  {len(results)} BAD vertices ({len(fails)} fail, {len(succs)} succ), {elapsed:.1f}s")

        # TEST 1: dim comparison
        print(f"\n  TEST 1: rk(ker(d3_T)_old) vs dim(ker(d3_Tv))")
        diffs = [r['rk_ker_T_old'] - r['dim_ker_Tv'] for r in results]
        diff_dist = Counter(diffs)
        print(f"    difference distribution: {dict(sorted(diff_dist.items()))}")

        # TEST 2: same span?
        same_count = sum(1 for r in results if r['same_span'])
        print(f"\n  TEST 2: span(ker(d3_T)_old) = span(emb ker(d3_Tv))?")
        print(f"    YES: {same_count}/{len(results)} ({100*same_count/max(len(results),1):.1f}%)")
        if not all(r['same_span'] for r in results):
            for r in results:
                if not r['same_span']:
                    print(f"      MISMATCH: rk_ker_T_old={r['rk_ker_T_old']}, "
                          f"rk_emb_old={r['rk_emb_old']}, rk_union={r['rk_union']}, "
                          f"rank_istar={r['rank_istar']}")
                    break

        # TEST 3: im(d_4^{T\v}) ⊂ im(d_4^T)_old?
        im_sub_count = sum(1 for r in results if r['im_Tv_in_im_T'])
        print(f"\n  TEST 3: im(d4_Tv) ⊂ im(d4_T)_old?")
        print(f"    YES: {im_sub_count}/{len(results)} ({100*im_sub_count/max(len(results),1):.1f}%)")
        if im_sub_count < len(results):
            for r in results:
                if not r['im_Tv_in_im_T']:
                    print(f"      FAIL: rk_im_T_old={r['rk_im_T_old']}, "
                          f"rk_emb_im_Tv={r['rk_emb_im_Tv']}, rk_union={r['rk_im_union']}, "
                          f"rank_istar={r['rank_istar']}")
                    break

        # Rank comparisons
        print(f"\n  Rank comparisons:")
        print(f"    rk(d4_T) mean={np.mean([r['rk_d4T'] for r in results]):.1f}")
        print(f"    rk(d4_Tv) mean={np.mean([r['rk_d4Tv'] for r in results]):.1f}")
        print(f"    rk(im_T_old) mean={np.mean([r['rk_im_T_old'] for r in results]):.1f}")
        print(f"    rk(emb_im_Tv) mean={np.mean([r['rk_emb_im_Tv'] for r in results]):.1f}")

        if fails:
            print(f"\n  FAILURE details:")
            for r in fails[:3]:
                print(f"    out={r['out_deg']}: rk_ker_T_old={r['rk_ker_T_old']}, "
                      f"dim_ker_Tv={r['dim_ker_Tv']}, same_span={r['same_span']}, "
                      f"im_Tv_in_im_T={r['im_Tv_in_im_T']}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
