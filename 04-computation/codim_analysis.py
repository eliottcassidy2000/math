"""
codim_analysis.py — Codimension analysis for i_*-injectivity

The key algebraic picture:
  ker(d_3^T) has dim = rk(d_4) + 1  (since b3=1)
  im(d_4) has codimension 1 in ker(d_3)
  The embedded H_3(T\v) gen is a SINGLE vector in A_3(T)
  rank(i_*)=0 iff this vector lies in the codim-1 subspace im(d_4)

Think of it as: im(d_4) is a "hyperplane" in ker(d_3).
A random vector in ker(d_3) misses it with probability (p-1)/p ~ 1.
But the embedded gen is NOT random — it's highly structured.

This script computes the "effective codimension" — how close is the
embedded gen to being in im(d_4)?

Specifically: project the embedded gen onto ker(d_3)/im(d_4) ≅ F_p.
The coefficient (0 or nonzero) determines i_*.

Also: look at the RESTRICTED picture — among the "old" part of ker(d_3),
what is the codimension of the old part of im(d_4)?

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


def restricted_codim_analysis(A, n, v):
    """Analyze the restricted (old-path) picture."""
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

    # Build boundary matrices
    idx2T = {p: i for i, p in enumerate(paths_2_T)}
    idx3T = {p: i for i, p in enumerate(paths_3_T)}
    idx2Tv = {p: i for i, p in enumerate(paths_2_Tv)}

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

    # Work in A_3(T) coordinates throughout
    # ker(d_3^T) via Omega
    d3o_T = bd3_T @ ob3_T.T % PRIME
    _, ker_vecs_T = _gauss_nullbasis_modp(
        d3o_T.astype(np.int32), d3o_T.shape[0], d3o_T.shape[1], PRIME) if d3o_T.size > 0 else (0, [])
    ker_d3T_A = (np.array(ker_vecs_T, dtype=np.int64) @ ob3_T % PRIME) if ker_vecs_T else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # im(d_4^T)
    im_d4_A = (bd4_T @ ob4_T.T % PRIME).T if ob4_T.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4_A.copy(), PRIME)) if im_d4_A.shape[0] > 0 else 0

    # Embedded ker(d_3^{T\v})
    d3o_Tv = bd3_Tv @ ob3_Tv.T % PRIME
    _, ker_vecs_Tv = _gauss_nullbasis_modp(
        d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME) if d3o_Tv.size > 0 else (0, [])
    ker_d3Tv_A = (np.array(ker_vecs_Tv, dtype=np.int64) @ ob3_Tv % PRIME) if ker_vecs_Tv else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)

    embed = np.zeros((len(paths_3_Tv), len(paths_3_T)), dtype=np.int64)
    for jv, path_v in enumerate(paths_3_Tv):
        mapped = tuple(remaining[x] for x in path_v)
        if mapped in idx3T:
            embed[jv, idx3T[mapped]] = 1

    embedded_ker = ker_d3Tv_A @ embed % PRIME if ker_d3Tv_A.shape[0] > 0 else np.zeros((0, len(paths_3_T)), dtype=np.int64)

    # rank(i_*)
    if im_d4_A.shape[0] > 0 and embedded_ker.shape[0] > 0:
        combined = np.vstack([im_d4_A, embedded_ker]) % PRIME
    elif embedded_ker.shape[0] > 0:
        combined = embedded_ker
    else:
        combined = im_d4_A
    rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
    rank_istar = rk_combined - rk_d4

    # === RESTRICTED ANALYSIS ===
    # Consider ONLY the "old" coordinate subspace (paths not through v)
    old_3_idx = [i for i, p in enumerate(paths_3_T) if v not in p]
    new_3_idx = [i for i, p in enumerate(paths_3_T) if v in p]

    # Restrict ker(d_3^T) to old coords
    # Each ker vector: project onto old coords only
    if ker_d3T_A.shape[0] > 0 and old_3_idx:
        ker_old_proj = ker_d3T_A[:, old_3_idx] % PRIME
        rk_ker_old_proj = int(_gauss_rank_np(ker_old_proj.copy(), PRIME))
    else:
        rk_ker_old_proj = 0

    # Restrict im(d_4^T) to old coords
    if im_d4_A.shape[0] > 0 and old_3_idx:
        im_old_proj = im_d4_A[:, old_3_idx] % PRIME
        rk_im_old_proj = int(_gauss_rank_np(im_old_proj.copy(), PRIME))
    else:
        rk_im_old_proj = 0

    # Restrict embedded ker to old coords (embedded gen is ONLY on old paths already)
    if embedded_ker.shape[0] > 0 and old_3_idx:
        emb_old_proj = embedded_ker[:, old_3_idx] % PRIME
        rk_emb_old = int(_gauss_rank_np(emb_old_proj.copy(), PRIME))
    else:
        rk_emb_old = 0

    # Is emb_old in span of im_old_proj?
    if im_old_proj.shape[0] > 0 and emb_old_proj.shape[0] > 0:
        test = np.vstack([im_d4_A[:, old_3_idx] % PRIME, emb_old_proj]) % PRIME
        rk_test = int(_gauss_rank_np(test.copy(), PRIME))
        emb_old_in_im_old = (rk_test == rk_im_old_proj)
    else:
        emb_old_in_im_old = False

    # Codim of im_old_proj in ker_old_proj
    # Need to compute: rank of [ker_old_proj rows] intersected with span of [im_old_proj rows]
    # Actually: project ker and im to same old-coord space, compute codim
    codim_old = rk_ker_old_proj - rk_im_old_proj

    # Also: what's the restriction of the H_3(T) generator to old coords?
    # The H_3(T) gen is the vector in ker(d3) NOT in im(d4)
    # Its old-coord projection tells us about its "old footprint"

    out_deg = int(sum(A[v]))

    return {
        'rank_istar': rank_istar,
        'out_deg': out_deg,
        'rk_d4': rk_d4,
        'dim_ker_d3T': ker_d3T_A.shape[0],
        'n_old_3': len(old_3_idx),
        'n_new_3': len(new_3_idx),
        'rk_ker_old_proj': rk_ker_old_proj,
        'rk_im_old_proj': rk_im_old_proj,
        'codim_old': codim_old,
        'emb_old_in_im_old': emb_old_in_im_old,
        'rk_emb_old': rk_emb_old,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"CODIMENSION ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        checked = 0
        target = 300 if n == 7 else 150
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
                r = restricted_codim_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

            if len(results) % 100 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)} BAD vertices, {elapsed:.1f}s")

        elapsed = time.time() - t0
        fails = [r for r in results if r['rank_istar'] == 0]
        succs = [r for r in results if r['rank_istar'] == 1]
        print(f"  {len(results)} BAD vertices, {len(fails)} failures, {elapsed:.1f}s")

        # OLD-COORDINATE CODIMENSION
        print(f"\n  --- Old-coordinate codimension ---")
        codims = [r['codim_old'] for r in results]
        codim_dist = Counter(codims)
        print(f"  codim(im_old_proj, ker_old_proj): {dict(sorted(codim_dist.items()))}")

        print(f"\n  rk_ker_old_proj: mean={np.mean([r['rk_ker_old_proj'] for r in results]):.1f}")
        print(f"  rk_im_old_proj:  mean={np.mean([r['rk_im_old_proj'] for r in results]):.1f}")

        # Is embedded gen's old projection in im_old_proj?
        emb_in = sum(1 for r in results if r['emb_old_in_im_old'])
        print(f"\n  emb_old IN im_old_proj: {emb_in}/{len(results)} ({100*emb_in/max(len(results),1):.1f}%)")

        if fails:
            print(f"\n  --- FAILURES ---")
            for r in fails[:5]:
                print(f"    out_deg={r['out_deg']}: codim_old={r['codim_old']}, "
                      f"rk_ker_old={r['rk_ker_old_proj']}, rk_im_old={r['rk_im_old_proj']}, "
                      f"emb_old_in_im={r['emb_old_in_im_old']}")

        if succs:
            print(f"\n  --- SUCCESSES (sample) ---")
            for r in succs[:5]:
                print(f"    out_deg={r['out_deg']}: codim_old={r['codim_old']}, "
                      f"rk_ker_old={r['rk_ker_old_proj']}, rk_im_old={r['rk_im_old_proj']}, "
                      f"emb_old_in_im={r['emb_old_in_im_old']}")

        # CODIM by success/fail
        if fails and succs:
            fail_codim = [r['codim_old'] for r in fails]
            succ_codim = [r['codim_old'] for r in succs]
            print(f"\n  codim_old — fail: {fail_codim[:10]}, succ mean: {np.mean(succ_codim):.2f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
