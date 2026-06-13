"""
ghost_cycle_structure.py — WHY are through-v-only cycles boundaries?

Key observation from ghost_cycle_dimensions.py:
  dim(tv-only ker(d_3)) = dim(tv-only im(d_4))  ALWAYS (in sample)

This is equivalent to: the projection of im(d_4) to tv-only coordinates
has the same rank as the projection of ker(d_3) to tv-only coordinates.

In other words: im(d_4) fills the tv-only subspace of ker(d_3) exactly.

Why? Consider the chain complex restricted to through-v paths:
  C_5^tv -> C_4^tv -d4-> C_3^tv -d3-> C_2^tv

The through-v boundary d_4(sigma) of a through-v 4-chain has both
tv and old components. The tv-only part is the projection pi_tv(d_4(sigma)).

Claim: pi_tv(im(d_4)) >= pi_tv(ker(d_3)) holds because:
1. The old part of d_4 is "generic" — it doesn't constrain the tv part
2. Every tv-only ker(d_3) element is the tv part of some boundary

To investigate: decompose d_4 into tv and old components:
  d_4 = d_4^{tv->tv} + d_4^{tv->old} + d_4^{old->tv} + d_4^{old->old}

For a through-v 4-chain w, d_4(w) has:
  - tv-to-tv: faces of w that are still through-v
  - tv-to-old: the v-deletion face (exactly one per single-face theorem)

For an old 4-chain w, d_4(w):
  - old-to-old: all faces (none use v)
  - old-to-tv: impossible (can't introduce v by deleting a vertex)

So d_4^{old->tv} = 0! This means the tv part of im(d_4) comes
ENTIRELY from through-v 4-chains.

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def decompose_boundary(A, n, v):
    """Decompose d_4 into tv/old block structure and analyze."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])
    paths_2 = ap.get(2, [])

    if not paths_3 or not paths_4:
        return None

    tv3 = [i for i, p in enumerate(paths_3) if v in p]
    old3 = [i for i, p in enumerate(paths_3) if v not in p]
    tv4 = [i for i, p in enumerate(paths_4) if v in p]
    old4 = [i for i, p in enumerate(paths_4) if v not in p]

    if not tv3 or not old3:
        return None

    def get_omega(ap, deg):
        paths = ap.get(deg, [])
        if not paths:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(paths), dtype=np.int64)

    ob3 = get_omega(ap, 3)
    ob4 = get_omega(ap, 4)

    # Build d_4 in A-path coords
    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    # d_4 in Omega coords: d4_omega = bd4 @ ob4^T
    d4_omega = bd4 @ ob4.T % PRIME if ob4.shape[0] > 0 else np.zeros((len(paths_3), 0), dtype=np.int64)

    # Extract blocks: d4_omega[tv3_rows, :] = tv part of image
    # and d4_omega[old3_rows, :] = old part of image
    tv3_set = set(tv3)
    old3_set = set(old3)

    # Verify: old->tv block is zero?
    # An old 4-path (no v) has faces all without v, so bd4[tv3_idx, old4_idx] should be 0
    old4_to_tv3_nonzero = 0
    for j_idx in old4:
        for i_idx in tv3:
            if bd4[i_idx, j_idx] % PRIME != 0:
                old4_to_tv3_nonzero += 1

    # Similarly check tv4->old3 (these are the v-deletion faces)
    tv4_to_old3_nonzero = 0
    for j_idx in tv4:
        for i_idx in old3:
            if bd4[i_idx, j_idx] % PRIME != 0:
                tv4_to_old3_nonzero += 1

    # Build d_3 in Omega coords
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME

    # ker(d_3) in A coords
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    ker_d3_A = ker_d3 @ ob3 % PRIME

    # im(d_4) in A coords
    im_d4 = d4_omega.T % PRIME if d4_omega.shape[1] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)

    # Project to tv coordinates
    tv3_sorted = sorted(tv3)
    old3_sorted = sorted(old3)

    # tv-only subspace of ker(d_3)
    ker_old_proj = ker_d3_A[:, old3_sorted] % PRIME
    rk_ker_old = int(_gauss_rank_np(ker_old_proj.copy(), PRIME))
    dim_tv_only_ker = ker_d3_A.shape[0] - rk_ker_old

    # tv projection of im(d_4)
    im_tv_proj = im_d4[:, tv3_sorted] % PRIME if im_d4.shape[0] > 0 else np.zeros((0, len(tv3_sorted)), dtype=np.int64)
    rk_im_tv = int(_gauss_rank_np(im_tv_proj.copy(), PRIME)) if im_tv_proj.shape[0] > 0 else 0

    # tv projection of ker(d_3)
    ker_tv_proj = ker_d3_A[:, tv3_sorted] % PRIME
    rk_ker_tv = int(_gauss_rank_np(ker_tv_proj.copy(), PRIME))

    # CRITICAL: d_3 restricted to tv-only coords
    # For a tv-only 3-chain z: d_3(z) has both tv and old 2-faces
    # The old 2-faces are the v-deletion faces
    tv2 = [i for i, p in enumerate(paths_2) if v in p]
    old2 = [i for i, p in enumerate(paths_2) if v not in p]

    # d_3 restricted: tv3 -> old2 (the v-deletion face map)
    d3_tv_to_old2 = bd3[np.ix_(old2, tv3)] % PRIME if old2 and tv3 else np.zeros((0, 0), dtype=np.int64)
    rk_d3_tv_old2 = int(_gauss_rank_np(d3_tv_to_old2.T.copy(), PRIME)) if d3_tv_to_old2.size > 0 else 0

    # d_3 restricted: tv3 -> tv2
    d3_tv_to_tv2 = bd3[np.ix_(tv2, tv3)] % PRIME if tv2 and tv3 else np.zeros((0, 0), dtype=np.int64)
    rk_d3_tv_tv2 = int(_gauss_rank_np(d3_tv_to_tv2.T.copy(), PRIME)) if d3_tv_to_tv2.size > 0 else 0

    return {
        'n_tv3': len(tv3), 'n_old3': len(old3),
        'n_tv4': len(tv4), 'n_old4': len(old4),
        'dim_omega3': ob3.shape[0], 'dim_omega4': ob4.shape[0],
        'old4_to_tv3_nonzero': old4_to_tv3_nonzero,
        'tv4_to_old3_nonzero': tv4_to_old3_nonzero,
        'dim_ker_d3': ker_d3_A.shape[0],
        'dim_im_d4': int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0,
        'dim_tv_only_ker': dim_tv_only_ker,
        'rk_im_tv_proj': rk_im_tv,
        'rk_ker_tv_proj': rk_ker_tv,
        'n_tv2': len(tv2), 'n_old2': len(old2),
        'rk_d3_tv_to_old2': rk_d3_tv_old2,
        'rk_d3_tv_to_tv2': rk_d3_tv_tv2,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"BOUNDARY DECOMPOSITION AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 300 if n <= 7 else 200

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = decompose_boundary(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} pairs, {elapsed:.1f}s")

        # Verify old->tv block is zero
        n_bad = sum(1 for r in results if r['old4_to_tv3_nonzero'] > 0)
        print(f"\n  old4->tv3 block nonzero: {n_bad}/{len(results)} (should be 0)")

        # tv4->old3 stats
        tv4_old3_vals = [r['tv4_to_old3_nonzero'] for r in results]
        print(f"  tv4->old3 nonzero entries: min={min(tv4_old3_vals)}, max={max(tv4_old3_vals)}, mean={np.mean(tv4_old3_vals):.1f}")

        # Key: rk(im_d4 tv-proj) vs rk(ker_d3 tv-proj) vs dim(tv-only ker)
        has_tv = [r for r in results if r['dim_tv_only_ker'] > 0]
        if has_tv:
            print(f"\n  --- Cases with tv-only cycles: {len(has_tv)} ---")
            for key in ['dim_tv_only_ker', 'rk_im_tv_proj', 'rk_ker_tv_proj']:
                vals = [r[key] for r in has_tv]
                print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

            # The tv-projection of im(d_4) vs tv-projection of ker(d_3)
            # If rk(im_tv_proj) >= rk(ker_tv_proj), then im(d_4) "covers" all of ker in tv direction
            gap = [r['rk_ker_tv_proj'] - r['rk_im_tv_proj'] for r in has_tv]
            print(f"\n  rk(ker_tv) - rk(im_tv): {dict(sorted(Counter(gap).items()))}")

            # d_3 restricted to tv -> old2
            # This is the map that must vanish for a tv-only cycle
            print(f"\n  d_3 restricted tv3->old2:")
            for key in ['rk_d3_tv_to_old2', 'rk_d3_tv_to_tv2', 'n_tv3', 'n_old2', 'n_tv2']:
                vals = [r[key] for r in has_tv]
                print(f"    {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

            # KEY RATIO: rk(d3_tv_to_old2) / n_tv3
            # This tells us how many "old constraints" there are on tv 3-chains
            old_constraint_ratios = [r['rk_d3_tv_to_old2'] / max(r['n_tv3'], 1) for r in has_tv]
            print(f"\n  rk(d3_tv_old2)/n_tv3: min={min(old_constraint_ratios):.3f}, "
                  f"max={max(old_constraint_ratios):.3f}, mean={np.mean(old_constraint_ratios):.3f}")

            # Dimension of tv 3-chains satisfying d3_tv_old2 = 0
            # This is dim(ker(d3_tv_old2)) = n_tv3 - rk(d3_tv_old2)
            tv_ker_dims = [r['n_tv3'] - r['rk_d3_tv_to_old2'] for r in has_tv]
            print(f"  dim(tv3 with zero old-boundary): min={min(tv_ker_dims)}, "
                  f"max={max(tv_ker_dims)}, mean={np.mean(tv_ker_dims):.1f}")

            # Compare with dim(tv-only ker(d_3))
            print(f"  dim(tv-only ker(d_3)): mean={np.mean([r['dim_tv_only_ker'] for r in has_tv]):.1f}")
            print(f"  Ratio tv_ker_dim / dim_tv_only_ker:")
            ratios = [tv_ker_dims[i] / max(has_tv[i]['dim_tv_only_ker'], 1) for i in range(len(has_tv))]
            print(f"    min={min(ratios):.2f}, max={max(ratios):.2f}, mean={np.mean(ratios):.2f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
