"""
tv_cycles_are_boundaries.py — Are through-v-only cycles always boundaries?

If every cycle in ker(d_3) supported entirely on through-v paths is in im(d_4),
then the H_3 generator can never be through-v-only, proving HYP-408.

This is equivalent to: the through-v-only part of ker(d_3) ∩ im(d_4) has the
same dimension as the through-v-only part of ker(d_3).

Test: compute ker(π|_{ker(d_3)}) and check if it's contained in im(d_4).

Author: opus-2026-03-09-S58
"""
import sys
import time
import numpy as np
from collections import Counter
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


def tv_cycles_boundary_test(A, n, v):
    """Check if all through-v-only cycles are boundaries."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])
    paths_2 = ap.get(2, [])

    if not paths_3:
        return None

    through_v = [i for i, p in enumerate(paths_3) if v in p]
    not_through_v = [i for i, p in enumerate(paths_3) if v not in p]

    if not through_v or not not_through_v:
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

    # ker(d_3) in A_3 coords
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    ker_d3_A = ker_d3 @ ob3 % PRIME

    if ker_d3_A.shape[0] == 0:
        return None

    # im(d_4) in A_3 coords
    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    im_d4 = (bd4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    # Find through-v-only cycles: ker(d_3) vectors that are zero on non-through-v coords
    ker_old_proj = ker_d3_A[:, not_through_v] % PRIME
    rk_old_proj = int(_gauss_rank_np(ker_old_proj.copy(), PRIME))
    n_tv_only = ker_d3_A.shape[0] - rk_old_proj  # dim of through-v-only cycle space

    if n_tv_only == 0:
        return {
            'n_tv_only': 0,
            'all_boundaries': True,
            'dim_ker_d3': ker_d3_A.shape[0],
            'rk_d4': rk_d4,
        }

    # Find a basis for the through-v-only subspace of ker(d_3)
    # We want vectors c in R^{dim_ker} such that c @ ker_d3_A has zero on not_through_v coords
    # i.e., c @ ker_old_proj = 0
    # This is the LEFT null space of ker_old_proj, = null space of ker_old_proj^T
    ker_old_proj_T = ker_old_proj.T % PRIME
    _, tv_basis = _gauss_nullbasis_modp(
        ker_old_proj_T.astype(np.int32), ker_old_proj_T.shape[0], ker_old_proj_T.shape[1], PRIME)
    tv_coeffs = np.array(tv_basis, dtype=np.int64) if tv_basis else np.zeros((0, ker_d3_A.shape[0]), dtype=np.int64)

    # Map back to A_3 coords
    tv_cycles_A = tv_coeffs @ ker_d3_A % PRIME if tv_coeffs.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)

    # Check: are all tv_cycles in im(d_4)?
    if tv_cycles_A.shape[0] > 0 and im_d4.shape[0] > 0:
        combined = np.vstack([im_d4, tv_cycles_A]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        all_in_im_d4 = (rk_combined == rk_d4)
        extra = rk_combined - rk_d4
    elif tv_cycles_A.shape[0] > 0:
        all_in_im_d4 = False
        extra = tv_cycles_A.shape[0]
    else:
        all_in_im_d4 = True
        extra = 0

    return {
        'n_tv_only': n_tv_only,
        'all_boundaries': all_in_im_d4,
        'extra_over_im_d4': extra,
        'dim_ker_d3': ker_d3_A.shape[0],
        'rk_d4': rk_d4,
    }


def main():
    for n in [6, 7, 8]:
        print(f"\n{'='*70}")
        print(f"THROUGH-v-ONLY CYCLES = BOUNDARIES? AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        target = 1000 if n <= 7 else 500
        t0 = time.time()

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) < 1:
                continue

            for v_cand in range(n):
                r = tv_cycles_boundary_test(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} pairs, {elapsed:.1f}s")

        has_tv = [r for r in results if r['n_tv_only'] > 0]
        print(f"  Has through-v-only cycles: {len(has_tv)}/{len(results)}")

        if has_tv:
            all_bdy = sum(1 for r in has_tv if r['all_boundaries'])
            print(f"  Of those, ALL are boundaries: {all_bdy}/{len(has_tv)}")

            if all_bdy < len(has_tv):
                print(f"  COUNTEREXAMPLE FOUND!")
                for r in has_tv:
                    if not r['all_boundaries']:
                        print(f"    n_tv_only={r['n_tv_only']}, extra={r['extra_over_im_d4']}, "
                              f"dim_ker={r['dim_ker_d3']}, rk_d4={r['rk_d4']}")
                        break
            else:
                print(f"  ✓ ALL through-v-only cycles are boundaries!")

        # Distribution of n_tv_only
        tvonly_dist = Counter(r['n_tv_only'] for r in results)
        print(f"\n  n_tv_only distribution: {dict(sorted(tvonly_dist.items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
