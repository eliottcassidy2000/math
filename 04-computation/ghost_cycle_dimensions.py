"""
ghost_cycle_dimensions.py — Dimensional analysis of the Ghost Cycle Theorem

At n=7, through-v-only cycles are ALWAYS boundaries.
At n=8, they USUALLY are but occasionally fail (~0.25%).

Key question: Is this a dimensional accident or a structural constraint?

Hypothesis: dim(tv-only ker(d_3)) is always small enough relative to
dim(im(d_4)) that containment is "forced" at n=7 but can fail at n=8.

We compute:
  - dim(ker(d_3))
  - dim(im(d_4)) = dim(Omega_4) - dim(ker(d_4))
  - codim(im(d_4), ker(d_3)) = beta_3
  - dim(tv-only subspace of ker(d_3))
  - dim(tv-only subspace of im(d_4))

If dim(tv-only ker(d_3)) <= dim(tv-only im(d_4)), that explains containment.

Also: the single-face theorem says each through-v 4-path contributes
exactly one old face. This constrains how im(d_4) relates to the
old/new decomposition.

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


def ghost_cycle_dimensions(A, n, v):
    """Compute all relevant dimensions for the Ghost Cycle Theorem."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths_3 = ap.get(3, [])
    paths_4 = ap.get(4, [])
    paths_2 = ap.get(2, [])
    paths_5 = ap.get(5, [])

    if not paths_3 or not paths_4:
        return None

    through_v_3 = [i for i, p in enumerate(paths_3) if v in p]
    not_through_v_3 = [i for i, p in enumerate(paths_3) if v not in p]
    through_v_4 = [i for i, p in enumerate(paths_4) if v in p]
    not_through_v_4 = [i for i, p in enumerate(paths_4) if v not in p]

    if not through_v_3 or not not_through_v_3:
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
    ob5 = get_omega(ap, 5)

    # d_3: Omega_3 -> A_2
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    d3o = bd3 @ ob3.T % PRIME

    # ker(d_3) in Omega coords
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    ker_d3_A = ker_d3 @ ob3 % PRIME

    if ker_d3_A.shape[0] == 0:
        return None

    # d_4: Omega_4 -> A_3
    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    im_d4_rows = (bd4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4_rows.copy(), PRIME)) if im_d4_rows.shape[0] > 0 else 0

    # d_5: Omega_5 -> A_4 (for ker(d_4))
    idx4 = {p: i for i, p in enumerate(paths_4)}
    bd5 = np.zeros((len(paths_4), len(paths_5)), dtype=np.int64) if paths_5 else np.zeros((len(paths_4), 0), dtype=np.int64)
    for j, path in enumerate(paths_5):
        for sign, face in boundary_faces(path):
            if face in idx4:
                bd5[idx4[face], j] = (bd5[idx4[face], j] + sign) % PRIME

    d4o = bd4 @ ob4.T % PRIME if ob4.shape[0] > 0 else np.zeros((len(paths_3), 0), dtype=np.int64)

    # Through-v-only subspace of ker(d_3)
    ker_old_proj = ker_d3_A[:, not_through_v_3] % PRIME
    rk_old_proj = int(_gauss_rank_np(ker_old_proj.copy(), PRIME))
    dim_tv_only_ker = ker_d3_A.shape[0] - rk_old_proj

    # Through-v-only subspace of im(d_4)
    if im_d4_rows.shape[0] > 0:
        im_old_proj = im_d4_rows[:, not_through_v_3] % PRIME
        rk_im_old_proj = int(_gauss_rank_np(im_old_proj.copy(), PRIME))
        dim_tv_only_im_d4 = rk_d4 - rk_im_old_proj
    else:
        dim_tv_only_im_d4 = 0

    # Check actual containment
    if dim_tv_only_ker > 0:
        ker_old_proj_T = ker_old_proj.T % PRIME
        _, tv_basis = _gauss_nullbasis_modp(
            ker_old_proj_T.astype(np.int32), ker_old_proj_T.shape[0], ker_old_proj_T.shape[1], PRIME)
        tv_coeffs = np.array(tv_basis, dtype=np.int64) if tv_basis else np.zeros((0, ker_d3_A.shape[0]), dtype=np.int64)
        tv_cycles_A = tv_coeffs @ ker_d3_A % PRIME if tv_coeffs.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)

        if tv_cycles_A.shape[0] > 0 and im_d4_rows.shape[0] > 0:
            combined = np.vstack([im_d4_rows, tv_cycles_A]) % PRIME
            rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
            all_in_im_d4 = (rk_combined == rk_d4)
        else:
            all_in_im_d4 = (tv_cycles_A.shape[0] == 0)
    else:
        all_in_im_d4 = True

    # Through-v stats for 4-paths in Omega
    dim_omega4 = ob4.shape[0]
    # How many omega_4 basis vectors are through-v only?
    if ob4.shape[0] > 0:
        ob4_old_proj = ob4[:, not_through_v_4] % PRIME if not_through_v_4 else np.zeros((ob4.shape[0], 0), dtype=np.int64)
        rk_ob4_old = int(_gauss_rank_np(ob4_old_proj.copy(), PRIME)) if ob4_old_proj.shape[1] > 0 else 0
        dim_tv_only_omega4 = dim_omega4 - rk_ob4_old
    else:
        dim_tv_only_omega4 = 0

    return {
        'dim_omega3': ob3.shape[0],
        'dim_omega4': dim_omega4,
        'dim_ker_d3': ker_d3_A.shape[0],
        'dim_im_d4': rk_d4,
        'beta3': ker_d3_A.shape[0] - rk_d4,
        'dim_tv_only_ker': dim_tv_only_ker,
        'dim_tv_only_im_d4': dim_tv_only_im_d4,
        'dim_tv_only_omega4': dim_tv_only_omega4,
        'n_tv_3': len(through_v_3),
        'n_old_3': len(not_through_v_3),
        'n_tv_4': len(through_v_4),
        'n_old_4': len(not_through_v_4),
        'all_in_im_d4': all_in_im_d4,
        'containment_gap': dim_tv_only_ker - dim_tv_only_im_d4,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"GHOST CYCLE DIMENSIONAL ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 500 if n <= 7 else 300

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = ghost_cycle_dimensions(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        has_tv = [r for r in results if r['dim_tv_only_ker'] > 0]
        no_tv = [r for r in results if r['dim_tv_only_ker'] == 0]
        print(f"\n  Has tv-only cycles: {len(has_tv)}/{len(results)}")

        if has_tv:
            # Dimension statistics
            print(f"\n  --- DIMENSIONS (cases with tv-only cycles) ---")
            for key in ['dim_omega3', 'dim_omega4', 'dim_ker_d3', 'dim_im_d4',
                        'dim_tv_only_ker', 'dim_tv_only_im_d4', 'dim_tv_only_omega4',
                        'containment_gap']:
                vals = [r[key] for r in has_tv]
                dist = Counter(vals)
                print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")
                if len(dist) <= 12:
                    print(f"    distribution: {dict(sorted(dist.items()))}")

            # Key test: is dim_tv_only_ker <= dim_tv_only_im_d4?
            contained = sum(1 for r in has_tv if r['dim_tv_only_ker'] <= r['dim_tv_only_im_d4'])
            print(f"\n  dim(tv-only ker) <= dim(tv-only im_d4): {contained}/{len(has_tv)}")

            # Check containment
            all_bdy = sum(1 for r in has_tv if r['all_in_im_d4'])
            print(f"  Actually contained in im(d_4): {all_bdy}/{len(has_tv)}")

            if all_bdy < len(has_tv):
                print(f"\n  COUNTEREXAMPLES:")
                for r in has_tv:
                    if not r['all_in_im_d4']:
                        print(f"    dim_tv_ker={r['dim_tv_only_ker']}, dim_tv_im={r['dim_tv_only_im_d4']}, "
                              f"gap={r['containment_gap']}, "
                              f"dim_ker={r['dim_ker_d3']}, dim_im={r['dim_im_d4']}")

            # Codim analysis: what fraction of ker(d_3) is im(d_4)?
            codim_vals = [r['dim_ker_d3'] - r['dim_im_d4'] for r in has_tv]
            print(f"\n  beta_3 (codim): {Counter(codim_vals)}")

            # Ratio analysis
            ratios = [r['dim_tv_only_ker'] / max(r['dim_ker_d3'], 1) for r in has_tv]
            print(f"\n  dim(tv-only ker)/dim(ker): min={min(ratios):.3f}, max={max(ratios):.3f}, mean={np.mean(ratios):.3f}")

            im_ratios = [r['dim_im_d4'] / max(r['dim_ker_d3'], 1) for r in has_tv]
            print(f"  dim(im_d4)/dim(ker): min={min(im_ratios):.3f}, max={max(im_ratios):.3f}, mean={np.mean(im_ratios):.3f}")

        if no_tv:
            print(f"\n  --- DIMENSIONS (cases WITHOUT tv-only cycles) ---")
            print(f"  dim_omega3: mean={np.mean([r['dim_omega3'] for r in no_tv]):.1f}")
            print(f"  dim_ker_d3: mean={np.mean([r['dim_ker_d3'] for r in no_tv]):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
