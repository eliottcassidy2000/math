"""
tv_cycle_structure.py — Structure of through-v-only cycles

A tv-only cycle z ∈ ker(d_3^T) satisfies:
  1. z ∈ Ω_3(T) (Omega constraints)
  2. d_3^{tv→tv}(z) = 0  (tv boundary vanishes)
  3. d_3^{tv→old}(z) = 0  (v-deletion faces cancel)

Similarly, a tv-only boundary b ∈ im(d_4^T) satisfies:
  b = d_4^{tv→tv}(w_tv) for some w_tv with d_4^{tv→old}(w_tv) ∈ im(d_4^{old→old})

The Ghost Cycle Theorem (K_tv = B_tv) means condition 2+3 forces z to be
in the image of d_4^{tv→tv} projected appropriately.

KEY IDEA: The "v-deletion face map" φ: tv 3-chains → old 2-chains
defined by φ(z) = [d_3^{tv→old}(z)] sends tv 3-paths to their v-deletion faces.

A tv-only cycle z satisfies φ(z) = 0. So K_tv ⊂ ker(φ) ∩ ker(d_3^{tv→tv}).

Similarly, there's a "v-deletion face map" ψ: tv 4-chains → old 3-chains
defined by ψ(w) = [d_4^{tv→old}(w)].

For a tv-only boundary b = d_4^{tv→tv}(w_tv), we need ψ(w_tv) ∈ im(d_4^{old→old}).

Study: dim(ker(φ) ∩ ker(d_3^{tv→tv})) vs dim(K_tv)
These should be equal if all Omega constraints are already accounted for.

Author: opus-2026-03-10-S59
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


def tv_cycle_block_analysis(A, n, v):
    """Analyze the block structure of boundary maps for tv-only cycles."""
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
    tv2 = [i for i, p in enumerate(paths_2) if v in p]
    old2 = [i for i, p in enumerate(paths_2) if v not in p]

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

    # Build d_3 and d_4 in full A-path coords
    idx2 = {p: i for i, p in enumerate(paths_2)}
    idx3 = {p: i for i, p in enumerate(paths_3)}

    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int64)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % PRIME

    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int64)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % PRIME

    # Extract block maps on A-path level (NOT Omega level)
    # d_3^{tv→old}: bd3[old2, tv3]
    d3_tv_old = bd3[np.ix_(old2, tv3)] % PRIME if old2 and tv3 else np.zeros((0, 0), dtype=np.int64)
    # d_3^{tv→tv}: bd3[tv2, tv3]
    d3_tv_tv = bd3[np.ix_(tv2, tv3)] % PRIME if tv2 and tv3 else np.zeros((0, 0), dtype=np.int64)
    # d_4^{tv→tv}: bd4[tv3, tv4]
    d4_tv_tv = bd4[np.ix_(tv3, tv4)] % PRIME if tv3 and tv4 else np.zeros((0, 0), dtype=np.int64)
    # d_4^{tv→old}: bd4[old3, tv4]
    d4_tv_old = bd4[np.ix_(old3, tv4)] % PRIME if old3 and tv4 else np.zeros((0, 0), dtype=np.int64)
    # d_4^{old→old}: bd4[old3, old4]
    d4_old_old = bd4[np.ix_(old3, old4)] % PRIME if old3 and old4 else np.zeros((0, 0), dtype=np.int64)

    # Omega-restricted versions: need to work in Omega_3^tv and Omega_4^tv
    # Omega_3^tv = vectors in Omega_3 supported on tv3 indices
    # Omega_4^tv = vectors in Omega_4 supported on tv4 indices

    # Project Omega basis to tv/old
    # ob3 rows are Omega_3 basis vectors in A_3 coords
    ob3_tv = ob3[:, tv3] % PRIME if tv3 else np.zeros((ob3.shape[0], 0), dtype=np.int64)
    ob3_old = ob3[:, old3] % PRIME if old3 else np.zeros((ob3.shape[0], 0), dtype=np.int64)

    # Omega_3 vectors supported only on tv3: ob3_old must be zero
    rk_ob3_old = int(_gauss_rank_np(ob3_old.copy(), PRIME)) if ob3_old.shape[1] > 0 else 0
    # Find basis for ker(ob3_old) = Omega vectors with zero old component
    if rk_ob3_old > 0:
        ob3_old_T = ob3_old.T % PRIME
        _, tv_omega_basis = _gauss_nullbasis_modp(
            ob3_old_T.astype(np.int32), ob3_old_T.shape[0], ob3_old_T.shape[1], PRIME)
        tv_omega = np.array(tv_omega_basis, dtype=np.int64) if tv_omega_basis else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        tv_omega = np.eye(ob3.shape[0], dtype=np.int64)

    dim_omega3_tv_only = tv_omega.shape[0]

    # Map tv-only Omega_3 vectors to A_3 tv-coords
    if dim_omega3_tv_only > 0:
        tv_omega_A_tv = (tv_omega @ ob3_tv) % PRIME
    else:
        tv_omega_A_tv = np.zeros((0, len(tv3)), dtype=np.int64)

    # Apply d_3^{tv→old} and d_3^{tv→tv} to tv-only Omega vectors
    if dim_omega3_tv_only > 0 and d3_tv_old.shape[0] > 0:
        phi_img = d3_tv_old @ tv_omega_A_tv.T % PRIME  # φ(z) = d_3^{tv→old}(z)
        rk_phi = int(_gauss_rank_np(phi_img.copy(), PRIME))
    else:
        rk_phi = 0

    if dim_omega3_tv_only > 0 and d3_tv_tv.shape[0] > 0:
        dtv_img = d3_tv_tv @ tv_omega_A_tv.T % PRIME  # d_3^{tv→tv}(z)
        rk_dtv = int(_gauss_rank_np(dtv_img.copy(), PRIME))
    else:
        rk_dtv = 0

    # Combined: d_3 restricted to tv-only Omega vectors
    if dim_omega3_tv_only > 0:
        bd3_on_tv_omega = bd3 @ (tv_omega @ ob3).T % PRIME  # full d_3 on tv-only Omega vectors
        rk_full_d3_tv = int(_gauss_rank_np(bd3_on_tv_omega.copy(), PRIME))
    else:
        rk_full_d3_tv = 0

    # ker(d_3) restricted to tv-only Omega: these ARE K_tv
    dim_K_tv = dim_omega3_tv_only - rk_full_d3_tv

    # Now do the same for im(d_4) → tv-only
    # B_tv = im(d_4) vectors that are tv-only
    # im(d_4) in A_3 coords
    im_d4 = (bd4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths_3)), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    # B_tv: im_d4 vectors with zero old3 coords
    if im_d4.shape[0] > 0:
        im_d4_old = im_d4[:, old3] % PRIME
        rk_im_d4_old = int(_gauss_rank_np(im_d4_old.copy(), PRIME))
        dim_B_tv = rk_d4 - rk_im_d4_old
    else:
        dim_B_tv = 0

    # Full ker(d_3)
    d3o = bd3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)
    dim_K = ker_d3.shape[0]
    beta3 = dim_K - rk_d4

    # Exact sequence check
    # d_3^2 = 0 check: d_3 ∘ d_4 = 0 (always, by construction)
    # ker(d_3^{tv→old} ∘ Ω_3^{tv}) ∩ ker(d_3^{tv→tv} ∘ Ω_3^{tv}) should give K_tv
    # But ker of d_3 FULL is what we need

    # KEY: compare rk_phi (rank of v-deletion face map on tv-only Omega)
    # with dim_omega3_tv_only
    # dim(ker(φ)) = dim_omega3_tv_only - rk_phi
    dim_ker_phi = dim_omega3_tv_only - rk_phi

    return {
        'n_tv3': len(tv3), 'n_old3': len(old3),
        'n_tv4': len(tv4), 'n_old4': len(old4),
        'n_tv2': len(tv2), 'n_old2': len(old2),
        'dim_omega3': ob3.shape[0],
        'dim_omega4': ob4.shape[0],
        'dim_omega3_tv_only': dim_omega3_tv_only,
        'dim_K': dim_K,
        'rk_d4': rk_d4,
        'beta3': beta3,
        'dim_K_tv': dim_K_tv,
        'dim_B_tv': dim_B_tv,
        'ghost_cycle': dim_K_tv == dim_B_tv,
        'rk_phi': rk_phi,  # rank of v-deletion face map on tv-only Omega_3
        'rk_d3_tv_tv': rk_dtv,  # rank of d_3^{tv→tv} on tv-only Omega_3
        'rk_full_d3_tv': rk_full_d3_tv,
        'dim_ker_phi': dim_ker_phi,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"TV-ONLY CYCLE BLOCK STRUCTURE AT n={n}")
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
                r = tv_cycle_block_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

        elapsed = time.time() - t0
        print(f"  {len(results)} (T,v) pairs, {elapsed:.1f}s")

        has_tv = [r for r in results if r['dim_omega3_tv_only'] > 0]
        print(f"\n  Has tv-only Omega_3: {len(has_tv)}/{len(results)}")

        if has_tv:
            for key in ['dim_omega3_tv_only', 'dim_K_tv', 'dim_B_tv',
                        'rk_phi', 'rk_d3_tv_tv', 'rk_full_d3_tv', 'dim_ker_phi']:
                vals = [r[key] for r in has_tv]
                dist = Counter(vals)
                print(f"  {key}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")
                if len(dist) <= 15:
                    print(f"    distribution: {dict(sorted(dist.items()))}")

            # Ghost Cycle check
            gc = sum(1 for r in has_tv if r['ghost_cycle'])
            print(f"\n  Ghost Cycle (K_tv = B_tv): {gc}/{len(has_tv)}")

            # KEY: relationship between rk_phi and rk_d3_tv_tv
            # The full d_3 rank on tv-only = rk_phi + rk_d3_tv_tv? No...
            # They might overlap. Let's check.
            print(f"\n  rk_phi + rk_d3_tv_tv vs rk_full_d3_tv:")
            sum_vs_full = [(r['rk_phi'] + r['rk_d3_tv_tv'], r['rk_full_d3_tv']) for r in has_tv]
            diffs = [s - f for s, f in sum_vs_full]
            print(f"    (rk_phi + rk_dtv) - rk_full: {dict(sorted(Counter(diffs).items()))}")
            # If diff = 0, the two maps have disjoint images
            # If diff > 0, they overlap

            # Dimension of tv-only cycles from block structure:
            # K_tv = dim(tv-only Omega_3) - rk(full d_3 on tv-only)
            # This should match dim_K_tv computed above
            k_tv_check = [r['dim_omega3_tv_only'] - r['rk_full_d3_tv'] for r in has_tv]
            k_tv_actual = [r['dim_K_tv'] for r in has_tv]
            mismatch = sum(1 for a, b in zip(k_tv_check, k_tv_actual) if a != b)
            print(f"\n  K_tv via block = K_tv direct: {len(has_tv) - mismatch}/{len(has_tv)}")

            # What fraction of tv-only Omega is killed by φ (v-deletion) vs d_3^{tv→tv}?
            phi_frac = [r['rk_phi'] / max(r['dim_omega3_tv_only'], 1) for r in has_tv]
            dtv_frac = [r['rk_d3_tv_tv'] / max(r['dim_omega3_tv_only'], 1) for r in has_tv]
            full_frac = [r['rk_full_d3_tv'] / max(r['dim_omega3_tv_only'], 1) for r in has_tv]
            print(f"\n  Fraction of tv-only Omega killed by:")
            print(f"    φ (v-del faces): min={min(phi_frac):.3f}, max={max(phi_frac):.3f}, mean={np.mean(phi_frac):.3f}")
            print(f"    d_3^{'{tv→tv}'}: min={min(dtv_frac):.3f}, max={max(dtv_frac):.3f}, mean={np.mean(dtv_frac):.3f}")
            print(f"    full d_3:     min={min(full_frac):.3f}, max={max(full_frac):.3f}, mean={np.mean(full_frac):.3f}")

        no_tv = [r for r in results if r['dim_omega3_tv_only'] == 0]
        if no_tv:
            print(f"\n  Cases with no tv-only Omega_3: {len(no_tv)}")
            # These trivially have K_tv = 0 = B_tv
            print(f"    dim_K_tv=0: {sum(1 for r in no_tv if r['dim_K_tv'] == 0)}/{len(no_tv)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
