"""
tv_bounding_chain.py — Extract the bounding 4-chain for through-v-only 3-cycles

For each tournament T with b3=1 at n=7 that has K_tv > 0 (nontrivial
through-v-only cycles), find the explicit 4-chain w with d_4(w) = z
for each basis element z of K_tv.

GOAL: Understand the STRUCTURE of the bounding chain. Specifically:
1. How much of w is old (4-paths not through v) vs tv (4-paths through v)?
2. Is w always mixed (both old and tv), or can it be pure tv?
3. Is there a pattern in how the old part of w cancels the old-spillover
   from the tv part?

Also tests the KEY THEOREM from LES analysis:
  H_3(T, T\v) = coker(i_*) since beta_2(T\v) = 0
  rank(i_*) relates to codim via:
    rank(i_*) = 0 and K_tv = 0 => codim = 1
    rank(i_*) = 1 => codim = 1 (LES)

Author: kind-pasteur-2026-03-10-S50
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


def analyze_bounding_chain(A, n, v):
    """Find and analyze bounding 4-chains for through-v-only 3-cycles."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1

    cc = full_chain_complex_modp(A, n, n - 1)
    if cc['bettis'].get(3, 0) != 1:
        return None

    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, n1 - 1)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    ap = enumerate_all_allowed(A, n, min(n-1, 6))

    paths = {deg: ap.get(deg, []) for deg in range(6)}

    old_idx = {deg: [i for i, p in enumerate(paths[deg]) if v not in p] for deg in range(6)}
    tv_idx = {deg: [i for i, p in enumerate(paths[deg]) if v in p] for deg in range(6)}

    def get_omega(deg):
        ps = paths[deg]
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    ob3 = get_omega(3)
    ob4 = get_omega(4)

    def build_bd(deg_dom, deg_cod):
        ps_dom = paths[deg_dom]
        ps_cod = paths[deg_cod]
        if not ps_dom or not ps_cod:
            return np.zeros((len(ps_cod), len(ps_dom)), dtype=np.int64)
        idx_cod = {p: i for i, p in enumerate(ps_cod)}
        D = np.zeros((len(ps_cod), len(ps_dom)), dtype=np.int64)
        for j, path in enumerate(ps_dom):
            for sign, face in boundary_faces(path):
                if face in idx_cod:
                    D[idx_cod[face], j] = (D[idx_cod[face], j] + sign) % PRIME
        return D

    D3 = build_bd(3, 2)
    D4 = build_bd(4, 3)

    # ker(d_3) in A_3 coords
    d3o = D3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    K = ker_d3 @ ob3 % PRIME

    # im(d_4) in A_3 coords
    im_d4 = (D4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    dim_K = K.shape[0]
    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    old3 = sorted(old_idx[3])
    tv3 = sorted(tv_idx[3])

    # Find K_tv basis (through-v-only cycles)
    pi_K = K[:, old3] % PRIME
    pi_K_T = pi_K.T % PRIME
    if pi_K_T.size > 0:
        _, tv_basis_coeffs = _gauss_nullbasis_modp(
            pi_K_T.astype(np.int32), pi_K_T.shape[0], pi_K_T.shape[1], PRIME)
        K_tv_coeffs = np.array(tv_basis_coeffs, dtype=np.int64) if tv_basis_coeffs else np.zeros((0, K.shape[0]), dtype=np.int64)
    else:
        K_tv_coeffs = np.zeros((0, K.shape[0]), dtype=np.int64)

    K_tv = K_tv_coeffs @ K % PRIME if K_tv_coeffs.shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    dim_Ktv = K_tv.shape[0]

    if dim_Ktv == 0:
        return {'dim_Ktv': 0, 'b3_Tv': b3_Tv}

    # For each K_tv basis element, find the bounding 4-chain
    # z = d_4(w) where w is in Omega_4
    # im_d4 rows are the boundary vectors. We need to express z as a combination.

    # Build the full boundary image matrix: each column of D4 @ ob4.T gives boundary in A_3
    bd_matrix = D4 @ ob4.T % PRIME  # shape: |A_3| x dim(Omega_4)

    results = []
    for i in range(dim_Ktv):
        z = K_tv[i] % PRIME  # through-v-only 3-cycle

        # Solve: bd_matrix @ w = z (mod prime)
        # bd_matrix: |A_3| x dim(Omega_4), z: |A_3|
        # Form [bd_matrix | -z] and find null space
        dim_omega4 = ob4.shape[0]
        neg_z = ((-z) % PRIME).reshape(-1, 1)
        ext_mat = np.hstack([bd_matrix, neg_z]) % PRIME
        _, nb = _gauss_nullbasis_modp(ext_mat.astype(np.int32), ext_mat.shape[0], ext_mat.shape[1], PRIME)

        if not nb:
            results.append({'in_image': False, 'dim_Ktv': dim_Ktv, 'idx': i})
            continue

        # Find a null vector with last coordinate = 1
        null_vecs = np.array(nb, dtype=np.int64)
        # We want null_vecs[:, -1] != 0 (the z coefficient)
        last_col = null_vecs[:, -1] % PRIME
        good = np.where(last_col != 0)[0]
        if len(good) == 0:
            results.append({'in_image': False, 'dim_Ktv': dim_Ktv, 'idx': i})
            continue

        # Take first good vector and normalize
        nv = null_vecs[good[0]].copy()
        inv_last = pow(int(nv[-1]), PRIME - 2, PRIME)
        w = nv[:-1] * inv_last % PRIME  # Omega_4 coefficients

        # Convert w to A_4 coordinates
        w_A4 = w @ ob4 % PRIME  # in A_4 raw path coordinates

        # Verify: D4 @ w_A4 should equal z
        check = D4 @ w_A4 % PRIME
        if not np.array_equal(check % PRIME, z % PRIME):
            # Verification failed
            results.append({'in_image': False, 'verify_fail': True, 'dim_Ktv': dim_Ktv, 'idx': i})
            continue

        # Analyze w_A4: old vs tv components
        old4 = sorted(old_idx[4])
        tv4 = sorted(tv_idx[4])

        w_old_count = sum(1 for j in old4 if w_A4[j] % PRIME != 0)
        w_tv_count = sum(1 for j in tv4 if w_A4[j] % PRIME != 0)
        w_old_total = len(old4)
        w_tv_total = len(tv4)

        # Check if w can be made pure tv (is there a solution with w_old = 0?)
        # This requires: D4[old3, tv4] @ w_tv = z[old3] = 0 (since z is tv-only)
        # AND D4[tv3, tv4] @ w_tv = z[tv3]
        # Since z is tv-only, z[old3] = 0 is automatic.
        # So we need: D4[old3, :] @ w_A4 = 0 on old coords
        # which is: D4_oo @ w_old + D4_ot @ w_tv = 0

        # The old-spillover from tv part: D4_ot @ w_tv
        D4_oo_block = D4[np.ix_(old3, old4)] if old3 and old4 else np.zeros((0,0), dtype=np.int64)
        D4_ot_block = D4[np.ix_(old3, tv4)] if old3 and tv4 else np.zeros((0,0), dtype=np.int64)

        w_tv_raw = w_A4[tv4] % PRIME if tv4 else np.zeros(0, dtype=np.int64)
        w_old_raw = w_A4[old4] % PRIME if old4 else np.zeros(0, dtype=np.int64)

        if D4_ot_block.size > 0 and len(w_tv_raw) > 0:
            spillover = D4_ot_block @ w_tv_raw % PRIME
            spillover_nonzero = int(np.count_nonzero(spillover))
        else:
            spillover_nonzero = 0

        results.append({
            'in_image': True,
            'dim_Ktv': dim_Ktv,
            'idx': i,
            'w_old_nonzero': w_old_count,
            'w_tv_nonzero': w_tv_count,
            'w_old_total': w_old_total,
            'w_tv_total': w_tv_total,
            'spillover_nonzero': spillover_nonzero,
            'b3_Tv': b3_Tv,
        })

    return {
        'dim_Ktv': dim_Ktv,
        'b3_Tv': b3_Tv,
        'chains': results,
    }


def main():
    n = 7
    print(f"{'='*70}")
    print(f"BOUNDING CHAIN ANALYSIS AT n={n}")
    print(f"{'='*70}")

    rng = np.random.RandomState(42)
    all_results = []
    total_chains = 0
    in_image_count = 0
    old_needed = 0
    pure_tv_possible = 0
    t0 = time.time()
    target = 100  # tournaments with K_tv > 0

    for trial in range(100000):
        if len(all_results) >= target:
            break
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        for v_cand in range(n):
            r = analyze_bounding_chain(A, n, v_cand)
            if r is None or r['dim_Ktv'] == 0:
                continue
            all_results.append(r)
            for c in r['chains']:
                total_chains += 1
                if c.get('in_image', False):
                    in_image_count += 1
                    if c['w_old_nonzero'] > 0:
                        old_needed += 1
                    else:
                        pure_tv_possible += 1
                    if c.get('spillover_nonzero', 0) > 0:
                        pass  # spillover exists
            if len(all_results) % 20 == 0:
                print(f"  Progress: {len(all_results)} pairs with K_tv>0, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    print(f"\n  Total: {len(all_results)} (T,v) with K_tv>0, {total_chains} basis vectors, {elapsed:.1f}s")

    # Summary statistics
    print(f"\n  IN IMAGE (z in im(d_4)): {in_image_count}/{total_chains}")
    print(f"  OLD COMPONENTS NEEDED: {old_needed}/{in_image_count}")
    print(f"  PURE TV SOLUTION EXISTS: {pure_tv_possible}/{in_image_count}")

    # Detailed breakdown
    print(f"\n  BOUNDING CHAIN STRUCTURE:")
    old_counts = [c['w_old_nonzero'] for r in all_results for c in r['chains'] if c.get('in_image')]
    tv_counts = [c['w_tv_nonzero'] for r in all_results for c in r['chains'] if c.get('in_image')]
    spillovers = [c['spillover_nonzero'] for r in all_results for c in r['chains'] if c.get('in_image')]

    if old_counts:
        print(f"    w_old nonzero entries: mean={np.mean(old_counts):.1f}, "
              f"min={min(old_counts)}, max={max(old_counts)}")
        print(f"    w_tv nonzero entries: mean={np.mean(tv_counts):.1f}, "
              f"min={min(tv_counts)}, max={max(tv_counts)}")
        print(f"    spillover nonzero: mean={np.mean(spillovers):.1f}, "
              f"min={min(spillovers)}, max={max(spillovers)}")

    # K_tv dimension distribution
    ktv_dims = [r['dim_Ktv'] for r in all_results]
    print(f"\n  K_tv dimension dist: {dict(Counter(ktv_dims))}")

    # b3(T\v) distribution
    b3tv_vals = [r['b3_Tv'] for r in all_results]
    print(f"  b3(T\\v) dist: {dict(Counter(b3tv_vals))}")

    # Key question: is w ALWAYS mixed (old + tv)?
    all_mixed = all(c.get('w_old_nonzero', 0) > 0 and c.get('w_tv_nonzero', 0) > 0
                     for r in all_results for c in r['chains'] if c.get('in_image'))
    print(f"\n  ALL bounding chains are mixed (old+tv): {all_mixed}")

    # Any verify failures?
    verify_fails = sum(1 for r in all_results for c in r['chains'] if c.get('verify_fail'))
    not_in_image = sum(1 for r in all_results for c in r['chains'] if not c.get('in_image'))
    print(f"  Verify failures: {verify_fails}")
    print(f"  NOT in image: {not_in_image}")

    # Ratio analysis: how much of w is old vs tv?
    if old_counts and tv_counts:
        ratios = [o / max(o + t, 1) for o, t in zip(old_counts, tv_counts)]
        print(f"\n  OLD FRACTION in bounding chain: mean={np.mean(ratios):.3f}, "
              f"min={min(ratios):.3f}, max={max(ratios):.3f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
