"""
tv_subcomplex_acyclicity.py — Test whether the through-v projected complex is acyclic

TWO APPROACHES:

1. MODIFIED TV COMPLEX: Take tv Omega chains, project boundary to tv part.
   Is d^{tv}_{p-1} o d^{tv}_p = 0? (chain complex?)
   If yes, does H_3^{tv} = 0?

2. RELATIVE COMPLEX: Q_p = Omega_p(T) / Omega_p(T\v).
   The proper relative homology H_3(T, T\v).
   Does K_tv map to zero in H_3(T, T\v)?

3. RANK(i_*) ANALYSIS: Does rank(i_*: H_3(T\v) -> H_3(T)) determine codim?
   LES gives: dim H_3(T, T\v) = beta_3(T) - rank(i_*) + beta_2(T\v)
   Since beta_2(T\v) = 0: dim H_3(T,T\v) = beta_3 - rank(i_*)

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


def tv_complex_analysis(A, n, v):
    """Analyze the through-v subcomplex structure."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1

    cc = full_chain_complex_modp(A, n, n - 1)
    if cc['bettis'].get(3, 0) != 1:
        return None

    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, n1 - 1)

    ap = enumerate_all_allowed(A, n, min(n-1, 7))

    paths = {deg: ap.get(deg, []) for deg in range(8)}
    old_idx = {deg: [i for i, p in enumerate(paths[deg]) if v not in p] for deg in range(8)}
    tv_idx = {deg: [i for i, p in enumerate(paths[deg]) if v in p] for deg in range(8)}

    def get_omega(deg):
        ps = paths[deg]
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

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

    # Get Omega bases and boundary maps for degrees 2-5
    ob = {}
    for deg in [2, 3, 4, 5]:
        ob[deg] = get_omega(deg)

    D = {}
    for deg in [3, 4, 5]:
        D[deg] = build_bd(deg, deg - 1)

    # ===== APPROACH 1: Modified TV complex =====
    # Omega_p in A_p coordinates: ob[p] rows give basis vectors

    # For each degree, extract tv-projected boundary
    # d_p^{mod}: C_p^{tv} -> C_{p-1}^{tv}
    # In Omega coords, the boundary of Omega basis vector b is D_p @ b
    # The tv projection keeps only tv coordinates

    def tv_projected_boundary(deg):
        """Get the tv-projected boundary matrix in Omega coords.
        Maps: Omega_deg (full) -> tv part of A_{deg-1}"""
        full_bd = D[deg] @ ob[deg].T % PRIME  # A_{deg-1} x dim(Omega_deg)
        tv_rows = sorted(tv_idx[deg - 1])
        if not tv_rows:
            return np.zeros((0, ob[deg].shape[0]), dtype=np.int64)
        return full_bd[tv_rows, :] % PRIME

    # Check: d^{tv}_{p-1} o d^{tv}_p = 0 ?
    # This is: (tv proj of D_{p-1}) @ Omega_{p-1} basis @ ... hmm, this isn't right
    # The modified complex acts on tv-ONLY Omega basis elements

    # Let me work in A_p coordinates directly.
    # The tv subspace of Omega_p: restrict Omega_p basis to tv coordinates only
    # omega_p^{tv} = ob[p][:, tv_idx[p]] (extract tv columns)
    # But this is NOT a basis for the tv subspace of Omega_p!
    # An Omega basis vector might have both old and tv components.

    # The tv-ONLY elements of Omega_p: those with zero old coords
    # These form a subspace V_p^{tv-only} of Omega_p.
    def tv_only_omega(deg):
        """Find basis for {z in Omega_deg : z has zero old coords}."""
        ob_deg = ob[deg]
        if ob_deg.shape[0] == 0:
            return np.zeros((0, len(paths[deg])), dtype=np.int64)
        old_cols = sorted(old_idx[deg])
        if not old_cols:
            return ob_deg.copy()  # all paths are tv
        # Find null space of ob_deg[:, old_cols]
        M = ob_deg[:, old_cols].T % PRIME  # old_cols x dim(Omega_deg)
        _, nb = _gauss_nullbasis_modp(M.astype(np.int32), M.shape[0], M.shape[1], PRIME)
        if not nb:
            return np.zeros((0, len(paths[deg])), dtype=np.int64)
        coeffs = np.array(nb, dtype=np.int64)
        return coeffs @ ob_deg % PRIME  # in A_deg coords

    V_tv = {}
    for deg in [2, 3, 4, 5]:
        V_tv[deg] = tv_only_omega(deg)

    # Modified boundary: D_p restricted to tv-only Omega_p -> tv-only part of A_{p-1}
    # The boundary of a tv-only element z has BOTH old and tv parts in general.
    # The tv part: D[p][tv_rows, :] @ z^T
    # For the MODIFIED complex, we project to tv part only.

    # Check d^{mod}_{p-1} o d^{mod}_p
    # d^{mod}_p(z) = tv_proj(D_p @ z)
    # d^{mod}_{p-1}(d^{mod}_p(z)) = tv_proj(D_{p-1} @ tv_proj(D_p @ z))
    # This involves projecting to tv at each step.

    # Let's check at p=4: d^{mod}_3(d^{mod}_4(z)) for z in V_tv[4]
    chain_complex_ok = True
    dd_failures = 0
    dd_total = 0

    for p in [4, 5]:
        if p > 5 or V_tv.get(p) is None or V_tv[p].shape[0] == 0:
            continue
        if D.get(p) is None or D.get(p-1) is None:
            continue

        tv_rows_pm1 = sorted(tv_idx[p - 1])
        tv_rows_pm2 = sorted(tv_idx[p - 2])

        if not tv_rows_pm1 or not tv_rows_pm2:
            continue

        for k in range(V_tv[p].shape[0]):
            z = V_tv[p][k]
            # d^{mod}_p(z)
            dz_full = D[p] @ z % PRIME
            dz_tv = dz_full[tv_rows_pm1] % PRIME  # tv projection

            # Now apply d^{mod}_{p-1}: need to find the Omega_{p-1} element
            # that matches dz_tv on tv coords and has 0 on old coords.
            # But dz_tv is in A_{p-1} coords (tv subset), NOT in Omega_{p-1}.
            # For d^{mod} to be well-defined, dz_tv must be expressible as
            # an element of Omega_{p-1} restricted to tv coords.

            # Actually, for the modified complex approach, we work with
            # A_p^{tv} (raw tv paths) directly, not Omega.
            # d^{mod}_p: A_p^{tv} -> A_{p-1}^{tv} via D_p[tv_rows, tv_cols]

            # Let's just check D_{p-1}[tv,tv] @ D_p[tv,tv] = 0
            pass

        D_p_tt = D[p][np.ix_(tv_rows_pm1, sorted(tv_idx[p]))] if tv_rows_pm1 and tv_idx[p] else np.zeros((0,0), dtype=np.int64)
        D_pm1_tt = D[p-1][np.ix_(tv_rows_pm2, tv_rows_pm1)] if tv_rows_pm2 and tv_rows_pm1 else np.zeros((0,0), dtype=np.int64)

        if D_pm1_tt.shape[1] > 0 and D_p_tt.shape[1] > 0 and D_pm1_tt.shape[1] == D_p_tt.shape[0]:
            dd = D_pm1_tt @ D_p_tt % PRIME
            dd_nnz = int(np.count_nonzero(dd))
            dd_total += 1
            if dd_nnz > 0:
                dd_failures += 1
                chain_complex_ok = False

    # ===== APPROACH 2: rank(i_*) analysis =====
    # Compute rank of i_*: H_3(T\v) -> H_3(T)

    # Embed H_3(T\v) generator into Omega_3(T)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(n1-1, 5))
    ob3_Tv = get_omega(3)  # This is for T, not T\v

    # Get T\v Omega_3 and chain complex
    paths_3_Tv = ap_Tv.get(3, [])
    if paths_3_Tv:
        P_Tv, nr_Tv, nc_Tv = _build_constraint_matrix(ap_Tv, 3, PRIME)
        if P_Tv is not None:
            _, nb_Tv = _gauss_nullbasis_modp(P_Tv, nr_Tv, nc_Tv, PRIME)
            ob3_Tv_actual = np.array(nb_Tv, dtype=np.int64) if nb_Tv else np.zeros((0, nc_Tv), dtype=np.int64)
        else:
            ob3_Tv_actual = np.eye(len(paths_3_Tv), dtype=np.int64)
    else:
        ob3_Tv_actual = np.zeros((0, 0), dtype=np.int64)

    b3_Tv = cc_Tv['bettis'].get(3, 0)

    # If b3_Tv = 0, rank(i_*) = 0 trivially
    if b3_Tv == 0:
        rank_istar = 0
    else:
        # Get H_3(T\v) generator in A_3(T\v) coords
        D3_Tv = build_bd(3, 2)  # This uses T's paths, not T\v's
        # Need to build D3 for T\v separately
        paths_2_Tv = ap_Tv.get(2, [])
        idx2_Tv = {p: i for i, p in enumerate(paths_2_Tv)}
        D3_Tv_mat = np.zeros((len(paths_2_Tv), len(paths_3_Tv)), dtype=np.int64)
        for j, path in enumerate(paths_3_Tv):
            for sign, face in boundary_faces(path):
                if face in idx2_Tv:
                    D3_Tv_mat[idx2_Tv[face], j] = (D3_Tv_mat[idx2_Tv[face], j] + sign) % PRIME

        d3o_Tv = D3_Tv_mat @ ob3_Tv_actual.T % PRIME if ob3_Tv_actual.shape[0] > 0 else np.zeros((1,1), dtype=np.int64)
        if d3o_Tv.size > 0:
            _, kv_Tv = _gauss_nullbasis_modp(d3o_Tv.astype(np.int32), d3o_Tv.shape[0], d3o_Tv.shape[1], PRIME)
            ker_d3_Tv = np.array(kv_Tv, dtype=np.int64) @ ob3_Tv_actual % PRIME if kv_Tv else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
        else:
            ker_d3_Tv = ob3_Tv_actual.copy()

        # Get im(d_4^{T\v})
        paths_4_Tv = ap_Tv.get(4, [])
        if paths_4_Tv:
            P4_Tv, nr4, nc4 = _build_constraint_matrix(ap_Tv, 4, PRIME)
            if P4_Tv is not None:
                _, nb4 = _gauss_nullbasis_modp(P4_Tv, nr4, nc4, PRIME)
                ob4_Tv = np.array(nb4, dtype=np.int64) if nb4 else np.zeros((0, nc4), dtype=np.int64)
            else:
                ob4_Tv = np.eye(len(paths_4_Tv), dtype=np.int64)
        else:
            ob4_Tv = np.zeros((0, 0), dtype=np.int64)

        idx3_Tv = {p: i for i, p in enumerate(paths_3_Tv)}
        D4_Tv_mat = np.zeros((len(paths_3_Tv), len(paths_4_Tv)), dtype=np.int64) if paths_4_Tv else np.zeros((len(paths_3_Tv), 0), dtype=np.int64)
        for j, path in enumerate(paths_4_Tv):
            for sign, face in boundary_faces(path):
                if face in idx3_Tv:
                    D4_Tv_mat[idx3_Tv[face], j] = (D4_Tv_mat[idx3_Tv[face], j] + sign) % PRIME

        im_d4_Tv = (D4_Tv_mat @ ob4_Tv.T % PRIME).T if ob4_Tv.shape[0] > 0 else np.zeros((0, len(paths_3_Tv)), dtype=np.int64)
        rk_d4_Tv = int(_gauss_rank_np(im_d4_Tv.copy(), PRIME)) if im_d4_Tv.shape[0] > 0 else 0

        # Embed into T's A_3: map T\v paths to T paths
        embed_3 = np.zeros((len(paths_3_Tv), len(paths[3])), dtype=np.int64)
        idx3_T = {p: i for i, p in enumerate(paths[3])}
        for j, path_v in enumerate(paths_3_Tv):
            mapped = tuple(remaining[x] for x in path_v)
            if mapped in idx3_T:
                embed_3[j, idx3_T[mapped]] = 1

        # H_3(T\v) generators in A_3(T) coords
        if ker_d3_Tv.shape[0] > 0:
            H3_Tv_in_T = ker_d3_Tv @ embed_3 % PRIME
        else:
            H3_Tv_in_T = np.zeros((0, len(paths[3])), dtype=np.int64)

        # Check rank of i_*: rank([im_d4_T | H3_Tv_in_T]) - rank(im_d4_T)
        im_d4_T = (D[4] @ ob[4].T % PRIME).T if ob[4].shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
        rk_d4_T = int(_gauss_rank_np(im_d4_T.copy(), PRIME)) if im_d4_T.shape[0] > 0 else 0

        if H3_Tv_in_T.shape[0] > 0 and im_d4_T.shape[0] > 0:
            combined = np.vstack([im_d4_T, H3_Tv_in_T]) % PRIME
        elif H3_Tv_in_T.shape[0] > 0:
            combined = H3_Tv_in_T
        else:
            combined = im_d4_T
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME)) if combined.shape[0] > 0 else 0
        rank_istar = rk_combined - rk_d4_T

    # ===== Collect results =====
    dim_K = K_dim = 0
    # Compute K_tv dimension
    full_bd = D[3] @ ob[3].T % PRIME
    if full_bd.size > 0:
        _, kv = _gauss_nullbasis_modp(full_bd.astype(np.int32), full_bd.shape[0], full_bd.shape[1], PRIME)
        ker_d3_T = np.array(kv, dtype=np.int64) @ ob[3] % PRIME if kv else np.zeros((0, len(paths[3])), dtype=np.int64)
    else:
        ker_d3_T = ob[3].copy()

    K = ker_d3_T
    old3 = sorted(old_idx[3])
    pi_K = K[:, old3] % PRIME if old3 else np.zeros((K.shape[0], 0), dtype=np.int64)
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME)) if pi_K.size > 0 else 0
    dim_Ktv = K.shape[0] - rk_pi_K

    return {
        'chain_complex_ok': chain_complex_ok,
        'dd_failures': dd_failures,
        'dd_total': dd_total,
        'dim_V_tv': {deg: V_tv[deg].shape[0] for deg in [2, 3, 4, 5]},
        'rank_istar': rank_istar,
        'b3_Tv': b3_Tv,
        'dim_Ktv': dim_Ktv,
        'dim_K': K.shape[0],
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"TV SUBCOMPLEX ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 200 if n == 7 else 80

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            for v_cand in range(min(n, 3)):  # Test first 3 vertices to save time
                r = tv_complex_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
                if len(results) % 50 == 0:
                    print(f"  Progress: {len(results)} pairs, {time.time()-t0:.1f}s")

        elapsed = time.time() - t0
        print(f"\n  Total: {len(results)} (T,v) pairs, {elapsed:.1f}s")

        if not results:
            continue

        # APPROACH 1: Chain complex check
        cc_ok = sum(1 for r in results if r['chain_complex_ok'])
        dd_fail = sum(r['dd_failures'] for r in results)
        print(f"\n  MODIFIED TV COMPLEX (d_tv o d_tv = 0?):")
        print(f"    Chain complex condition: {cc_ok}/{len(results)} OK")
        print(f"    Total dd failures: {dd_fail}")

        # TV-only Omega dimensions
        print(f"\n  TV-ONLY OMEGA DIMENSIONS:")
        for deg in [2, 3, 4, 5]:
            vals = [r['dim_V_tv'][deg] for r in results]
            if vals:
                print(f"    V_tv[{deg}]: mean={np.mean(vals):.1f}, "
                      f"min={min(vals)}, max={max(vals)}")

        # APPROACH 2: rank(i_*) analysis
        print(f"\n  RANK(i_*) ANALYSIS:")
        istar_dist = Counter(r['rank_istar'] for r in results)
        print(f"    rank(i_*) distribution: {dict(sorted(istar_dist.items()))}")
        print(f"    b3(T\\v) distribution: {dict(Counter(r['b3_Tv'] for r in results))}")

        # K_tv dimension
        ktv_dist = Counter(r['dim_Ktv'] for r in results)
        print(f"    K_tv dimension distribution: {dict(sorted(ktv_dist.items()))}")

        # KEY: rank(i_*) vs K_tv
        # When rank(i_*) = 0 and K_tv > 0, we need a non-LES argument
        istar0_ktv_pos = sum(1 for r in results if r['rank_istar'] == 0 and r['dim_Ktv'] > 0)
        istar1_ktv_pos = sum(1 for r in results if r['rank_istar'] == 1 and r['dim_Ktv'] > 0)
        print(f"\n  CRITICAL CASES:")
        print(f"    rank(i_*)=0 AND K_tv>0: {istar0_ktv_pos}")
        print(f"    rank(i_*)=1 AND K_tv>0: {istar1_ktv_pos}")
        print(f"    rank(i_*)=0 AND K_tv=0: {sum(1 for r in results if r['rank_istar'] == 0 and r['dim_Ktv'] == 0)}")
        print(f"    rank(i_*)=1 AND K_tv=0: {sum(1 for r in results if r['rank_istar'] == 1 and r['dim_Ktv'] == 0)}")

        # Cross-tabulation
        print(f"\n  CROSS-TABULATION (rank(i_*) x b3(T\\v)):")
        for istar in sorted(istar_dist.keys()):
            for b3tv in sorted(set(r['b3_Tv'] for r in results)):
                count = sum(1 for r in results if r['rank_istar'] == istar and r['b3_Tv'] == b3tv)
                if count > 0:
                    print(f"    i_*={istar}, b3_Tv={b3tv}: {count}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
