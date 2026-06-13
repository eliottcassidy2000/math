"""
block_decomposition_hyp408.py — Block decomposition approach to proving HYP-408

KEY INSIGHT: Both the constraint matrix P and the boundary matrix D have
block-triangular structure when we decompose paths into old (no v) and
through-v (tv) parts:

  P = [ P_oo  P_ot ]   and   D = [ D_oo  D_ot ]
      [  0    P_tt ]             [  0    D_tt ]

This forces a fiber structure on Omega_p and on ker(d_p).

GOAL: Track all block ranks at n=7 for b3=1 tournaments to find the
algebraic identity that forces codim(pi_old(im d_4), pi_old(ker d_3)) = 1.

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
    enumerate_all_allowed, enumerate_allowed_paths,
    build_adj_lists,
    _build_constraint_matrix, _gauss_rank_np, _gauss_nullbasis_modp,
    full_chain_complex_modp, boundary_faces,
    RANK_PRIME
)

PRIME = RANK_PRIME


def block_analysis(A, n, v, verbose=False):
    """Full block decomposition of the chain complex around degree 3."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1

    # Check Betti numbers
    cc = full_chain_complex_modp(A, n, n - 1)
    if cc['bettis'].get(3, 0) != 1:
        return None

    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, n1 - 1)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    # Enumerate all allowed paths in T
    ap = enumerate_all_allowed(A, n, min(n-1, 6))

    paths = {}
    for deg in range(6):
        paths[deg] = ap.get(deg, [])

    # Classify paths as old/tv
    old_idx = {}
    tv_idx = {}
    for deg in range(6):
        old_idx[deg] = [i for i, p in enumerate(paths[deg]) if v not in p]
        tv_idx[deg] = [i for i, p in enumerate(paths[deg]) if v in p]

    # Build Omega bases for degrees 2,3,4,5
    def get_omega_basis(deg):
        ps = paths[deg]
        if not ps:
            return np.zeros((0, 0), dtype=np.int64)
        P, nr, nc = _build_constraint_matrix(ap, deg, PRIME)
        if P is not None:
            _, nb = _gauss_nullbasis_modp(P, nr, nc, PRIME)
            return np.array(nb, dtype=np.int64) if nb else np.zeros((0, nc), dtype=np.int64)
        return np.eye(len(ps), dtype=np.int64)

    ob = {}
    for deg in [2, 3, 4, 5]:
        ob[deg] = get_omega_basis(deg)

    # Build RAW boundary matrices (before Omega projection)
    def build_raw_boundary(deg_domain, deg_codomain):
        """Build boundary matrix D: A_{deg_codomain} <- A_{deg_domain}"""
        ps_dom = paths[deg_domain]
        ps_cod = paths[deg_codomain]
        if not ps_dom or not ps_cod:
            return np.zeros((len(ps_cod), len(ps_dom)), dtype=np.int64)
        idx_cod = {p: i for i, p in enumerate(ps_cod)}
        D = np.zeros((len(ps_cod), len(ps_dom)), dtype=np.int64)
        for j, path in enumerate(ps_dom):
            for sign, face in boundary_faces(path):
                if face in idx_cod:
                    D[idx_cod[face], j] = (D[idx_cod[face], j] + sign) % PRIME
        return D

    # Boundary in Omega coordinates: d_p = D_p @ Omega_p^T (rows of result = A_{p-1} coords)
    # Then project to Omega_{p-1}: result is in Omega_{p-1} basis if we left-multiply by Omega_{p-1}

    # For our purposes, we work in A_p coordinates (raw path indices)
    # d_p sends Omega_p elements (in A_p) to A_{p-1}, and the result must lie in Omega_{p-1}

    D3 = build_raw_boundary(3, 2)  # A_2 <- A_3
    D4 = build_raw_boundary(4, 3)  # A_3 <- A_4

    # Block decomposition of D3 and D4
    # D3 blocks: rows=A_2 (old/tv), cols=A_3 (old/tv)
    D3_oo = D3[np.ix_(old_idx[2], old_idx[3])] if old_idx[2] and old_idx[3] else np.zeros((0,0), dtype=np.int64)
    D3_ot = D3[np.ix_(old_idx[2], tv_idx[3])] if old_idx[2] and tv_idx[3] else np.zeros((0,0), dtype=np.int64)
    D3_to = D3[np.ix_(tv_idx[2], old_idx[3])] if tv_idx[2] and old_idx[3] else np.zeros((0,0), dtype=np.int64)
    D3_tt = D3[np.ix_(tv_idx[2], tv_idx[3])] if tv_idx[2] and tv_idx[3] else np.zeros((0,0), dtype=np.int64)

    D4_oo = D4[np.ix_(old_idx[3], old_idx[4])] if old_idx[3] and old_idx[4] else np.zeros((0,0), dtype=np.int64)
    D4_ot = D4[np.ix_(old_idx[3], tv_idx[4])] if old_idx[3] and tv_idx[4] else np.zeros((0,0), dtype=np.int64)
    D4_to = D4[np.ix_(tv_idx[3], old_idx[4])] if tv_idx[3] and old_idx[4] else np.zeros((0,0), dtype=np.int64)
    D4_tt = D4[np.ix_(tv_idx[3], tv_idx[4])] if tv_idx[3] and tv_idx[4] else np.zeros((0,0), dtype=np.int64)

    # Verify block triangularity: D_to should be 0
    d3_to_nnz = int(np.count_nonzero(D3_to)) if D3_to.size else 0
    d4_to_nnz = int(np.count_nonzero(D4_to)) if D4_to.size else 0

    # Now compute things in Omega coordinates
    # ker(d_3) in A_3 coords
    d3_omega = D3 @ ob[3].T % PRIME  # A_2 rows, Omega_3 cols
    if d3_omega.size > 0:
        _, kv = _gauss_nullbasis_modp(d3_omega.astype(np.int32), d3_omega.shape[0], d3_omega.shape[1], PRIME)
        ker_d3_omega = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob[3].shape[0]), dtype=np.int64)
    else:
        ker_d3_omega = np.eye(ob[3].shape[0], dtype=np.int64)

    K = ker_d3_omega @ ob[3] % PRIME  # ker(d_3) in A_3 coords, rows are cycle vectors

    # im(d_4) in A_3 coords
    im_d4_omega = (D4 @ ob[4].T % PRIME).T if ob[4].shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    rk_im_d4 = int(_gauss_rank_np(im_d4_omega.copy(), PRIME)) if im_d4_omega.shape[0] > 0 else 0

    dim_K = K.shape[0]
    beta3 = dim_K - rk_im_d4
    if beta3 != 1:
        return None

    # Old/tv projections
    old3 = sorted(old_idx[3])
    tv3 = sorted(tv_idx[3])

    # pi_old(K) and pi_old(im d_4)
    pi_K = K[:, old3] % PRIME if old3 else np.zeros((K.shape[0], 0), dtype=np.int64)
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME)) if pi_K.size > 0 else 0

    pi_B = im_d4_omega[:, old3] % PRIME if im_d4_omega.shape[0] > 0 and old3 else np.zeros((0, len(old3)), dtype=np.int64)
    rk_pi_B = int(_gauss_rank_np(pi_B.copy(), PRIME)) if pi_B.size > 0 else 0

    codim = rk_pi_K - rk_pi_B

    # Also compute: pi_tv(K) and pi_tv(im d_4)
    pi_K_tv = K[:, tv3] % PRIME if tv3 else np.zeros((K.shape[0], 0), dtype=np.int64)
    rk_pi_K_tv = int(_gauss_rank_np(pi_K_tv.copy(), PRIME)) if pi_K_tv.size > 0 else 0

    pi_B_tv = im_d4_omega[:, tv3] % PRIME if im_d4_omega.shape[0] > 0 and tv3 else np.zeros((0, len(tv3)), dtype=np.int64)
    rk_pi_B_tv = int(_gauss_rank_np(pi_B_tv.copy(), PRIME)) if pi_B_tv.size > 0 else 0

    # Omega block structure: how many Omega_3 basis vectors are pure old, pure tv, mixed?
    pure_old_omega = 0
    pure_tv_omega = 0
    mixed_omega = 0
    for row in ob[3]:
        has_old = any(row[i] != 0 for i in old_idx[3])
        has_tv = any(row[i] != 0 for i in tv_idx[3])
        if has_old and has_tv:
            mixed_omega += 1
        elif has_old:
            pure_old_omega += 1
        elif has_tv:
            pure_tv_omega += 1

    # Key dimension data
    n_old = {deg: len(old_idx[deg]) for deg in range(6)}
    n_tv = {deg: len(tv_idx[deg]) for deg in range(6)}

    # Constraint matrix block structure for degree 3
    P3, nr3, nc3 = _build_constraint_matrix(ap, 3, PRIME)
    if P3 is not None:
        P3_np = np.array(P3, dtype=np.int64) if isinstance(P3, list) else P3
        # Identify which constraint rows correspond to old/tv non-allowed faces
        # We need to reconstruct the non-allowed faces
        apm1_set = set(ap.get(2, []))
        na_faces = []
        seen = {}
        for path in paths[3]:
            for sign, face in boundary_faces(path):
                if len(set(face)) == len(face) and face not in apm1_set:
                    if face not in seen:
                        seen[face] = len(na_faces)
                        na_faces.append(face)

        na_old = [i for i, f in enumerate(na_faces) if v not in f]
        na_tv = [i for i, f in enumerate(na_faces) if v in f]

        # Block ranks of P3
        if na_old and old_idx[3]:
            P3_oo_rk = int(_gauss_rank_np(P3_np[np.ix_(na_old, old_idx[3])].copy(), PRIME))
        else:
            P3_oo_rk = 0
        if na_old and tv_idx[3]:
            P3_ot_rk = int(_gauss_rank_np(P3_np[np.ix_(na_old, tv_idx[3])].copy(), PRIME))
        else:
            P3_ot_rk = 0
        if na_tv and old_idx[3]:
            P3_to = P3_np[np.ix_(na_tv, old_idx[3])]
            P3_to_nnz = int(np.count_nonzero(P3_to))
        else:
            P3_to_nnz = 0
        if na_tv and tv_idx[3]:
            P3_tt_rk = int(_gauss_rank_np(P3_np[np.ix_(na_tv, tv_idx[3])].copy(), PRIME))
        else:
            P3_tt_rk = 0
    else:
        na_old = []
        na_tv = []
        P3_oo_rk = 0
        P3_ot_rk = 0
        P3_to_nnz = 0
        P3_tt_rk = 0

    # D\v chain complex info
    dim_omega3_Tv = cc_Tv.get('dims', {}).get(3, 0) if isinstance(cc_Tv, dict) else 0

    # Omega_3(T\v) dimension from full_chain_complex output
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(n1-1, 5))
    ob3_Tv = get_omega_basis(3) if False else None  # lazy
    # Just get dim from cc_Tv
    dims_Tv = cc_Tv.get('dims', {})

    result = {
        'codim': codim,
        'beta3': beta3,
        'b3_Tv': b3_Tv,
        'dim_K': dim_K,
        'rk_im_d4': rk_im_d4,
        'rk_pi_K': rk_pi_K,
        'rk_pi_B': rk_pi_B,
        'rk_pi_K_tv': rk_pi_K_tv,
        'rk_pi_B_tv': rk_pi_B_tv,
        'd3_to_nnz': d3_to_nnz,
        'd4_to_nnz': d4_to_nnz,
        'n_old_3': n_old[3],
        'n_tv_3': n_tv[3],
        'n_old_4': n_old.get(4, 0),
        'n_tv_4': n_tv.get(4, 0),
        'dim_omega3': ob[3].shape[0],
        'pure_old_omega3': pure_old_omega,
        'pure_tv_omega3': pure_tv_omega,
        'mixed_omega3': mixed_omega,
        'P3_oo_rk': P3_oo_rk,
        'P3_to_nnz': P3_to_nnz,
        'P3_tt_rk': P3_tt_rk,
        'na_old_3': len(na_old),
        'na_tv_3': len(na_tv),
        # Kernel decomposition
        'dim_K_tv': dim_K - rk_pi_K,  # = dim of through-v-only cycles
        'dim_B_tv': rk_im_d4 - rk_pi_B,  # = dim of through-v-only boundaries
    }

    if verbose and codim != 1:
        print(f"  CODIM={codim} ANOMALY:")
        print(f"    dim_K={dim_K}, rk_im_d4={rk_im_d4}")
        print(f"    rk_pi_K={rk_pi_K}, rk_pi_B={rk_pi_B}")
        print(f"    K_tv={dim_K - rk_pi_K}, B_tv={rk_im_d4 - rk_pi_B}")

    return result


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"BLOCK DECOMPOSITION ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        t0 = time.time()
        target = 300 if n == 7 else 100

        for trial in range(50000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            for v_cand in range(n):
                r = block_analysis(A, n, v_cand, verbose=(n==8))
                if r is None:
                    continue
                results.append(r)
                if len(results) % 100 == 0:
                    print(f"  Progress: {len(results)} pairs, {time.time()-t0:.1f}s")

        elapsed = time.time() - t0
        print(f"\n  Total: {len(results)} (T,v) pairs, {elapsed:.1f}s")

        if not results:
            print("  No b3=1 tournaments found!")
            continue

        # === CORE RESULT: codim distribution ===
        codim_dist = Counter(r['codim'] for r in results)
        print(f"\n  CODIM DISTRIBUTION: {dict(sorted(codim_dist.items()))}")

        # === BLOCK TRIANGULARITY VERIFICATION ===
        d3_to_fails = sum(1 for r in results if r['d3_to_nnz'] > 0)
        d4_to_fails = sum(1 for r in results if r['d4_to_nnz'] > 0)
        print(f"\n  BLOCK TRIANGULARITY:")
        print(f"    D3_to nonzero: {d3_to_fails}/{len(results)} (should be 0)")
        print(f"    D4_to nonzero: {d4_to_fails}/{len(results)} (should be 0)")

        # === CONSTRAINT MATRIX STRUCTURE ===
        p3_to_fails = sum(1 for r in results if r['P3_to_nnz'] > 0)
        print(f"\n  CONSTRAINT MATRIX P3:")
        print(f"    P3_to nonzero: {p3_to_fails}/{len(results)} (block triangularity?)")
        print(f"    P3_oo rank: mean={np.mean([r['P3_oo_rk'] for r in results]):.1f}")
        print(f"    P3_tt rank: mean={np.mean([r['P3_tt_rk'] for r in results]):.1f}")

        # === OMEGA3 DECOMPOSITION ===
        print(f"\n  OMEGA_3 DECOMPOSITION:")
        print(f"    dim(Omega_3): mean={np.mean([r['dim_omega3'] for r in results]):.1f}")
        print(f"    pure old: mean={np.mean([r['pure_old_omega3'] for r in results]):.1f}")
        print(f"    pure tv: mean={np.mean([r['pure_tv_omega3'] for r in results]):.1f}")
        print(f"    mixed: mean={np.mean([r['mixed_omega3'] for r in results]):.1f}")

        # === RAW PATH COUNTS ===
        print(f"\n  RAW PATH COUNTS (deg 3):")
        print(f"    old: mean={np.mean([r['n_old_3'] for r in results]):.1f}")
        print(f"    tv: mean={np.mean([r['n_tv_3'] for r in results]):.1f}")
        print(f"  RAW PATH COUNTS (deg 4):")
        print(f"    old: mean={np.mean([r['n_old_4'] for r in results]):.1f}")
        print(f"    tv: mean={np.mean([r['n_tv_4'] for r in results]):.1f}")

        # === KEY: TV projection rank comparison ===
        print(f"\n  KEY RANKS:")
        print(f"    rk(pi_old(K)): mean={np.mean([r['rk_pi_K'] for r in results]):.1f}")
        print(f"    rk(pi_old(B)): mean={np.mean([r['rk_pi_B'] for r in results]):.1f}")
        print(f"    rk(pi_tv(K)): mean={np.mean([r['rk_pi_K_tv'] for r in results]):.1f}")
        print(f"    rk(pi_tv(B)): mean={np.mean([r['rk_pi_B_tv'] for r in results]):.1f}")

        # === K_tv vs B_tv ===
        K_tv_vals = [r['dim_K_tv'] for r in results]
        B_tv_vals = [r['dim_B_tv'] for r in results]
        print(f"\n  THROUGH-V-ONLY DIMENSIONS:")
        print(f"    dim(K_tv): dist={dict(Counter(K_tv_vals))}")
        print(f"    dim(B_tv): dist={dict(Counter(B_tv_vals))}")
        gaps = [r['dim_K_tv'] - r['dim_B_tv'] for r in results]
        print(f"    K_tv - B_tv: dist={dict(Counter(gaps))}")

        # === CRUCIAL: Look for rank identity ===
        # codim = rk_pi_K - rk_pi_B = 1 - (K_tv - B_tv)
        # We want to find WHY K_tv = B_tv
        print(f"\n  RANK IDENTITY SEARCH:")
        # Compare dim_K with Tv-related quantities
        print(f"    b3(T\\v): dist={dict(Counter(r['b3_Tv'] for r in results))}")

        # Does rk_pi_K = dim_K - dim(K_tv)?
        check1 = all(r['rk_pi_K'] == r['dim_K'] - r['dim_K_tv'] for r in results)
        print(f"    rk_pi_K = dim_K - K_tv: {check1}")

        # Does rk_pi_B = rk_im_d4 - dim(B_tv)?
        check2 = all(r['rk_pi_B'] == r['rk_im_d4'] - r['dim_B_tv'] for r in results)
        print(f"    rk_pi_B = rk_im_d4 - B_tv: {check2}")

        # Does rk_pi_K_tv + rk_pi_K = dim_K? (sanity)
        check3 = all(r['rk_pi_K_tv'] + r['rk_pi_K'] >= r['dim_K'] for r in results)
        print(f"    rk_pi_K_tv + rk_pi_K >= dim_K: {check3}")

        # Look for new identities involving mixed vs pure omega counts
        print(f"\n  OMEGA STRUCTURE vs CODIM:")
        for codim_val in sorted(codim_dist.keys()):
            subset = [r for r in results if r['codim'] == codim_val]
            if subset:
                print(f"    codim={codim_val} ({len(subset)} cases):")
                print(f"      dim_K={np.mean([r['dim_K'] for r in subset]):.1f}, "
                      f"rk_im_d4={np.mean([r['rk_im_d4'] for r in subset]):.1f}")
                print(f"      K_tv={np.mean([r['dim_K_tv'] for r in subset]):.1f}, "
                      f"B_tv={np.mean([r['dim_B_tv'] for r in subset]):.1f}")
                print(f"      mixed_omega3={np.mean([r['mixed_omega3'] for r in subset]):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
