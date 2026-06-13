"""
istar_n8_investigation.py — rank(i_*) and relative Betti at n=8

Critical test: consecutive seesaw (HYP-394) FAILS at n=8 (beta_3=beta_4=1 exists).
Does the proof architecture still work?

Key questions:
1. Is rank(i_*^3) still {0,1} at n=8? (i_*-injectivity)
2. Is H_3(T,T\\v) still ≤ 1?
3. Is H_4(T,T\\v) = 0 at n=8? (This is NO LONGER guaranteed by seesaw)
4. What happens for beta_3=beta_4=1 tournaments specifically?

Author: opus-2026-03-09-S55
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed, boundary_faces,
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    full_chain_complex_modp,
    RANK_PRIME
)

PRIME = RANK_PRIME


def compute_istar_all_degrees(A, n, v, max_p=7):
    """Compute rank(i_*^p) for all degrees p. Returns (rank_istar, bT, bTv)."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    data_T = full_chain_complex_modp(A, n, max_p)
    data_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    bT = data_T['bettis']
    bTv = data_Tv['bettis']

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    rank_istar = {}

    for p in range(max_p + 1):
        bp_T = bT.get(p, 0)
        bp_Tv = bTv.get(p, 0)

        if bp_Tv == 0 or bp_T == 0:
            rank_istar[p] = 0
            continue

        paths_T_p = ap_T.get(p, [])
        paths_Tv_p = ap_Tv.get(p, [])
        paths_T_pp1 = ap_T.get(p + 1, [])
        paths_Tv_pm1 = ap_Tv.get(p - 1, []) if p > 0 else []

        if not paths_T_p or not paths_Tv_p:
            rank_istar[p] = 0
            continue

        # Omega basis for T at degree p
        P_T, na_rows_T, na_cols_T = _build_constraint_matrix(ap_T, p, PRIME)
        if P_T is not None:
            rank_P_T, nbasis_T = _gauss_nullbasis_modp(P_T, na_rows_T, na_cols_T, PRIME)
            omega_basis_T = np.array(nbasis_T, dtype=np.int64) if nbasis_T else None
        else:
            omega_basis_T = np.eye(len(paths_T_p), dtype=np.int64)

        if omega_basis_T is None or len(omega_basis_T) == 0:
            rank_istar[p] = 0
            continue

        # Omega basis for T\v at degree p
        P_Tv, na_rows_Tv, na_cols_Tv = _build_constraint_matrix(ap_Tv, p, PRIME)
        if P_Tv is not None:
            rank_P_Tv, nbasis_Tv = _gauss_nullbasis_modp(P_Tv, na_rows_Tv, na_cols_Tv, PRIME)
            omega_basis_Tv = np.array(nbasis_Tv, dtype=np.int64) if nbasis_Tv else None
        else:
            omega_basis_Tv = np.eye(len(paths_Tv_p), dtype=np.int64)

        if omega_basis_Tv is None or len(omega_basis_Tv) == 0:
            rank_istar[p] = 0
            continue

        # Boundary im(d_{p+1}^T) in A_p(T) coords
        idx_Tp = {tuple(q): i for i, q in enumerate(paths_T_p)}
        if paths_T_pp1:
            bd_T_pp1 = np.zeros((len(paths_T_p), len(paths_T_pp1)), dtype=np.int64)
            for j, sigma in enumerate(paths_T_pp1):
                for sign, face in boundary_faces(sigma):
                    face_t = tuple(face)
                    if face_t in idx_Tp:
                        bd_T_pp1[idx_Tp[face_t], j] = (bd_T_pp1[idx_Tp[face_t], j] + sign) % PRIME

            P_Tpp1, na_rows_pp1, na_cols_pp1 = _build_constraint_matrix(ap_T, p + 1, PRIME)
            if P_Tpp1 is not None:
                r_pp1, nb_pp1 = _gauss_nullbasis_modp(P_Tpp1, na_rows_pp1, na_cols_pp1, PRIME)
                omega_basis_T_pp1 = np.array(nb_pp1, dtype=np.int64) if nb_pp1 else None
            else:
                omega_basis_T_pp1 = np.eye(len(paths_T_pp1), dtype=np.int64)

            if omega_basis_T_pp1 is not None and len(omega_basis_T_pp1) > 0:
                im_dp1 = bd_T_pp1 @ omega_basis_T_pp1.T % PRIME
            else:
                im_dp1 = np.zeros((len(paths_T_p), 1), dtype=np.int64)
        else:
            im_dp1 = np.zeros((len(paths_T_p), 1), dtype=np.int64)

        # Embed Omega_p(T\v) into A_p(T)
        embed_matrix = np.zeros((len(paths_T_p), len(paths_Tv_p)), dtype=np.int64)
        for j, path_tv in enumerate(paths_Tv_p):
            embedded = embed_path(path_tv)
            if embedded in idx_Tp:
                embed_matrix[idx_Tp[embedded], j] = 1

        embedded_omega_Tv = embed_matrix @ omega_basis_Tv.T % PRIME

        # Build cycle space Z_p(T\v)
        if p > 0 and paths_Tv_pm1:
            idx_Tvpm1 = {tuple(q): i for i, q in enumerate(paths_Tv_pm1)}
            bd_Tv_p = np.zeros((len(paths_Tv_pm1), len(paths_Tv_p)), dtype=np.int64)
            for j, sigma in enumerate(paths_Tv_p):
                for sign, face in boundary_faces(sigma):
                    face_t = tuple(face)
                    if face_t in idx_Tvpm1:
                        bd_Tv_p[idx_Tvpm1[face_t], j] = (bd_Tv_p[idx_Tvpm1[face_t], j] + sign) % PRIME

            d_p_Tv_omega = bd_Tv_p @ omega_basis_Tv.T % PRIME
            nrows_dTv = d_p_Tv_omega.shape[0]
            ncols_dTv = d_p_Tv_omega.shape[1]
            M_list = [[int(x) % PRIME for x in d_p_Tv_omega[r]] for r in range(nrows_dTv)]
            rank_dTv, zycles_Tv = _gauss_nullbasis_modp(M_list, nrows_dTv, ncols_dTv, PRIME)
            if zycles_Tv:
                zycles_Tv = np.array(zycles_Tv, dtype=np.int64)
            else:
                rank_istar[p] = 0
                continue
        else:
            zycles_Tv = omega_basis_Tv

        # Embedded cycles in A_p(T) coords
        embedded_cycles = embedded_omega_Tv @ zycles_Tv.T % PRIME

        # rank(i_*^p) = rank([im_dp1 | embedded_cycles]) - rank(im_dp1)
        rank_B = _gauss_rank_np(im_dp1.copy() % PRIME, PRIME)
        combined = np.hstack([im_dp1, embedded_cycles]) % PRIME
        rank_combined = _gauss_rank_np(combined, PRIME)
        rank_istar[p] = rank_combined - rank_B

    return rank_istar, bT, bTv


def compute_relative_bettis_from_les(rank_istar, bT, bTv, max_p=7):
    """Derive H_p(T,T\\v) from LES."""
    rel_bettis = {}
    for p in range(max_p + 1):
        coker_ip = bT.get(p, 0) - rank_istar.get(p, 0)
        ker_ip_m1 = bTv.get(p - 1, 0) - rank_istar.get(p - 1, 0) if p > 0 else 0
        rel_bettis[p] = coker_ip + ker_ip_m1
    return rel_bettis


def main():
    print("=" * 70)
    print("RANK(i_*) AT ALL DEGREES — n=8 INVESTIGATION")
    print("=" * 70)

    n = 8
    max_p = n - 1
    rng = np.random.RandomState(42)

    print(f"\n--- Phase 1: beta_3=1 tournaments at n={n} ---")
    found = 0
    found_b3b4 = 0
    target = 30  # beta_3=1 tournaments
    t0 = time.time()

    istar_data = []
    rel_betti_profiles = Counter()

    trials = 0
    while found < target:
        A = random_tournament(n, rng)
        trials += 1
        cc = full_chain_complex_modp(A, n, max_p)
        b3 = cc['bettis'].get(3, 0)
        if b3 != 1:
            continue
        found += 1
        b4 = cc['bettis'].get(4, 0)
        if b4 > 0:
            found_b3b4 += 1

        profile_T = tuple(cc['bettis'].get(p, 0) for p in range(max_p + 1))

        # Check a sample of vertices (not all 8 — too slow)
        for v in range(min(n, 4)):
            ri, bT, bTv = compute_istar_all_degrees(A, n, v, max_p)
            b3_sub = bTv.get(3, 0)
            b4_sub = bTv.get(4, 0)
            vtype = "BAD" if b3_sub == 1 else "GOOD"

            rel = compute_relative_bettis_from_les(ri, bT, bTv, max_p)
            istar_data.append((vtype, ri, bT, bTv, rel, b4 > 0))

            rel_profile = tuple(rel.get(p, 0) for p in range(max_p + 1))
            rel_betti_profiles[(vtype, rel_profile, b4 > 0)] += 1

        if found % 5 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target} beta_3=1 found ({trials} trials), "
                  f"{found_b3b4} with b4>0, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments ({trials} trials), "
          f"{len(istar_data)} vertex records, {elapsed:.1f}s")
    print(f"  beta_3=beta_4=1 tournaments: {found_b3b4}")

    # Analyze
    for vtype in ["GOOD", "BAD"]:
        recs = [(ri, bT, bTv, rel, b4flag) for vt, ri, bT, bTv, rel, b4flag in istar_data if vt == vtype]
        print(f"\n  {vtype} vertices ({len(recs)}):")
        for p in range(max_p + 1):
            vals = Counter(ri.get(p, 0) for ri, _, _, _, _ in recs)
            print(f"    p={p}: rank(i_*)={dict(vals)}")

    # Relative Betti profiles
    print(f"\n  Relative Betti profiles H_*(T,T\\v):")
    for (vtype, profile, b4flag), count in rel_betti_profiles.most_common(30):
        b4str = " [b4>0]" if b4flag else ""
        print(f"    {vtype}: {profile}: {count}{b4str}")

    # Claim verification
    print(f"\n  CLAIM VERIFICATION:")
    h4_violations = 0
    h3_gt1 = 0
    for vtype, ri, bT, bTv, rel, b4flag in istar_data:
        if rel.get(4, 0) != 0:
            h4_violations += 1
            print(f"    H_4^rel VIOLATION: {vtype}, betti_T={tuple(bT.get(p,0) for p in range(8))}, "
                  f"rel={tuple(rel.get(p,0) for p in range(8))}")
        if rel.get(3, 0) > 1:
            h3_gt1 += 1
    print(f"    H_4(T,T\\v) = 0 violations: {h4_violations}")
    print(f"    H_3(T,T\\v) > 1 violations: {h3_gt1}")

    # Universality check
    print(f"\n  RANK(i_*) UNIVERSALITY:")
    for p in range(max_p + 1):
        good_vals = set(ri.get(p, 0) for vt, ri, _, _, _, _ in istar_data if vt == "GOOD")
        bad_vals = set(ri.get(p, 0) for vt, ri, _, _, _, _ in istar_data if vt == "BAD")
        print(f"    p={p}: GOOD i_* in {good_vals}, BAD i_* in {bad_vals}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
