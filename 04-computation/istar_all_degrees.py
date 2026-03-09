"""
istar_all_degrees.py — Compute rank(i_*^p) for ALL degrees p

For beta_3=1 tournaments, computes rank of the inclusion-induced map
  i_*^p: H_p(T\\v) -> H_p(T)
at every degree, using mod-p arithmetic.

From the LES, once we know all rank(i_*^p), we can derive all H_p(T,T\\v).

Method: For each p, build the cycle space Z_p(T\\v) embedded in Omega_p(T),
compute its image in H_p(T) = Z_p(T) / B_p(T), which is:
  rank(i_*^p) = rank([B_p(T) | Z_p(T\\v)_embedded]) - rank(B_p(T))

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


def compute_istar_all_degrees(A, n, v, max_p=6):
    """Compute rank(i_*^p) for all degrees p.

    Returns dict: {p: rank_istar} for p = 0, ..., max_p.
    Also returns bettis_T and bettis_Tv for LES computation.
    """
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    # Full chain complex data for T and T\v
    data_T = full_chain_complex_modp(A, n, max_p)
    data_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    bT = data_T['bettis']
    bTv = data_Tv['bettis']

    # For rank(i_*^p), we need:
    # 1. The cycle space Z_p(T\v) in Omega_p(T\v) coords
    # 2. The embedding Omega_p(T\v) -> Omega_p(T) in A-path coords
    # 3. The boundary space B_p(T) = im(d_{p+1}) in Omega_p(T) coords
    # 4. rank(i_*^p) = rank([B_p(T) columns | embedded Z_p(T\v) columns]) - rank(B_p(T))

    # We need A-path level data for both complexes
    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    rank_istar = {}

    for p in range(max_p + 1):
        bp_T = bT.get(p, 0)
        bp_Tv = bTv.get(p, 0)

        # Trivial cases
        if bp_Tv == 0:
            rank_istar[p] = 0
            continue
        if bp_T == 0:
            rank_istar[p] = 0
            continue

        # Both H_p(T) and H_p(T\v) are nonzero.
        # Need to compute rank of the induced map.

        paths_T_p = ap_T.get(p, [])
        paths_Tv_p = ap_Tv.get(p, [])
        paths_T_pm1 = ap_T.get(p - 1, []) if p > 0 else []
        paths_T_pp1 = ap_T.get(p + 1, [])
        paths_Tv_pm1 = ap_Tv.get(p - 1, []) if p > 0 else []
        paths_Tv_pp1 = ap_Tv.get(p + 1, [])

        if not paths_T_p or not paths_Tv_p:
            rank_istar[p] = 0
            continue

        # Build Omega basis for T at degree p
        P_T, na_rows_T, na_cols_T = _build_constraint_matrix(ap_T, p, PRIME)
        if P_T is not None:
            rank_P_T, nbasis_T = _gauss_nullbasis_modp(P_T, na_rows_T, na_cols_T, PRIME)
            omega_basis_T = np.array(nbasis_T, dtype=np.int64) if nbasis_T else None
        else:
            omega_basis_T = np.eye(len(paths_T_p), dtype=np.int64)

        if omega_basis_T is None or len(omega_basis_T) == 0:
            rank_istar[p] = 0
            continue

        # Build Omega basis for T\v at degree p
        P_Tv, na_rows_Tv, na_cols_Tv = _build_constraint_matrix(ap_Tv, p, PRIME)
        if P_Tv is not None:
            rank_P_Tv, nbasis_Tv = _gauss_nullbasis_modp(P_Tv, na_rows_Tv, na_cols_Tv, PRIME)
            omega_basis_Tv = np.array(nbasis_Tv, dtype=np.int64) if nbasis_Tv else None
        else:
            omega_basis_Tv = np.eye(len(paths_Tv_p), dtype=np.int64)

        if omega_basis_Tv is None or len(omega_basis_Tv) == 0:
            rank_istar[p] = 0
            continue

        # Build boundary matrix d_p: A_p(T) -> A_{p-1}(T)
        if p > 0 and paths_T_pm1:
            idx_Tpm1 = {tuple(q): i for i, q in enumerate(paths_T_pm1)}
            bd_T_p = np.zeros((len(paths_T_pm1), len(paths_T_p)), dtype=np.int64)
            for j, sigma in enumerate(paths_T_p):
                for sign, face in boundary_faces(sigma):
                    face_t = tuple(face)
                    if face_t in idx_Tpm1:
                        bd_T_p[idx_Tpm1[face_t], j] = (bd_T_p[idx_Tpm1[face_t], j] + sign) % PRIME
        else:
            bd_T_p = None

        # d_p in Omega coords: bd_T_p @ omega_basis_T^T
        if bd_T_p is not None:
            d_p_omega = bd_T_p @ omega_basis_T.T % PRIME
        else:
            d_p_omega = None

        # Cycle space Z_p(T) in Omega coords = ker(d_p_omega)
        # Z_p(T) vectors are in omega_basis_T space: z such that d_p_omega @ z^T = 0
        # Actually d_p_omega is (|A_{p-1}| x omega_dim_T_p), so ker is the null space.
        if d_p_omega is not None:
            # Use Gauss to find null space of d_p_omega (as map from R^omega_dim to R^|A_{p-1}|)
            # Transpose: (omega_dim x |A_{p-1}|) and find left null space... no.
            # d_p_omega is (|A_{p-1}| x omega_dim_T). Null space = vectors v with d_p_omega @ v = 0.
            # That's the right null space of d_p_omega.
            pass  # Will use a different approach

        # Build boundary matrix d_{p+1}: A_{p+1}(T) -> A_p(T)
        if paths_T_pp1:
            idx_Tp = {tuple(q): i for i, q in enumerate(paths_T_p)}
            bd_T_pp1 = np.zeros((len(paths_T_p), len(paths_T_pp1)), dtype=np.int64)
            for j, sigma in enumerate(paths_T_pp1):
                for sign, face in boundary_faces(sigma):
                    face_t = tuple(face)
                    if face_t in idx_Tp:
                        bd_T_pp1[idx_Tp[face_t], j] = (bd_T_pp1[idx_Tp[face_t], j] + sign) % PRIME

            # d_{p+1} in Omega coords of degree p:
            # im(d_{p+1}) in Omega_p(T) coords needs omega basis for p+1 too
            P_Tpp1, na_rows_pp1, na_cols_pp1 = _build_constraint_matrix(ap_T, p + 1, PRIME)
            if P_Tpp1 is not None:
                rank_pp1, nbasis_pp1 = _gauss_nullbasis_modp(P_Tpp1, na_rows_pp1, na_cols_pp1, PRIME)
                omega_basis_T_pp1 = np.array(nbasis_pp1, dtype=np.int64) if nbasis_pp1 else None
            else:
                omega_basis_T_pp1 = np.eye(len(paths_T_pp1), dtype=np.int64)

            if omega_basis_T_pp1 is not None and len(omega_basis_T_pp1) > 0:
                # im(d_{p+1}) in A_p coords = bd_T_pp1 @ omega_basis_T_pp1^T
                im_dp1_Acoords = bd_T_pp1 @ omega_basis_T_pp1.T % PRIME  # (|A_p| x omega_dim_{p+1})
            else:
                im_dp1_Acoords = np.zeros((len(paths_T_p), 1), dtype=np.int64)
        else:
            im_dp1_Acoords = np.zeros((len(paths_T_p), 1), dtype=np.int64)

        # Embed Omega_p(T\v) into A_p(T) coordinates
        # omega_basis_Tv rows are in A_p(T\v) coords. Embed to A_p(T) coords.
        idx_Tp = {tuple(q): i for i, q in enumerate(paths_T_p)}
        embed_matrix = np.zeros((len(paths_T_p), len(paths_Tv_p)), dtype=np.int64)
        for j, path_tv in enumerate(paths_Tv_p):
            embedded = embed_path(path_tv)
            if embedded in idx_Tp:
                embed_matrix[idx_Tp[embedded], j] = 1

        # Omega_p(T\v) embedded in A_p(T) coords:
        # columns = embed_matrix @ omega_basis_Tv^T  (|A_p(T)| x omega_dim_Tv)
        embedded_omega_Tv = embed_matrix @ omega_basis_Tv.T % PRIME

        # Build boundary d_p for T\v in A_p(T\v) coords
        if p > 0 and paths_Tv_pm1:
            idx_Tvpm1 = {tuple(q): i for i, q in enumerate(paths_Tv_pm1)}
            bd_Tv_p = np.zeros((len(paths_Tv_pm1), len(paths_Tv_p)), dtype=np.int64)
            for j, sigma in enumerate(paths_Tv_p):
                for sign, face in boundary_faces(sigma):
                    face_t = tuple(face)
                    if face_t in idx_Tvpm1:
                        bd_Tv_p[idx_Tvpm1[face_t], j] = (bd_Tv_p[idx_Tvpm1[face_t], j] + sign) % PRIME

            # d_p(T\v) in Omega coords: bd_Tv_p @ omega_basis_Tv^T  (rows in A_{p-1}(T\v))
            d_p_Tv_omega = bd_Tv_p @ omega_basis_Tv.T % PRIME

            # Cycle space Z_p(T\v) in Omega_p(T\v) coords = null space of d_p_Tv_omega
            # d_p_Tv_omega is (|A_{p-1}(T\v)| x omega_dim_Tv_p), null space = right null
            nrows_dTv = d_p_Tv_omega.shape[0]
            ncols_dTv = d_p_Tv_omega.shape[1]
            M_list = [[int(x) % PRIME for x in d_p_Tv_omega[r]] for r in range(nrows_dTv)]
            rank_dTv, zycles_Tv = _gauss_nullbasis_modp(
                M_list, nrows_dTv, ncols_dTv, PRIME
            )
            if zycles_Tv:
                zycles_Tv = np.array(zycles_Tv, dtype=np.int64)  # (ker_dim x omega_dim_Tv_p)
            else:
                rank_istar[p] = 0
                continue
        else:
            # p=0: all of Omega_0(T\v) is cycles
            zycles_Tv = omega_basis_Tv

        # Embedded cycles: embed_matrix @ omega_basis_Tv^T @ zycles_Tv^T
        # = embedded_omega_Tv @ zycles_Tv^T  (|A_p(T)| x ker_dim)
        embedded_cycles = embedded_omega_Tv @ zycles_Tv.T % PRIME

        # Now rank(i_*^p) = rank([im_dp1_Acoords | embedded_cycles]) - rank(im_dp1_Acoords)
        # Both are in A_p(T) coordinates.
        rank_B = _gauss_rank_np(im_dp1_Acoords.copy() % PRIME, PRIME)
        combined = np.hstack([im_dp1_Acoords, embedded_cycles]) % PRIME
        rank_combined = _gauss_rank_np(combined, PRIME)
        rank_istar[p] = rank_combined - rank_B

    return rank_istar, bT, bTv


def compute_relative_bettis_from_les(rank_istar, bT, bTv, max_p=6):
    """Use the LES to compute H_p(T,T\\v) from rank(i_*).

    LES: ... H_p(T\\v) -i_*-> H_p(T) -j_*-> H_p(T,T\\v) -delta-> H_{p-1}(T\\v) ...

    Exactness gives:
      dim H_p(T,T\\v) = dim(coker i_*^p) + dim(ker i_*^{p-1})
                       = (b_p(T) - rank(i_*^p)) + (b_{p-1}(T\\v) - rank(i_*^{p-1}))
    """
    rel_bettis = {}
    for p in range(max_p + 1):
        coker_ip = bT.get(p, 0) - rank_istar.get(p, 0)
        ker_ip_m1 = bTv.get(p - 1, 0) - rank_istar.get(p - 1, 0) if p > 0 else 0
        rel_bettis[p] = coker_ip + ker_ip_m1
    return rel_bettis


def main():
    print("=" * 70)
    print("RANK(i_*) AT ALL DEGREES — BETA_3=1 TOURNAMENTS")
    print("=" * 70)

    n = 7
    max_p = n - 1
    num_samples = 80

    print(f"\n--- n={n}, sampling {num_samples} beta_3=1 tournaments ---")

    found = 0
    t0 = time.time()

    # Collect: for each (vertex_type, p): list of rank_istar values
    istar_data = []  # list of (vtype, rank_istar_dict, bT, bTv)
    rel_betti_profiles = Counter()

    while found < num_samples:
        A = random_tournament(n)
        cc = full_chain_complex_modp(A, n, max_p)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            ri, bT, bTv = compute_istar_all_degrees(A, n, v, max_p)
            b3_sub = bTv.get(3, 0)
            vtype = "BAD" if b3_sub == 1 else "GOOD"

            rel = compute_relative_bettis_from_les(ri, bT, bTv, max_p)
            istar_data.append((vtype, ri, bT, bTv, rel))

            rel_profile = tuple(rel.get(p, 0) for p in range(max_p + 1))
            rel_betti_profiles[(vtype, rel_profile)] += 1

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{num_samples} found, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tournaments, {len(istar_data)} vertex records, {elapsed:.1f}s")

    # Analyze rank(i_*) by vertex type and degree
    for vtype in ["GOOD", "BAD"]:
        recs = [(ri, bT, bTv, rel) for vt, ri, bT, bTv, rel in istar_data if vt == vtype]
        print(f"\n  {vtype} vertices ({len(recs)}):")

        for p in range(max_p + 1):
            vals = [ri.get(p, 0) for ri, _, _, _ in recs]
            bT_vals = [bT.get(p, 0) for _, bT, _, _ in recs]
            bTv_vals = [bTv.get(p, 0) for _, _, bTv, _ in recs]
            dist = Counter(vals)
            print(f"    p={p}: rank(i_*)={dict(dist)}, "
                  f"b_p(T) range=[{min(bT_vals)},{max(bT_vals)}], "
                  f"b_p(T\\v) range=[{min(bTv_vals)},{max(bTv_vals)}]")

    # Show relative Betti profiles
    print(f"\n  Relative Betti profiles H_*(T,T\\v):")
    for (vtype, profile), count in rel_betti_profiles.most_common(20):
        print(f"    {vtype}: {profile}: {count}")

    # Check specific claims
    print(f"\n  CLAIM VERIFICATION:")
    h4_violations = 0
    h3_gt1 = 0
    for vtype, ri, bT, bTv, rel in istar_data:
        if rel.get(4, 0) != 0:
            h4_violations += 1
        if rel.get(3, 0) > 1:
            h3_gt1 += 1
    print(f"    H_4(T,T\\v) = 0 violations: {h4_violations}")
    print(f"    H_3(T,T\\v) > 1 violations: {h3_gt1}")

    # Check if rank(i_*) has universal values
    print(f"\n  RANK(i_*) UNIVERSALITY CHECK:")
    for p in range(max_p + 1):
        good_vals = set(ri.get(p, 0) for vt, ri, _, _, _ in istar_data if vt == "GOOD")
        bad_vals = set(ri.get(p, 0) for vt, ri, _, _, _ in istar_data if vt == "BAD")
        print(f"    p={p}: GOOD i_* in {good_vals}, BAD i_* in {bad_vals}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
