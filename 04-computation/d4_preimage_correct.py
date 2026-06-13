"""
d4_preimage_correct.py — Corrected: use actual H_3 generator, not just any cycle

Bug in previous script: ker_d3Tv[0] might be a boundary of T\\v.
The correct z is a cycle NOT in im(d_4^{T\\v}).

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


def get_omega_basis(ap, deg, prime):
    paths = ap.get(deg, [])
    if not paths:
        return None
    P, nr, nc = _build_constraint_matrix(ap, deg, prime)
    if P is not None:
        r, nb = _gauss_nullbasis_modp(P, nr, nc, prime)
        return np.array(nb, dtype=np.int64) if nb else None
    else:
        return np.eye(len(paths), dtype=np.int64)


def find_h3_generator(ap, ob3, ob4, paths3, paths2, paths4):
    """Find a cycle z in ker(d_3) that is NOT in im(d_4).
    Returns z in A_3 coordinates."""

    if ob3 is None or ob3.shape[0] == 0:
        return None

    # d_3 in Omega coords: bd_3 @ ob3.T
    idx2 = {tuple(q): i for i, q in enumerate(paths2)}
    bd3 = np.zeros((len(paths2), len(paths3)), dtype=np.int64)
    for j, sigma in enumerate(paths3):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx2:
                bd3[idx2[ft], j] = (bd3[idx2[ft], j] + sign) % PRIME

    d3_omega = bd3 @ ob3.T % PRIME  # A_2 x dim_omega3

    # ker(d_3) in Omega_3 coords
    M_d3 = [[int(x) % PRIME for x in d3_omega[r]] for r in range(d3_omega.shape[0])]
    _, ker_d3 = _gauss_nullbasis_modp(M_d3, d3_omega.shape[0], d3_omega.shape[1], PRIME)
    if not ker_d3:
        return None
    ker_d3 = np.array(ker_d3, dtype=np.int64)  # (ker_dim x dim_omega3)

    # im(d_4) in A_3 coords
    if ob4 is not None and ob4.shape[0] > 0 and paths4:
        idx3 = {tuple(q): i for i, q in enumerate(paths3)}
        bd4 = np.zeros((len(paths3), len(paths4)), dtype=np.int64)
        for j, sigma in enumerate(paths4):
            for sign, face in boundary_faces(sigma):
                ft = tuple(face)
                if ft in idx3:
                    bd4[idx3[ft], j] = (bd4[idx3[ft], j] + sign) % PRIME
        im_d4 = bd4 @ ob4.T % PRIME  # A_3 x dim_omega4
    else:
        im_d4 = np.zeros((len(paths3), 1), dtype=np.int64)

    rank_im = _gauss_rank_np(im_d4.copy() % PRIME, PRIME)

    # Find a cycle NOT in im(d_4)
    # Convert each ker_d3 vector to A_3 coords: ob3.T @ ker_d3[i]
    for i in range(ker_d3.shape[0]):
        z_omega = ker_d3[i]
        z_Apath = (ob3.T @ z_omega.reshape(-1, 1)).flatten() % PRIME

        # Check if z is in im(d_4)
        z_col = z_Apath.reshape(-1, 1) % PRIME
        combined = np.hstack([im_d4, z_col]) % PRIME
        rank_combined = _gauss_rank_np(combined, PRIME)
        if rank_combined > rank_im:
            # z is NOT in im(d_4) — this is an H_3 generator!
            return z_Apath

    # All cycles are boundaries?! beta_3 = 0
    return None


def analyze_d4_correct(A, n, v, verbose=False):
    """Correct analysis: use actual H_3 generator."""
    max_p = n - 1
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]

    cc_T = full_chain_complex_modp(A, n, max_p)
    cc_Tv = full_chain_complex_modp(A_sub, n1, min(max_p, n1 - 1))

    if cc_T['bettis'].get(3, 0) != 1 or cc_Tv['bettis'].get(3, 0) != 1:
        return None

    ap_T = enumerate_all_allowed(A, n, max_p)
    ap_Tv = enumerate_all_allowed(A_sub, n1, min(max_p, n1 - 1))

    def embed_path(path_tv):
        return tuple(remaining[i] for i in path_tv)

    ob_T3 = get_omega_basis(ap_T, 3, PRIME)
    ob_T4 = get_omega_basis(ap_T, 4, PRIME)
    ob_Tv3 = get_omega_basis(ap_Tv, 3, PRIME)
    ob_Tv4 = get_omega_basis(ap_Tv, 4, PRIME)

    paths_T_3 = ap_T.get(3, [])
    paths_T_4 = ap_T.get(4, [])
    paths_Tv_2 = ap_Tv.get(2, [])
    paths_Tv_3 = ap_Tv.get(3, [])
    paths_Tv_4 = ap_Tv.get(4, [])

    # Find ACTUAL H_3(T\v) generator (not just any cycle!)
    z_Tv = find_h3_generator(ap_Tv, ob_Tv3, ob_Tv4, paths_Tv_3, paths_Tv_2, paths_Tv_4)
    if z_Tv is None:
        return None

    # Embed z in A_3(T) coords
    idx_T3 = {tuple(q): i for i, q in enumerate(paths_T_3)}
    z_T = np.zeros(len(paths_T_3), dtype=np.int64)
    for j, path_tv in enumerate(paths_Tv_3):
        if z_Tv[j] % PRIME != 0:
            embedded = embed_path(path_tv)
            if embedded in idx_T3:
                z_T[idx_T3[embedded]] = z_Tv[j] % PRIME

    # im(d_4^T) in A_3(T) coords
    bd_T4 = np.zeros((len(paths_T_3), len(paths_T_4)), dtype=np.int64)
    for j, sigma in enumerate(paths_T_4):
        for sign, face in boundary_faces(sigma):
            ft = tuple(face)
            if ft in idx_T3:
                bd_T4[idx_T3[ft], j] = (bd_T4[idx_T3[ft], j] + sign) % PRIME

    if ob_T4 is not None and ob_T4.shape[0] > 0:
        im_d4_T = bd_T4 @ ob_T4.T % PRIME
    else:
        im_d4_T = np.zeros((len(paths_T_3), 1), dtype=np.int64)

    rank_im_T = _gauss_rank_np(im_d4_T.copy() % PRIME, PRIME)

    # Is z_T in im(d_4^T)?
    z_col = z_T.reshape(-1, 1) % PRIME
    combined = np.hstack([im_d4_T, z_col]) % PRIME
    rank_combined = _gauss_rank_np(combined, PRIME)
    z_in_im = (rank_combined == rank_im_T)

    # Also find H_3(T) generator and check alignment
    z_T_gen = find_h3_generator(ap_T, ob_T3, ob_T4, paths_T_3, ap_T.get(2, []), paths_T_4)

    # Check if embedded z_T is proportional to H_3(T) generator mod im(d_4^T)
    # Since beta_3 = 1, H_3(T) is 1-dim. If z_T is nonzero in H_3(T),
    # then [z_T] = c * [z_T_gen] for some scalar c.
    # Proportional modulo im(d_4^T): z_T - c * z_T_gen ∈ im(d_4^T)
    if z_T_gen is not None:
        # Check: is [im_d4 | z_T | z_T_gen] rank = rank(im_d4) + 1?
        # If so, z_T and z_T_gen represent the same 1-dim space modulo im(d_4)
        z_gen_col = z_T_gen.reshape(-1, 1) % PRIME
        three_cols = np.hstack([im_d4_T, z_col, z_gen_col]) % PRIME
        rank_three = _gauss_rank_np(three_cols, PRIME)
        # rank should be rank_im + 1 (both z's in same 1-dim complement)
        proportional = (rank_three == rank_im_T + 1)
    else:
        proportional = False

    result = {
        'z_in_im_d4_T': z_in_im,
        'z_proportional_to_H3T': proportional,
        'rank_im_d4': rank_im_T,
        'ker_d3': cc_T['kers'].get(3, 0),
    }

    if verbose:
        print(f"    z (true H_3 gen) in im(d_4^T)? {z_in_im}")
        print(f"    z proportional to H_3(T) gen? {proportional}")
        print(f"    rank(d_4) = {rank_im_T}, ker(d_3) = {cc_T['kers'].get(3, 0)}")

    return result


def main():
    print("=" * 70)
    print("CORRECTED d_4 PREIMAGE — USING TRUE H_3 GENERATORS")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found = 0
    target = 60
    t0 = time.time()
    results = []

    while found < target:
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        found += 1

        for v in range(n):
            r = analyze_d4_correct(A, n, v, verbose=(found <= 2))
            if r is not None:
                results.append(r)

        if found % 10 == 0:
            elapsed = time.time() - t0
            print(f"  {found}/{target}, {len(results)} bad verts, {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Done: {found} tours, {len(results)} bad verts, {elapsed:.1f}s")

    print(f"\n  CORRECTED RESULTS:")
    z_in_im = sum(1 for r in results if r['z_in_im_d4_T'])
    z_prop = sum(1 for r in results if r['z_proportional_to_H3T'])
    print(f"    z (true H_3(T\\v) gen) in im(d_4^T): {z_in_im}/{len(results)}")
    print(f"    z proportional to H_3(T) gen: {z_prop}/{len(results)}")

    if z_in_im == 0:
        print(f"\n    *** i_*-INJECTIVITY CONFIRMED: z never becomes a boundary ***")
    if z_prop == len(results):
        print(f"    *** EMBEDDED z IS the H_3(T) generator (up to scalar) ***")


if __name__ == '__main__':
    main()
    print("\nDONE.")
