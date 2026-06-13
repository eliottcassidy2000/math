"""
ghost_cycle_deep_debug.py — Deep debug of K_tv computation discrepancy

Focuses on trial=58, v=2 at n=7 where Ktv_A=3 but Ktv_B=2.

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import numpy as np
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


def deep_debug(A, n, v):
    """Deep debug of K_tv computation."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths = {deg: ap.get(deg, []) for deg in range(7)}

    old3 = sorted([i for i, p in enumerate(paths[3]) if v not in p])
    tv3 = sorted([i for i, p in enumerate(paths[3]) if v in p])

    print(f"  n_paths_3={len(paths[3])}, n_old3={len(old3)}, n_tv3={len(tv3)}")

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
    print(f"  dim_Omega3={ob3.shape[0]}, dim_Omega4={ob4.shape[0]}")

    idx2 = {p: i for i, p in enumerate(paths[2])}
    D3 = np.zeros((len(paths[2]), len(paths[3])), dtype=np.int64)
    for j, path in enumerate(paths[3]):
        for sign, face in boundary_faces(path):
            if face in idx2:
                D3[idx2[face], j] = (D3[idx2[face], j] + sign) % PRIME

    idx3 = {p: i for i, p in enumerate(paths[3])}
    D4 = np.zeros((len(paths[3]), len(paths[4])), dtype=np.int64)
    for j, path in enumerate(paths[4]):
        for sign, face in boundary_faces(path):
            if face in idx3:
                D4[idx3[face], j] = (D4[idx3[face], j] + sign) % PRIME

    # METHOD A
    print(f"\n  === METHOD A ===")
    d3o = D3 @ ob3.T % PRIME
    print(f"  d3o shape: {d3o.shape}")
    _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
    ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    K = ker_d3 @ ob3 % PRIME
    dim_K = K.shape[0]
    print(f"  dim_K={dim_K}")

    im_d4 = (D4 @ ob4.T % PRIME).T
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME))
    beta3 = dim_K - rk_d4
    print(f"  rk_d4={rk_d4}, beta3={beta3}")

    pi_K = K[:, old3] % PRIME
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME))
    Ktv_A = dim_K - rk_pi_K
    print(f"  rk_pi_K={rk_pi_K}, Ktv_A={Ktv_A}")

    # Verify: find the K_tv basis vectors explicitly
    pi_K_T = pi_K.T % PRIME
    _, ktv_basis = _gauss_nullbasis_modp(pi_K_T.astype(np.int32), pi_K_T.shape[0], pi_K_T.shape[1], PRIME)
    ktv_coeffs = np.array(ktv_basis, dtype=np.int64) if ktv_basis else np.zeros((0, dim_K), dtype=np.int64)
    Ktv_explicit = ktv_coeffs @ K % PRIME
    print(f"  Ktv explicit basis size: {Ktv_explicit.shape[0]}")

    # Check that Ktv_explicit vectors are indeed tv-only and cycles
    for i in range(Ktv_explicit.shape[0]):
        z = Ktv_explicit[i]
        old_nonzero = sum(1 for j in old3 if z[j] % PRIME != 0)
        bd = D3 @ z % PRIME
        bd_nonzero = int(np.count_nonzero(bd % PRIME))
        print(f"    K_tv vector {i}: old_nonzero={old_nonzero}, boundary_nonzero={bd_nonzero}")

    # METHOD B
    print(f"\n  === METHOD B ===")
    ob3_old = ob3[:, old3] % PRIME
    rk_ob3_old = int(_gauss_rank_np(ob3_old.copy(), PRIME))
    print(f"  rk_ob3_old={rk_ob3_old}")

    ob3_old_T = ob3_old.T % PRIME
    _, tv_omega_basis = _gauss_nullbasis_modp(
        ob3_old_T.astype(np.int32), ob3_old_T.shape[0], ob3_old_T.shape[1], PRIME)
    tv_omega = np.array(tv_omega_basis, dtype=np.int64) if tv_omega_basis else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    dim_omega3_tv = tv_omega.shape[0]
    print(f"  dim_omega3_tv={dim_omega3_tv}")

    # Check: tv_omega vectors should have zero old component
    tv_omega_A = tv_omega @ ob3 % PRIME
    for i in range(min(dim_omega3_tv, 3)):
        old_nonzero = sum(1 for j in old3 if tv_omega_A[i][j] % PRIME != 0)
        print(f"    tv_omega vec {i}: old_nonzero={old_nonzero}")

    bd3_on_tv = D3 @ tv_omega_A.T % PRIME
    rk_full_d3_tv = int(_gauss_rank_np(bd3_on_tv.copy(), PRIME))
    Ktv_B = dim_omega3_tv - rk_full_d3_tv
    print(f"  rk_full_d3_tv={rk_full_d3_tv}, Ktv_B={Ktv_B}")

    # KEY CHECK: are the Ktv_A explicit vectors all in Omega_3^{tv-only}?
    print(f"\n  === RECONCILIATION ===")
    print(f"  Ktv_A={Ktv_A}, Ktv_B={Ktv_B}, diff={Ktv_A - Ktv_B}")

    # Each explicit K_tv vector should be expressible as tv_omega @ ob3 linear combination
    # Check: is Ktv_explicit in the span of tv_omega_A?
    if Ktv_explicit.shape[0] > 0 and tv_omega_A.shape[0] > 0:
        combined = np.vstack([tv_omega_A, Ktv_explicit]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        rk_tv_omega = int(_gauss_rank_np(tv_omega_A.copy(), PRIME))
        new_dirs = rk_combined - rk_tv_omega
        print(f"  K_tv vectors NOT in Omega_3^{'{tv-only}'}: {new_dirs} new directions")

        if new_dirs > 0:
            # Find which K_tv vectors are NOT in tv-only Omega
            for i in range(Ktv_explicit.shape[0]):
                z = Ktv_explicit[i:i+1]
                test = np.vstack([tv_omega_A, z]) % PRIME
                rk_test = int(_gauss_rank_np(test.copy(), PRIME))
                in_span = (rk_test == rk_tv_omega)
                z_old = [z[0, j] % PRIME for j in old3]
                z_old_nonzero = sum(1 for x in z_old if x != 0)
                print(f"    K_tv vec {i}: in tv-only Omega span = {in_span}, "
                      f"old_nonzero = {z_old_nonzero}")

                # This vector has z_old = 0 (it's in K_tv) but is NOT in Omega_3^{tv-only}
                # How? z is in Omega_3 (it's a linear combination of ob3 rows) with z_old = 0
                # So z SHOULD be in Omega_3^{tv-only}

                # Double-check: can we express z as c @ ob3?
                # z = ktv_coeffs[i] @ K = ktv_coeffs[i] @ ker_d3 @ ob3
                # So z = (ktv_coeffs[i] @ ker_d3) @ ob3
                # The coefficient vector is ktv_coeffs[i] @ ker_d3 in Omega coords
                omega_coeffs = ktv_coeffs[i] @ ker_d3 % PRIME
                # Verify: omega_coeffs @ ob3 should equal z
                z_check = omega_coeffs @ ob3 % PRIME
                matches = np.array_equal(z_check % PRIME, z[0] % PRIME)
                print(f"      omega_coeffs valid: {matches}")

                # Now: omega_coeffs @ ob3_old should be zero
                old_check = omega_coeffs @ ob3_old % PRIME
                old_check_nonzero = int(np.count_nonzero(old_check % PRIME))
                print(f"      omega_coeffs @ ob3_old nonzero: {old_check_nonzero}")

                if old_check_nonzero > 0:
                    print(f"      BUG! This K_tv vector has z_old = 0 in A_3 but")
                    print(f"      omega_coeffs @ ob3_old ≠ 0. Must be an arithmetic error.")
                    # Check if it's a modular arithmetic issue
                    z_old_exact = omega_coeffs @ ob3_old
                    z_old_mod = z_old_exact % PRIME
                    print(f"      Raw old values (first 5): {z_old_exact[:5]}")
                    print(f"      Mod old values (first 5): {z_old_mod[:5]}")
                else:
                    # omega_coeffs has ob3_old component = 0
                    # So it SHOULD be in tv_omega span
                    # Check if it's in tv_omega span
                    coeff_test = np.vstack([tv_omega, omega_coeffs.reshape(1, -1)]) % PRIME
                    rk_ct = int(_gauss_rank_np(coeff_test.copy(), PRIME))
                    rk_tv_o = int(_gauss_rank_np(tv_omega.copy(), PRIME))
                    print(f"      In tv_omega span (Omega coords): {rk_ct == rk_tv_o}")
    else:
        print(f"  No K_tv vectors to check")


def main():
    n = 7
    rng = np.random.RandomState(42)

    # Generate tournament at trial=58
    trial_target = 58
    b3_1_count = 0
    for trial in range(50000):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue
        if trial == trial_target:
            print(f"Found target tournament at trial={trial}")
            v = 2
            print(f"Debugging v={v}")
            deep_debug(A, n, v)
            break


if __name__ == '__main__':
    main()
    print("\nDONE.")
