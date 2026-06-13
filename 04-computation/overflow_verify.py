"""
overflow_verify.py — Verify that int64 overflow affects K_tv computation

RANK_PRIME = 2^31 - 1. Matrix entries up to PRIME-1 ≈ 2^31.
Product of two entries: up to (2^31)^2 = 2^62.
Sum of 50 products: up to 50 * 2^62 ≈ 2^67.6 > 2^63 = int64 max.

This means numpy's @ operator gives WRONG results for
    tv_omega @ ob3
when the inner dimension is large enough.

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


def safe_matmul_mod(A, B, prime):
    """Safe matrix multiplication A @ B mod prime, avoiding int64 overflow.

    With prime ≈ 2^31, each product ≈ 2^62. Safe to sum at most
    floor(2^63 / 2^62) = 2 products before taking mod.
    We chunk the inner dimension accordingly.
    """
    n, k = A.shape
    _, m = B.shape
    # Max inner-dim elements per chunk before overflow risk
    max_per_chunk = max(1, (2**63 - 1) // ((prime - 1) * (prime - 1)))
    C = np.zeros((n, m), dtype=np.int64)
    for start in range(0, k, max_per_chunk):
        end = min(start + max_per_chunk, k)
        C = (C + A[:, start:end] @ B[start:end, :]) % prime
    return C


def compare_with_safe(A, n, v):
    """Compare unsafe vs safe matrix multiplication for K_tv."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths = {deg: ap.get(deg, []) for deg in range(7)}

    old3 = sorted([i for i, p in enumerate(paths[3]) if v not in p])
    tv3 = sorted([i for i, p in enumerate(paths[3]) if v in p])

    if not old3 or not tv3 or not paths[4]:
        return None

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

    # ---- METHOD A (my standard) ----
    d3o_unsafe = D3 @ ob3.T % PRIME
    d3o_safe = safe_matmul_mod(D3, ob3.T, PRIME)
    d3o_match = np.array_equal(d3o_unsafe, d3o_safe)

    if not d3o_match:
        # Count mismatches
        mismatches = np.sum(d3o_unsafe != d3o_safe)
        print(f"  d3o MISMATCH! {mismatches} entries differ")

    _, kv_safe = _gauss_nullbasis_modp(d3o_safe.astype(np.int32), d3o_safe.shape[0], d3o_safe.shape[1], PRIME)
    ker_d3_safe = np.array(kv_safe, dtype=np.int64) if kv_safe else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    K_safe = safe_matmul_mod(ker_d3_safe, ob3, PRIME)
    dim_K = K_safe.shape[0]

    _, kv_unsafe = _gauss_nullbasis_modp(d3o_unsafe.astype(np.int32), d3o_unsafe.shape[0], d3o_unsafe.shape[1], PRIME)
    ker_d3_unsafe = np.array(kv_unsafe, dtype=np.int64) if kv_unsafe else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    K_unsafe = ker_d3_unsafe @ ob3 % PRIME

    im_d4_safe = safe_matmul_mod(D4, ob4.T, PRIME).T
    rk_d4 = int(_gauss_rank_np(im_d4_safe.copy(), PRIME))
    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    # K_tv via safe
    pi_K_safe = K_safe[:, old3] % PRIME
    rk_pi_K_safe = int(_gauss_rank_np(pi_K_safe.copy(), PRIME))
    Ktv_safe = dim_K - rk_pi_K_safe

    pi_K_unsafe = K_unsafe[:, old3] % PRIME
    rk_pi_K_unsafe = int(_gauss_rank_np(pi_K_unsafe.copy(), PRIME))
    Ktv_unsafe = K_unsafe.shape[0] - rk_pi_K_unsafe

    # B_tv via safe
    pi_B_safe = im_d4_safe[:, old3] % PRIME
    rk_pi_B_safe = int(_gauss_rank_np(pi_B_safe.copy(), PRIME))
    Btv_safe = rk_d4 - rk_pi_B_safe

    # ---- METHOD B (opus's, both safe and unsafe) ----
    ob3_old = ob3[:, old3] % PRIME
    ob3_old_T = ob3_old.T % PRIME
    _, tv_omega_basis = _gauss_nullbasis_modp(
        ob3_old_T.astype(np.int32), ob3_old_T.shape[0], ob3_old_T.shape[1], PRIME)
    tv_omega = np.array(tv_omega_basis, dtype=np.int64) if tv_omega_basis else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    dim_omega3_tv = tv_omega.shape[0]

    # Unsafe
    tv_omega_A_unsafe = tv_omega @ ob3 % PRIME
    bd3_on_tv_unsafe = D3 @ tv_omega_A_unsafe.T % PRIME
    rk_unsafe = int(_gauss_rank_np(bd3_on_tv_unsafe.copy(), PRIME))
    Ktv_B_unsafe = dim_omega3_tv - rk_unsafe

    # Safe
    tv_omega_A_safe = safe_matmul_mod(tv_omega, ob3, PRIME)
    bd3_on_tv_safe = safe_matmul_mod(D3, tv_omega_A_safe.T, PRIME)
    rk_safe_B = int(_gauss_rank_np(bd3_on_tv_safe.copy(), PRIME))
    Ktv_B_safe = dim_omega3_tv - rk_safe_B

    # Check overflow detection
    tv_omega_A_overflow = not np.array_equal(tv_omega_A_unsafe, tv_omega_A_safe)

    return {
        'd3o_match': d3o_match,
        'Ktv_A_safe': Ktv_safe,
        'Ktv_A_unsafe': Ktv_unsafe,
        'Ktv_B_safe': Ktv_B_safe,
        'Ktv_B_unsafe': Ktv_B_unsafe,
        'Btv_safe': Btv_safe,
        'gc_safe': Ktv_safe == Btv_safe,
        'gc_unsafe_A': Ktv_unsafe == Btv_safe,
        'gc_unsafe_B': Ktv_B_unsafe == Btv_safe,
        'tv_omega_overflow': tv_omega_A_overflow,
        'dim_omega3_tv': dim_omega3_tv,
        'inner_dim': ob3.shape[0],
    }


def main():
    for n in [7]:
        print(f"{'='*70}")
        print(f"OVERFLOW VERIFICATION AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        overflow_count = 0
        ktv_mismatch = 0

        for trial in range(50000):
            if len(results) >= 600:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = compare_with_safe(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)

                if r['tv_omega_overflow']:
                    overflow_count += 1
                if r['Ktv_B_safe'] != r['Ktv_B_unsafe']:
                    ktv_mismatch += 1
                    print(f"  trial={trial}, v={v_cand}: "
                          f"Ktv_B safe={r['Ktv_B_safe']} unsafe={r['Ktv_B_unsafe']} "
                          f"overflow={r['tv_omega_overflow']} inner_dim={r['inner_dim']}")

        print(f"\n  Total: {len(results)} pairs")
        print(f"  tv_omega @ ob3 OVERFLOW: {overflow_count}/{len(results)}")
        print(f"  Ktv_B mismatch (safe vs unsafe): {ktv_mismatch}/{len(results)}")

        # Ghost Cycle check
        gc_safe = sum(1 for r in results if r['gc_safe'])
        gc_unsafe_B = sum(1 for r in results if r['gc_unsafe_B'])
        print(f"\n  Ghost Cycle (SAFE): {gc_safe}/{len(results)}")
        print(f"  Ghost Cycle (UNSAFE method B): {gc_unsafe_B}/{len(results)}")

        # d3o match
        d3o_fail = sum(1 for r in results if not r['d3o_match'])
        print(f"  d3o overflow: {d3o_fail}/{len(results)}")

        # Ktv_A safe vs unsafe
        ktv_A_mismatch = sum(1 for r in results if r['Ktv_A_safe'] != r['Ktv_A_unsafe'])
        print(f"  Ktv_A mismatch (safe vs unsafe): {ktv_A_mismatch}/{len(results)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
