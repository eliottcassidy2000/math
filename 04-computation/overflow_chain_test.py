"""
overflow_chain_test.py — Test chained matrix multiply overflow

The suspected bug: D3 @ (tv_omega @ ob3).T % PRIME
- tv_omega @ ob3 doesn't overflow int64 (entries < 2^63)
- BUT entries can be >> PRIME (they're correct but unreduced)
- D3 entries are up to PRIME-1 ≈ 2^31
- D3 @ X.T: products D3[i,k] * X[k,j] up to 2^31 * 2^62 = 2^93 > 2^63 → OVERFLOW!

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


def check_chain_overflow(A, n, v):
    """Check for overflow in D3 @ (tv_omega @ ob3).T."""
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

    d3o = D3 @ ob3.T % PRIME
    _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
    ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    K = ker_d3 @ ob3 % PRIME
    dim_K = K.shape[0]

    im_d4 = (D4 @ ob4.T % PRIME).T
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME))
    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    ob3_old = ob3[:, old3] % PRIME
    rk_ob3_old = int(_gauss_rank_np(ob3_old.copy(), PRIME))
    ob3_old_T = ob3_old.T % PRIME
    _, tv_omega_basis = _gauss_nullbasis_modp(
        ob3_old_T.astype(np.int32), ob3_old_T.shape[0], ob3_old_T.shape[1], PRIME)
    tv_omega = np.array(tv_omega_basis, dtype=np.int64) if tv_omega_basis else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    dim_omega3_tv = tv_omega.shape[0]

    if dim_omega3_tv == 0:
        return None

    # The unreduced intermediate
    X = tv_omega @ ob3  # NOT reduced mod PRIME!
    max_X = int(np.max(np.abs(X)))

    # Max entry of D3
    max_D3 = int(np.max(D3))

    # Potential product overflow
    # D3 max entry * X max entry
    max_product = max_D3 * max_X

    # Does this overflow int64?
    overflow = max_product > (2**63 - 1)

    # Compute both ways
    bd3_chained = D3 @ X.T % PRIME  # potentially overflowed
    X_reduced = X % PRIME
    bd3_safe = D3 @ X_reduced.T % PRIME

    mismatch = not np.array_equal(bd3_chained, bd3_safe)

    rk_chained = int(_gauss_rank_np(bd3_chained.copy(), PRIME))
    rk_safe = int(_gauss_rank_np(bd3_safe.copy(), PRIME))

    Ktv_chained = dim_omega3_tv - rk_chained
    Ktv_safe = dim_omega3_tv - rk_safe

    return {
        'max_X': max_X,
        'max_D3': max_D3,
        'max_product': max_product,
        'overflow': overflow,
        'mismatch': mismatch,
        'Ktv_chained': Ktv_chained,
        'Ktv_safe': Ktv_safe,
        'dim_omega3_tv': dim_omega3_tv,
    }


def main():
    n = 7
    print(f"{'='*70}")
    print(f"CHAIN OVERFLOW TEST AT n={n}")
    print(f"{'='*70}")

    rng = np.random.RandomState(42)
    results = []
    overflow_count = 0
    mismatch_count = 0
    ktv_mismatch_count = 0

    for trial in range(50000):
        if len(results) >= 600:
            break
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue

        for v_cand in range(n):
            r = check_chain_overflow(A, n, v_cand)
            if r is None:
                continue
            results.append(r)
            if r['overflow']:
                overflow_count += 1
            if r['mismatch']:
                mismatch_count += 1
            if r['Ktv_chained'] != r['Ktv_safe']:
                ktv_mismatch_count += 1
                print(f"  trial={trial}, v={v_cand}: Ktv chained={r['Ktv_chained']} safe={r['Ktv_safe']} "
                      f"max_X={r['max_X']:.2e} max_prod={r['max_product']:.2e} overflow={r['overflow']}")

    print(f"\n  Total: {len(results)} pairs")
    print(f"  D3 @ X.T overflow possible: {overflow_count}/{len(results)}")
    print(f"  D3 @ X.T actual mismatch: {mismatch_count}/{len(results)}")
    print(f"  Ktv mismatch: {ktv_mismatch_count}/{len(results)}")

    # Max values
    max_max_X = max(r['max_X'] for r in results)
    max_max_product = max(r['max_product'] for r in results)
    print(f"  Max X entry: {max_max_X:.4e} (int64 max: {2**63-1:.4e})")
    print(f"  Max D3*X product: {max_max_product:.4e}")
    print(f"  Overflow threshold: {2**63-1:.4e}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
