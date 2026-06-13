"""
ghost_cycle_debug.py — Debug the K_tv/B_tv discrepancy between two methods

Method A (kind-pasteur): K_tv = dim(K) - rk(K[:, old3])
Method B (opus): K_tv = dim(Omega_3^{tv-only}) - rk(d_3 on Omega_3^{tv-only})

Both SHOULD be equivalent. This script computes both and finds where they disagree.

Author: kind-pasteur-2026-03-10-S50
"""
import sys
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


def compare_methods(A, n, v):
    """Compare two methods for computing K_tv and B_tv."""
    max_p = min(n - 1, 6)
    ap = enumerate_all_allowed(A, n, max_p)

    paths = {deg: ap.get(deg, []) for deg in range(7)}

    old3 = sorted([i for i, p in enumerate(paths[3]) if v not in p])
    tv3 = sorted([i for i, p in enumerate(paths[3]) if v in p])
    old4 = sorted([i for i, p in enumerate(paths[4]) if v not in p])
    tv4 = sorted([i for i, p in enumerate(paths[4]) if v in p])

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

    # === METHOD A (kind-pasteur): K_tv via ker(d_3) projected ===
    d3o = D3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    K = ker_d3 @ ob3 % PRIME
    dim_K = K.shape[0]

    im_d4 = (D4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    # K_tv method A
    pi_K = K[:, old3] % PRIME
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME)) if pi_K.size > 0 else 0
    Ktv_A = dim_K - rk_pi_K

    # B_tv method A
    pi_B = im_d4[:, old3] % PRIME if im_d4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64)
    rk_pi_B = int(_gauss_rank_np(pi_B.copy(), PRIME)) if pi_B.shape[0] > 0 else 0
    Btv_A = rk_d4 - rk_pi_B

    # === METHOD B (opus): K_tv via tv-only Omega ===
    ob3_old = ob3[:, old3] % PRIME
    ob3_tv = ob3[:, tv3] % PRIME

    rk_ob3_old = int(_gauss_rank_np(ob3_old.copy(), PRIME)) if ob3_old.shape[1] > 0 else 0
    if rk_ob3_old > 0:
        ob3_old_T = ob3_old.T % PRIME
        _, tv_omega_basis = _gauss_nullbasis_modp(
            ob3_old_T.astype(np.int32), ob3_old_T.shape[0], ob3_old_T.shape[1], PRIME)
        tv_omega = np.array(tv_omega_basis, dtype=np.int64) if tv_omega_basis else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        tv_omega = np.eye(ob3.shape[0], dtype=np.int64)

    dim_omega3_tv = tv_omega.shape[0]

    # d_3 on tv-only Omega vectors
    if dim_omega3_tv > 0:
        bd3_on_tv = D3 @ (tv_omega @ ob3).T % PRIME
        rk_full_d3_tv = int(_gauss_rank_np(bd3_on_tv.copy(), PRIME))
    else:
        rk_full_d3_tv = 0

    Ktv_B = dim_omega3_tv - rk_full_d3_tv

    # B_tv method B (same as A)
    Btv_B = Btv_A  # Same computation

    return {
        'Ktv_A': Ktv_A, 'Btv_A': Btv_A,
        'Ktv_B': Ktv_B, 'Btv_B': Btv_B,
        'agree_Ktv': Ktv_A == Ktv_B,
        'agree_Btv': Btv_A == Btv_B,
        'gc_A': Ktv_A == Btv_A,
        'gc_B': Ktv_B == Btv_B,
        'dim_K': dim_K,
        'rk_d4': rk_d4,
        'rk_pi_K': rk_pi_K,
        'rk_pi_B': rk_pi_B,
        'dim_omega3_tv': dim_omega3_tv,
        'rk_full_d3_tv': rk_full_d3_tv,
    }


def main():
    n = 7
    print(f"{'='*70}")
    print(f"K_tv/B_tv METHOD COMPARISON AT n={n}")
    print(f"{'='*70}")

    rng = np.random.RandomState(42)
    results = []
    disagree = []

    for trial in range(50000):
        if len(results) >= 600:
            break
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        if cc['bettis'].get(3, 0) != 1:
            continue

        for v_cand in range(n):
            r = compare_methods(A, n, v_cand)
            if r is None:
                continue
            results.append(r)
            if not r['agree_Ktv']:
                disagree.append((trial, v_cand, r))
                print(f"  DISAGREE at trial={trial}, v={v_cand}: "
                      f"Ktv_A={r['Ktv_A']}, Ktv_B={r['Ktv_B']}")

    print(f"\n  Total: {len(results)} pairs")
    print(f"  K_tv agree: {sum(1 for r in results if r['agree_Ktv'])}/{len(results)}")
    print(f"  Ghost Cycle (method A): {sum(1 for r in results if r['gc_A'])}/{len(results)}")
    print(f"  Ghost Cycle (method B): {sum(1 for r in results if r['gc_B'])}/{len(results)}")

    if disagree:
        print(f"\n  DISAGREEMENTS ({len(disagree)}):")
        for trial, v, r in disagree[:5]:
            print(f"    trial={trial}, v={v}")
            print(f"      Ktv_A={r['Ktv_A']} (dim_K={r['dim_K']}, rk_pi_K={r['rk_pi_K']})")
            print(f"      Ktv_B={r['Ktv_B']} (dim_omega3_tv={r['dim_omega3_tv']}, rk_d3_tv={r['rk_full_d3_tv']})")
    else:
        print(f"\n  All K_tv computations AGREE between methods.")
        print(f"  The 'failures' in opus's output must come from B_tv computation difference.")

        # Check B_tv marginals
        Ktv_dist = Counter(r['Ktv_A'] for r in results)
        Btv_dist = Counter(r['Btv_A'] for r in results)
        print(f"\n  K_tv dist: {dict(sorted(Ktv_dist.items()))}")
        print(f"  B_tv dist: {dict(sorted(Btv_dist.items()))}")
        print(f"  K_tv - B_tv: {dict(sorted(Counter(r['Ktv_A'] - r['Btv_A'] for r in results).items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
