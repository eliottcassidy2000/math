"""
h3_generator_structure.py — Analyze the H_3 generator structure for HYP-408

KEY THEOREM (kind-pasteur-S50):
  For tournaments with beta_3(T) = 1:
    rank(i_*: H_3(T\v) -> H_3(T)) = beta_3(T\v)

  This gives a CLEAN proof of HYP-408 for Case 1 (b3(T\v) = 1):
    When b3(T\v) = 1, rank(i_*) = 1, so H_3(T) is spanned by the embedded
    H_3(T\v) generator, which is purely old. Therefore K_tv maps to 0 in H_3(T).

  The HARD case is Case 2: b3(T\v) = 0.
    Here rank(i_*) = 0 and the H_3 generator is NOT an embedded T\v cycle.
    Yet codim = 1 still holds (K_tv = B_tv).
    WHY does the generator have nontrivial old component?

THIS SCRIPT: Extracts the H_3(T) generator for Case 2, analyzes:
1. The old/tv support decomposition of the generator
2. Whether the old component is a BOUNDARY in T\v
3. Whether the tv component alone determines the class
4. Whether ALL vertices v with b3(T\v) = 0 give the same generator structure

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
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


def generator_analysis(A, n, v):
    """Extract and analyze the H_3(T) generator structure."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1

    cc = full_chain_complex_modp(A, n, n - 1)
    if cc['bettis'].get(3, 0) != 1:
        return None

    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, n1 - 1)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

    # Only interested in Case 2: b3(T\v) = 0
    if b3_Tv != 0:
        return {'case': 1, 'b3_Tv': b3_Tv}

    ap = enumerate_all_allowed(A, n, min(n-1, 6))
    paths = {deg: ap.get(deg, []) for deg in range(7)}
    old_idx3 = [i for i, p in enumerate(paths[3]) if v not in p]
    tv_idx3 = [i for i, p in enumerate(paths[3]) if v in p]

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

    # Compute ker(d_3) and im(d_4) in A_3 coords
    d3o = D3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        K_coeffs = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        K_coeffs = np.eye(ob3.shape[0], dtype=np.int64)

    K = K_coeffs @ ob3 % PRIME

    im_d4 = (D4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    dim_K = K.shape[0]
    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    # Find the H_3 generator: a vector in K not in im(d_4)
    # Stack im_d4 and K, find the new direction
    if im_d4.shape[0] > 0:
        combined = np.vstack([im_d4, K]) % PRIME
    else:
        combined = K.copy()
    # Row-reduce to find the generator (last vector not in im_d4 span)
    # Actually, just use the complementary basis approach
    gen_found = False
    generator = None
    for i in range(K.shape[0]):
        z = K[i:i+1]
        if im_d4.shape[0] > 0:
            test = np.vstack([im_d4, z]) % PRIME
            rk = int(_gauss_rank_np(test.copy(), PRIME))
            if rk > rk_d4:
                generator = z[0]
                gen_found = True
                break
        else:
            generator = z[0]
            gen_found = True
            break

    if not gen_found:
        return None

    # Analyze generator structure
    old_support = sum(1 for j in old_idx3 if generator[j] % PRIME != 0)
    tv_support = sum(1 for j in tv_idx3 if generator[j] % PRIME != 0)
    total_support = old_support + tv_support

    # Is old component zero?
    old_is_zero = (old_support == 0)

    # What fraction is old vs tv?
    old_fraction = old_support / max(total_support, 1)

    # Is old component a cycle in T\v?
    z_old_in_Tv = np.zeros(len(paths[3]), dtype=np.int64)
    for j in old_idx3:
        z_old_in_Tv[j] = generator[j]
    d3_z_old = D3 @ z_old_in_Tv % PRIME
    old_2_indices = [i for i, p in enumerate(paths[2]) if v not in p]
    tv_2_indices = [i for i, p in enumerate(paths[2]) if v in p]

    old_component_is_cycle = all(d3_z_old[j] == 0 for j in old_2_indices)
    # The old 2-faces of z_old should vanish (since z_old only has old 3-paths)
    # D3_to = 0, so old 3-paths only produce old 2-faces.
    # So d3(z_old) has only old 2-face components.
    old_bd_in_old2 = np.array([d3_z_old[j] for j in old_2_indices])
    old_bd_nonzero = int(np.count_nonzero(old_bd_in_old2 % PRIME))

    # Check tv component: D_tt(z_tv) and D_ot(z_tv)
    z_tv = np.zeros(len(paths[3]), dtype=np.int64)
    for j in tv_idx3:
        z_tv[j] = generator[j]
    d3_z_tv = D3 @ z_tv % PRIME
    tv_bd_in_old2 = sum(1 for j in old_2_indices if d3_z_tv[j] % PRIME != 0)
    tv_bd_in_tv2 = sum(1 for j in tv_2_indices if d3_z_tv[j] % PRIME != 0)

    # K_tv dimension
    old3 = sorted(old_idx3)
    pi_K = K[:, old3] % PRIME
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME)) if pi_K.size > 0 else 0
    dim_Ktv = dim_K - rk_pi_K

    # Position of v in support paths
    v_positions = Counter()
    for j in tv_idx3:
        if generator[j] % PRIME != 0:
            path = paths[3][j]
            v_pos = path.index(v)
            v_positions[v_pos] += 1

    return {
        'case': 2,
        'b3_Tv': 0,
        'old_support': old_support,
        'tv_support': tv_support,
        'total_support': total_support,
        'old_is_zero': old_is_zero,
        'old_fraction': old_fraction,
        'old_component_is_cycle_in_Tv': old_bd_nonzero == 0,
        'old_bd_nonzero': old_bd_nonzero,
        'tv_bd_in_old2': tv_bd_in_old2,
        'tv_bd_in_tv2': tv_bd_in_tv2,
        'dim_Ktv': dim_Ktv,
        'v_positions': dict(v_positions),
        'out_deg': int(sum(A[v])),
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"H_3 GENERATOR STRUCTURE ANALYSIS AT n={n} (Case 2: b3(T\\v)=0)")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        case1_count = 0
        case2_results = []
        t0 = time.time()
        target = 300 if n == 7 else 100

        for trial in range(50000):
            if len(case2_results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue
            for v_cand in range(n):
                r = generator_analysis(A, n, v_cand)
                if r is None:
                    continue
                if r['case'] == 1:
                    case1_count += 1
                else:
                    case2_results.append(r)
                    if len(case2_results) % 50 == 0:
                        print(f"  Progress: {len(case2_results)} Case 2 pairs, {time.time()-t0:.1f}s")

        elapsed = time.time() - t0
        print(f"\n  Total: {len(case2_results)} Case 2, {case1_count} Case 1, {elapsed:.1f}s")

        if not case2_results:
            print("  No Case 2 found!")
            continue

        # KEY: Is old component ever zero?
        old_zero = sum(1 for r in case2_results if r['old_is_zero'])
        print(f"\n  OLD COMPONENT IS ZERO: {old_zero}/{len(case2_results)} "
              f"({'NEVER' if old_zero == 0 else f'{old_zero} CASES!'})")

        # Old fraction statistics
        fracs = [r['old_fraction'] for r in case2_results]
        print(f"\n  OLD FRACTION: mean={np.mean(fracs):.3f}, "
              f"min={min(fracs):.3f}, max={max(fracs):.3f}")

        # Support sizes
        print(f"\n  SUPPORT SIZES:")
        print(f"    old: mean={np.mean([r['old_support'] for r in case2_results]):.1f}")
        print(f"    tv: mean={np.mean([r['tv_support'] for r in case2_results]):.1f}")

        # Is old component a cycle in T\v?
        old_cycle = sum(1 for r in case2_results if r['old_component_is_cycle_in_Tv'])
        print(f"\n  OLD COMPONENT IS CYCLE IN T\\v: {old_cycle}/{len(case2_results)}")
        if old_cycle > 0:
            # This would be surprising — since b3(T\v) = 0, any cycle in T\v is a boundary
            print(f"    (Since b3(T\\v)=0, these old-component cycles are BOUNDARIES in T\\v)")

        # Old boundary nonzero count
        old_bd_vals = [r['old_bd_nonzero'] for r in case2_results]
        print(f"    Old component boundary nonzero entries: mean={np.mean(old_bd_vals):.1f}")

        # TV boundary spillover
        tv_old_spill = [r['tv_bd_in_old2'] for r in case2_results]
        tv_tv_bd = [r['tv_bd_in_tv2'] for r in case2_results]
        print(f"\n  TV COMPONENT BOUNDARY:")
        print(f"    spillover to old 2-faces: mean={np.mean(tv_old_spill):.1f}")
        print(f"    tv 2-face boundary: mean={np.mean(tv_tv_bd):.1f}")
        # Key check: tv spillover should equal negative of old boundary
        # (since d_3(z) = d_3(z_old) + d_3(z_tv) = 0, so d_3(z_old)_old = -d_3(z_tv)_old)

        # V position distribution
        all_positions = Counter()
        for r in case2_results:
            for pos, cnt in r['v_positions'].items():
                all_positions[pos] += cnt
        print(f"\n  V POSITION IN SUPPORT PATHS: {dict(sorted(all_positions.items()))}")

        # Out-degree distribution
        out_degs = [r['out_deg'] for r in case2_results]
        print(f"  OUT-DEGREE: dist={dict(Counter(out_degs))}")

        # K_tv dimension
        ktv_vals = [r['dim_Ktv'] for r in case2_results]
        print(f"  K_tv dimension: dist={dict(Counter(ktv_vals))}")

        # Cross: K_tv > 0 with generator analysis
        ktv_pos = [r for r in case2_results if r['dim_Ktv'] > 0]
        print(f"\n  CASES WITH K_tv > 0: {len(ktv_pos)}/{len(case2_results)}")
        if ktv_pos:
            print(f"    Old fraction: mean={np.mean([r['old_fraction'] for r in ktv_pos]):.3f}")
            print(f"    Old support: mean={np.mean([r['old_support'] for r in ktv_pos]):.1f}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
