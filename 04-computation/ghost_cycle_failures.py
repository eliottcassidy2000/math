"""
ghost_cycle_failures.py — Find and analyze cases where K_tv ≠ B_tv

From opus-S59's tv_cycle_structure.py: Ghost Cycle (K_tv = B_tv) fails in 14/504 at n=7.
By the Ghost Cycle ⟺ HYP-408 equivalence (opus-S58):
  K_tv ≠ B_tv ⟺ codim(π_old(B), π_old(K)) ≠ 1

This script:
1. Reproduces the opus finding
2. For each failure case, computes codim directly
3. Checks if beta_3(T) is still 1 (sanity)
4. Checks what beta_3(T\v) is for the failing vertex
5. Characterizes the failing tournaments

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


def full_analysis(A, n, v):
    """Compute K_tv, B_tv, and codim for a (T,v) pair."""
    remaining = [i for i in range(n) if i != v]
    n1 = n - 1

    # Beta_3(T)
    cc = full_chain_complex_modp(A, n, n - 1)
    b3_T = cc['bettis'].get(3, 0)
    if b3_T != 1:
        return None

    # Beta_3(T\v)
    A_sub = [[A[remaining[i]][remaining[j]] for j in range(n1)] for i in range(n1)]
    cc_Tv = full_chain_complex_modp(A_sub, n1, n1 - 1)
    b3_Tv = cc_Tv['bettis'].get(3, 0)

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

    # ker(d_3) in Omega coords
    d3o = D3 @ ob3.T % PRIME
    if d3o.size > 0:
        _, kv = _gauss_nullbasis_modp(d3o.astype(np.int32), d3o.shape[0], d3o.shape[1], PRIME)
        ker_d3 = np.array(kv, dtype=np.int64) if kv else np.zeros((0, ob3.shape[0]), dtype=np.int64)
    else:
        ker_d3 = np.eye(ob3.shape[0], dtype=np.int64)

    K = ker_d3 @ ob3 % PRIME
    dim_K = K.shape[0]

    # im(d_4) in A_3 coords
    im_d4 = (D4 @ ob4.T % PRIME).T if ob4.shape[0] > 0 else np.zeros((0, len(paths[3])), dtype=np.int64)
    rk_d4 = int(_gauss_rank_np(im_d4.copy(), PRIME)) if im_d4.shape[0] > 0 else 0

    beta3 = dim_K - rk_d4
    if beta3 != 1:
        return None

    # K_tv = ker(π_old|_K)
    pi_K = K[:, old3] % PRIME
    rk_pi_K = int(_gauss_rank_np(pi_K.copy(), PRIME)) if pi_K.size > 0 else 0
    dim_Ktv = dim_K - rk_pi_K

    # B_tv = ker(π_old|_B)
    if im_d4.shape[0] > 0:
        pi_B = im_d4[:, old3] % PRIME
        rk_pi_B = int(_gauss_rank_np(pi_B.copy(), PRIME))
        dim_Btv = rk_d4 - rk_pi_B
    else:
        dim_Btv = 0

    # codim = dim(π_old(K)) - dim(π_old(B))
    # But we should compute codim properly: dim(π_old(K)/π_old(B))
    # = rank of π_old(K) relative to π_old(B)
    if rk_pi_B > 0:
        combined = np.vstack([pi_B if im_d4.shape[0] > 0 else np.zeros((0, len(old3)), dtype=np.int64), pi_K]) % PRIME
        rk_combined = int(_gauss_rank_np(combined.copy(), PRIME))
        codim = rk_combined - rk_pi_B
    else:
        codim = rk_pi_K

    # Actually: codim = rk(π_old(K)) - rk(π_old(B)) when π_old(B) ⊂ π_old(K)
    # Since B ⊂ K, we have π_old(B) ⊂ π_old(K), so this is correct.
    # But need to check: rk_pi_K - rk_pi_B or the combined rank approach.
    # Since π_old(B) ⊂ π_old(K), rk(combined) = rk_pi_K, so codim = rk_pi_K - rk_pi_B.
    codim_direct = rk_pi_K - rk_pi_B

    # Verify: codim = beta_3 - (dim_Ktv - dim_Btv)
    codim_from_formula = beta3 - (dim_Ktv - dim_Btv)

    out_deg = int(sum(A[v]))

    return {
        'b3_T': b3_T,
        'b3_Tv': b3_Tv,
        'dim_K': dim_K,
        'rk_d4': rk_d4,
        'beta3': beta3,
        'dim_Ktv': dim_Ktv,
        'dim_Btv': dim_Btv,
        'codim_direct': codim_direct,
        'codim_from_formula': codim_from_formula,
        'ghost_cycle': dim_Ktv == dim_Btv,
        'rk_pi_K': rk_pi_K,
        'rk_pi_B': rk_pi_B,
        'out_deg': out_deg,
    }


def main():
    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"GHOST CYCLE FAILURE ANALYSIS AT n={n}")
        print(f"{'='*70}")

        rng = np.random.RandomState(42)
        results = []
        failures = []
        t0 = time.time()
        target = 600 if n == 7 else 400

        for trial in range(80000):
            if len(results) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            if cc['bettis'].get(3, 0) != 1:
                continue

            for v_cand in range(n):
                r = full_analysis(A, n, v_cand)
                if r is None:
                    continue
                results.append(r)
                if not r['ghost_cycle']:
                    failures.append((trial, v_cand, r, A.copy()))
                    print(f"  FAILURE at trial={trial}, v={v_cand}: "
                          f"K_tv={r['dim_Ktv']}, B_tv={r['dim_Btv']}, "
                          f"codim={r['codim_direct']}, b3_Tv={r['b3_Tv']}")

            if len(results) % 100 == 0 and len(results) > 0:
                elapsed = time.time() - t0
                print(f"  Progress: {len(results)} pairs, {len(failures)} failures, {elapsed:.1f}s")

        elapsed = time.time() - t0
        gc_count = sum(1 for r in results if r['ghost_cycle'])
        print(f"\n  Total: {len(results)} pairs, Ghost Cycle: {gc_count}/{len(results)}, {elapsed:.1f}s")

        # Summary
        print(f"\n  FAILURES: {len(failures)}/{len(results)}")

        if failures:
            print(f"\n  FAILURE DETAILS:")
            for trial, v, r, A_fail in failures[:10]:
                print(f"    trial={trial}, v={v}: K_tv={r['dim_Ktv']}, B_tv={r['dim_Btv']}, "
                      f"codim={r['codim_direct']}, b3_Tv={r['b3_Tv']}, out_deg={r['out_deg']}")
                print(f"      codim_formula={r['codim_from_formula']}")

            # b3_Tv distribution for failures
            b3_fail = Counter(r['b3_Tv'] for _, _, r, _ in failures)
            print(f"\n  Failure b3(T\\v) dist: {dict(b3_fail)}")

            # b3_Tv distribution for successes
            b3_pass = Counter(r['b3_Tv'] for r in results if r['ghost_cycle'])
            print(f"  Success b3(T\\v) dist: {dict(b3_pass)}")

            # codim distribution
            codim_dist = Counter(r['codim_direct'] for r in results)
            print(f"\n  Codim distribution (ALL): {dict(sorted(codim_dist.items()))}")
            codim_fail = Counter(r['codim_direct'] for _, _, r, _ in failures)
            print(f"  Codim distribution (failures): {dict(sorted(codim_fail.items()))}")

            # dim_Ktv - dim_Btv
            diff_dist = Counter(r['dim_Ktv'] - r['dim_Btv'] for r in results)
            print(f"\n  K_tv - B_tv distribution: {dict(sorted(diff_dist.items()))}")

        # Cross-tabulation: codim vs b3_Tv
        print(f"\n  CROSS-TAB: codim x b3_Tv")
        for b3v in sorted(set(r['b3_Tv'] for r in results)):
            sub = [r for r in results if r['b3_Tv'] == b3v]
            codim_dist = Counter(r['codim_direct'] for r in sub)
            gc_sub = sum(1 for r in sub if r['ghost_cycle'])
            print(f"    b3_Tv={b3v}: codim={dict(sorted(codim_dist.items()))}, "
                  f"ghost_cycle={gc_sub}/{len(sub)}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
