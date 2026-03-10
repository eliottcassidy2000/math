"""
seesaw_mechanism_n7.py — Why does b3=1 force b4=0 at n=7?

The seesaw b3*b4=0 holds at n=7 (exhaustive) but fails at n=8.
What's special about n=7?

Key structural question: When b3(T)=1 at n=7, what prevents rank(d_5) from
being < ker(d_4)? (Since beta_4 = ker(d_4) - rank(d_5))

Approach: For b3=1 tournaments at n=7, analyze:
1. dim(Omega_p) for all p — the chain complex dimensions
2. ker(d_4) and im(d_5) — the exact ranks
3. Whether ker(d_4) is "forced" to be zero by dimensional constraints

If dim(Omega_5) >= ker(d_4) and d_5 has full rank onto ker(d_4),
then b4=0.

Author: opus-2026-03-09-S56
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
    _gauss_rank_np, _gauss_nullbasis_modp,
    _build_constraint_matrix,
    full_chain_complex_modp,
    RANK_PRIME
)


def analyze_seesaw(A, n, verbose=False):
    """Analyze why b4=0 when b3=1."""
    max_p = n - 1
    cc = full_chain_complex_modp(A, n, max_p)
    bettis = {p: cc['bettis'].get(p, 0) for p in range(n)}

    if bettis.get(3, 0) != 1:
        return None

    # Get detailed chain complex info
    dims = {}
    ranks = {}
    for p in range(n):
        dims[p] = cc.get(f'dim_omega_{p}', 0)
        ranks[p] = cc.get(f'rank_d_{p}', 0)

    # Recompute from scratch for full info
    ap = enumerate_all_allowed(A, n, max_p)

    dim_A = {}
    dim_Omega = {}
    for p in range(n):
        paths = ap.get(p, [])
        dim_A[p] = len(paths)
        if not paths:
            dim_Omega[p] = 0
            continue
        P, nr, nc = _build_constraint_matrix(ap, p, RANK_PRIME)
        if P is not None:
            r, nb = _gauss_nullbasis_modp(P, nr, nc, RANK_PRIME)
            dim_Omega[p] = len(nb) if nb else 0
        else:
            dim_Omega[p] = len(paths)

    result = {
        'bettis': bettis,
        'dim_Omega': dim_Omega,
        'dim_A': dim_A,
    }

    # Compute ranks from the chain complex data
    # rank(d_p) = dim(Omega_p) - ker(d_p) = dim(Omega_p) - (beta_p + rank(d_{p+1}))
    # Or directly: rank(d_p) is in cc

    if verbose:
        print(f"  bettis: {[bettis.get(p, 0) for p in range(n)]}")
        print(f"  dim(Omega_p): {[dim_Omega.get(p, 0) for p in range(n)]}")
        print(f"  dim(A_p): {[dim_A.get(p, 0) for p in range(n)]}")

    return result


def main():
    print("=" * 70)
    print("SEESAW MECHANISM AT n=7: WHY b3=1 => b4=0")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)

    found_b3_1 = 0
    found_b3_0 = 0
    target = 200
    t0 = time.time()
    results_b3_1 = []
    results_b3_0 = []

    for trial in range(20000):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        b3 = cc['bettis'].get(3, 0)

        if b3 == 1 and found_b3_1 < target:
            found_b3_1 += 1
            r = analyze_seesaw(A, n, verbose=(found_b3_1 <= 2))
            if r:
                results_b3_1.append(r)
        elif b3 == 0 and found_b3_0 < 50:
            found_b3_0 += 1
            r = analyze_seesaw(A, n)
            # r is None if b3 != 1, compute directly
            ap = enumerate_all_allowed(A, n, n - 1)
            dim_Omega = {}
            for p in range(n):
                paths = ap.get(p, [])
                if not paths:
                    dim_Omega[p] = 0
                    continue
                P, nr, nc = _build_constraint_matrix(ap, p, RANK_PRIME)
                if P is not None:
                    _, nb = _gauss_nullbasis_modp(P, nr, nc, RANK_PRIME)
                    dim_Omega[p] = len(nb) if nb else 0
                else:
                    dim_Omega[p] = len(paths)
            results_b3_0.append({
                'bettis': {p: cc['bettis'].get(p, 0) for p in range(n)},
                'dim_Omega': dim_Omega,
            })

        if found_b3_1 >= target and found_b3_0 >= 50:
            break

    elapsed = time.time() - t0
    print(f"\n  Collected {found_b3_1} b3=1 and {found_b3_0} b3=0, {elapsed:.1f}s")

    # Compare Omega dimensions between b3=1 and b3=0
    print(f"\n  Omega dimensions for b3=1 tournaments:")
    for p in range(n):
        vals = [r['dim_Omega'][p] for r in results_b3_1]
        print(f"    Omega_{p}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

    print(f"\n  Omega dimensions for b3=0 tournaments:")
    for p in range(n):
        vals = [r['dim_Omega'][p] for r in results_b3_0]
        print(f"    Omega_{p}: min={min(vals)}, max={max(vals)}, mean={np.mean(vals):.1f}")

    # Key: Omega_4 and Omega_5 for b3=1
    print(f"\n  KEY: Omega_4 vs Omega_5 for b3=1:")
    for r in results_b3_1[:10]:
        o4 = r['dim_Omega'][4]
        o5 = r['dim_Omega'][5]
        o3 = r['dim_Omega'][3]
        print(f"    Omega_3={o3}, Omega_4={o4}, Omega_5={o5}, ratio_5/4={o5/o4:.2f}" if o4 > 0 else
              f"    Omega_3={o3}, Omega_4={o4}, Omega_5={o5}")

    # Euler characteristic consistency
    print(f"\n  Euler characteristic for b3=1:")
    chi_vals = []
    for r in results_b3_1:
        chi = sum((-1)**p * r['dim_Omega'].get(p, 0) for p in range(n))
        chi_vals.append(chi)
    chi_dist = Counter(chi_vals)
    print(f"    chi distribution: {dict(sorted(chi_dist.items()))}")

    # At n=7 with b3=1: chi = b0 - b1 + b2 - b3 + b4 - b5 + b6
    # = 1 - 0 + 0 - 1 + 0 - 0 + 0 = 0 (since b1=0 for b3=1 by seesaw)
    # If chi = 0 and we know b0=1, b1=0, b2=0, b3=1, b5=0, b6=0:
    # then b4 = 0 (to make sum = 0)!
    print(f"\n  INSIGHT: If chi=0 and all other bettis known,")
    print(f"  then b4 = chi - 1 + 0 - 0 + 1 - 0 + 0 = chi")
    print(f"  So b4 = 0 iff chi = 0 for b3=1 tournaments!")


if __name__ == '__main__':
    main()
    print("\nDONE.")
