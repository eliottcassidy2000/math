"""
beta_n2_investigation.py — Investigate beta_{n-2} being generically nonzero at n=8

Opus-S59 found: beta_6 (= beta_{n-2}) is nonzero for 449/500 = 89.8% at n=8.
For n <= 7, beta_{n-2} = 0 always.

Questions:
1. What is the structural mechanism? beta_{n-2} involves (n-1)-vertex paths = vertex-deletion paths
2. Is beta_{n-2} related to the Hamiltonian path structure?
3. Does the onset shift: beta_{n-2} appears at n=8 for first time?
4. What's dim(Omega_{n-2}) and dim(Omega_{n-1}) at n=8?

The chain complex at top degrees:
  Omega_{n-1} --d_{n-1}--> Omega_{n-2} --d_{n-2}--> Omega_{n-3}
  beta_{n-2} = dim(ker d_{n-2}) - dim(im d_{n-1})

Omega_{n-1} = Hamiltonian paths (using all n vertices)
Omega_{n-2} = chains on (n-1)-vertex paths with all faces allowed

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp, RANK_PRIME,
    enumerate_all_allowed, _build_constraint_matrix,
    _gauss_rank_np, _gauss_nullbasis_modp
)


def investigate_top_degrees(n, num_samples=200):
    """Investigate beta_{n-2} and related dimensions."""
    print(f"\n{'='*60}")
    print(f"BETA_{{n-2}} = BETA_{n-2} AT n={n}")
    print(f"{'='*60}")

    results = []
    t0 = time.time()

    for trial in range(num_samples):
        A = random_tournament(n)
        # Compute up to degree n-1 (full complex)
        max_deg = n - 1
        cc = full_chain_complex_modp(A, n, max_deg)
        bettis = cc['bettis']
        dims = cc.get('omega_dims', {})

        b_nm2 = bettis.get(n - 2, 0)
        b_nm1 = bettis.get(n - 1, 0)
        b_nm3 = bettis.get(n - 3, 0)

        dim_nm1 = dims.get(n - 1, 0)
        dim_nm2 = dims.get(n - 2, 0)
        dim_nm3 = dims.get(n - 3, 0)

        results.append({
            'bettis': dict(sorted(bettis.items())),
            'b_nm2': b_nm2,
            'b_nm1': b_nm1,
            'b_nm3': b_nm3,
            'dim_nm1': dim_nm1,
            'dim_nm2': dim_nm2,
            'dim_nm3': dim_nm3,
        })

    elapsed = time.time() - t0
    print(f"  {num_samples} tournaments, {elapsed:.1f}s")

    # Statistics
    b_nm2_vals = [r['b_nm2'] for r in results]
    b_nm1_vals = [r['b_nm1'] for r in results]
    b_nm3_vals = [r['b_nm3'] for r in results]
    dim_nm1_vals = [r['dim_nm1'] for r in results]
    dim_nm2_vals = [r['dim_nm2'] for r in results]

    from collections import Counter
    print(f"\n  beta_{n-2} distribution: {dict(sorted(Counter(b_nm2_vals).items()))}")
    print(f"  beta_{n-1} distribution: {dict(sorted(Counter(b_nm1_vals).items()))}")
    print(f"  beta_{n-3} distribution: {dict(sorted(Counter(b_nm3_vals).items()))}")

    nonzero_rate = sum(1 for v in b_nm2_vals if v > 0) / num_samples
    print(f"\n  beta_{n-2} nonzero rate: {nonzero_rate:.1%} ({sum(1 for v in b_nm2_vals if v > 0)}/{num_samples})")
    if b_nm2_vals:
        print(f"  beta_{n-2} mean: {np.mean(b_nm2_vals):.2f}, max: {max(b_nm2_vals)}")

    print(f"\n  dim(Omega_{n-1}): min={min(dim_nm1_vals)}, max={max(dim_nm1_vals)}, mean={np.mean(dim_nm1_vals):.1f}")
    print(f"  dim(Omega_{n-2}): min={min(dim_nm2_vals)}, max={max(dim_nm2_vals)}, mean={np.mean(dim_nm2_vals):.1f}")

    # Check relationship between beta_{n-2} and Hamiltonian path count
    # dim(Omega_{n-1}) = number of Hamiltonian paths
    print(f"\n  Correlation: beta_{n-2} vs dim(Omega_{n-1}) (HP count):")
    hp_betas = [(r['dim_nm1'], r['b_nm2']) for r in results]
    hp_groups = {}
    for hp, b in hp_betas:
        hp_groups.setdefault(hp, []).append(b)
    for hp in sorted(hp_groups.keys())[:10]:
        vals = hp_groups[hp]
        mean_b = np.mean(vals)
        nonzero = sum(1 for v in vals if v > 0)
        print(f"    HP={hp}: mean_beta_{n-2}={mean_b:.2f}, nonzero={nonzero}/{len(vals)}")

    # Print a few full Betti sequences
    print(f"\n  Sample Betti sequences:")
    for i in range(min(5, num_samples)):
        print(f"    {results[i]['bettis']}")

    return results


def main():
    print("BETA_{n-2} INVESTIGATION")
    print("=" * 60)

    random.seed(42)
    np.random.seed(42)

    # n=5,6,7: should all have beta_{n-2}=0
    for n in [5, 6, 7]:
        investigate_top_degrees(n, num_samples=200)

    # n=8: the interesting case
    investigate_top_degrees(8, num_samples=300)


if __name__ == '__main__':
    main()
    print("\nDONE.")
