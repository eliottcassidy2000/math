"""
top_exactness_test.py — Test exactness at degree n-2 for tournaments

Finding from beta6_n8_test: ker(d_{n-2}) = rk(d_{n-1}) = dim(Omega_{n-1}) for ALL 50 n=8 tournaments.
This means the top complex is EXACT at degree n-2:
  Omega_{n-1} --d_{n-1}--> Omega_{n-2} --d_{n-2}--> Omega_{n-3}
  ker(d_{n-2}) = im(d_{n-1})  =>  β_{n-2} = 0

Stronger: dim(Omega_{n-1}) = ker(d_{n-2}), meaning d_{n-1} is INJECTIVE on Omega_{n-1}.

Note: Omega_{n-1} = allowed (n-1)-paths with all faces allowed = Hamiltonian paths
      whose all subpaths are allowed.

Test this at n=5,6,7,8 with larger samples and verify both:
1. β_{n-2} = 0 always
2. d_{n-1} is injective (rk(d_{n-1}) = dim(Omega_{n-1}))

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp


def test_top_exactness(n, num_samples=200):
    """Test exactness at degree n-2."""
    print(f"\n{'='*60}")
    print(f"TOP EXACTNESS AT n={n}: degree {n-2}")
    print(f"{'='*60}")

    rng = np.random.RandomState(42)
    t0 = time.time()

    ker_nm2_vals = []
    rk_nm1_vals = []
    omega_nm1_vals = []
    beta_nm2_vals = []
    beta_nm1_vals = []
    injectivity_count = 0
    exactness_count = 0

    for trial in range(num_samples):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p=n-1)

        ker_nm2 = cc['kers'].get(n-2, 0)
        rk_nm1 = cc['ranks'].get(n-1, 0)
        omega_nm1 = cc['omega_dims'].get(n-1, 0)
        beta_nm2 = cc['bettis'].get(n-2, 0)
        beta_nm1 = cc['bettis'].get(n-1, 0)

        ker_nm2_vals.append(ker_nm2)
        rk_nm1_vals.append(rk_nm1)
        omega_nm1_vals.append(omega_nm1)
        beta_nm2_vals.append(beta_nm2)
        beta_nm1_vals.append(beta_nm1)

        if rk_nm1 == omega_nm1:
            injectivity_count += 1
        if ker_nm2 == rk_nm1:
            exactness_count += 1

    elapsed = time.time() - t0
    print(f"  {num_samples} tournaments, {elapsed:.1f}s")

    print(f"\n  beta_{n-2}: {dict(sorted(Counter(beta_nm2_vals).items()))}")
    print(f"  beta_{n-1}: {dict(sorted(Counter(beta_nm1_vals).items()))}")
    print(f"  dim(Omega_{n-1}) range: [{min(omega_nm1_vals)}, {max(omega_nm1_vals)}], mean={np.mean(omega_nm1_vals):.1f}")
    print(f"  rk(d_{n-1}) range: [{min(rk_nm1_vals)}, {max(rk_nm1_vals)}], mean={np.mean(rk_nm1_vals):.1f}")
    print(f"  ker(d_{n-2}) range: [{min(ker_nm2_vals)}, {max(ker_nm2_vals)}], mean={np.mean(ker_nm2_vals):.1f}")

    print(f"\n  d_{n-1} INJECTIVE: {injectivity_count}/{num_samples}")
    print(f"  EXACT at degree {n-2}: {exactness_count}/{num_samples}")
    print(f"  beta_{n-2} = 0: {sum(1 for v in beta_nm2_vals if v == 0)}/{num_samples}")

    # Additional: check if omega_{n-1} = ker(d_{n-2}) always
    match = sum(1 for k, o in zip(ker_nm2_vals, omega_nm1_vals) if k == o)
    print(f"  ker(d_{n-2}) = dim(Omega_{n-1}): {match}/{num_samples}")

    # Check for non-trivial Omega_{n-1}
    nontrivial = sum(1 for o in omega_nm1_vals if o > 0)
    print(f"  Omega_{n-1} nontrivial: {nontrivial}/{num_samples}")


def main():
    print("TOP-DEGREE EXACTNESS TEST")
    print("=" * 60)

    for n in [5, 6, 7, 8]:
        samples = 500 if n <= 7 else 200
        test_top_exactness(n, samples)


if __name__ == '__main__':
    main()
    print("\nDONE.")
