"""
Euler characteristic test: chi = sum (-1)^k dim(Omega_k) for random tournaments.

For each n in {3,4,5,6,7,8}, sample tournaments and compute:
  chi_omega = sum_{k=0}^{n-1} (-1)^k * omega_dims[k]
  chi_betti = sum_{k=0}^{n-1} (-1)^k * bettis[k]

Verify they match (standard homological algebra identity) and check
whether chi is universal (same for all tournaments at given n).
"""

import sys
import numpy as np
from collections import Counter

sys.path.insert(0, '04-computation')
from tournament_utils import random_tournament, full_chain_complex_modp


def compute_euler_chars(A, n):
    """Compute chi from omega_dims and from bettis, return both."""
    result = full_chain_complex_modp(A, n, max_p=n-1)
    omega_dims = result['omega_dims']
    bettis = result['bettis']

    chi_omega = 0
    chi_betti = 0
    for k in range(n):
        sign = (-1) ** k
        chi_omega += sign * omega_dims.get(k, 0)
        chi_betti += sign * bettis.get(k, 0)

    return chi_omega, chi_betti, omega_dims, bettis


def main():
    print("=" * 70)
    print("EULER CHARACTERISTIC TEST")
    print("chi = sum_{k=0}^{n-1} (-1)^k dim(Omega_k)")
    print("Verify: chi_omega == chi_betti for all tournaments")
    print("Check: Is chi universal (same for all T at given n)?")
    print("=" * 70)

    rng = np.random.RandomState(42)

    configs = {
        3: 500,
        4: 500,
        5: 500,
        6: 500,
        7: 500,
        8: 200,
    }

    for n in sorted(configs.keys()):
        num_samples = configs[n]
        print(f"\n{'='*60}")
        print(f"n = {n}, samples = {num_samples}, max_p = {n-1}")
        print(f"{'='*60}")

        chi_counter = Counter()
        mismatch_count = 0
        all_omega_profiles = []
        all_betti_profiles = []

        for trial in range(num_samples):
            A = random_tournament(n, rng)
            chi_o, chi_b, od, bt = compute_euler_chars(A, n)

            chi_counter[chi_o] += 1

            if chi_o != chi_b:
                mismatch_count += 1
                if mismatch_count <= 3:
                    print(f"  MISMATCH trial {trial}: chi_omega={chi_o}, chi_betti={chi_b}")
                    print(f"    omega_dims={dict(od)}")
                    print(f"    bettis={dict(bt)}")

            all_omega_profiles.append(tuple(od.get(k, 0) for k in range(n)))
            all_betti_profiles.append(tuple(bt.get(k, 0) for k in range(n)))

            if trial == 0:
                # Print detailed info for first tournament
                print(f"  Example (first tournament):")
                print(f"    omega_dims = {dict(od)}")
                print(f"    bettis     = {dict(bt)}")
                print(f"    chi_omega  = {chi_o}")
                print(f"    chi_betti  = {chi_b}")

        # Report
        print(f"\n  --- Results for n={n} ---")
        print(f"  Mismatches (chi_omega != chi_betti): {mismatch_count}/{num_samples}")

        is_universal = (len(chi_counter) == 1)
        print(f"  chi is UNIVERSAL: {is_universal}")
        print(f"  chi distribution:")
        for val in sorted(chi_counter.keys()):
            count = chi_counter[val]
            pct = 100.0 * count / num_samples
            print(f"    chi = {val:4d}: {count:5d} / {num_samples} ({pct:6.2f}%)")

        if is_universal:
            unique_val = list(chi_counter.keys())[0]
            print(f"  Universal value: chi = {unique_val}")
            if unique_val == n:
                print(f"  ** chi == n (matches n={n}) **")
            elif unique_val == 1:
                print(f"  ** chi == 1 **")
            else:
                print(f"  ** chi == {unique_val} (not n, not 1) **")

        # Omega profile diversity
        omega_unique = len(set(all_omega_profiles))
        betti_unique = len(set(all_betti_profiles))
        print(f"  Distinct omega_dim profiles: {omega_unique}")
        print(f"  Distinct betti profiles: {betti_unique}")

        # Show a few omega profiles if not too many
        if omega_unique <= 20:
            oc = Counter(all_omega_profiles)
            print(f"  Omega_dim profile distribution (k=0..{n-1}):")
            for prof in sorted(oc.keys(), key=lambda x: -oc[x]):
                print(f"    {prof}: {oc[prof]} times")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print("="*70)
    print("For each n, chi_omega should always equal chi_betti (homological identity).")
    print("If chi is universal, it depends only on n, not on the tournament.")
    print("Done.")


if __name__ == '__main__':
    main()
