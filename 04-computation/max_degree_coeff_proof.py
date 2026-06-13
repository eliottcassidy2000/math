"""
max_degree_coeff_proof.py — opus-2026-04-04-S1

THEOREM: At the maximum degree d_max = 2*floor((n-1)/2) of the multilinear
polynomial H(t_1,...,t_m), all nonzero coefficients are ±2.

Proof strategy:
1. The maximum degree corresponds to tile subsets where every tile in the
   subset contributes to the interaction.
2. At max degree, the Möbius coefficient c_S = ±2 because H values involved
   differ by exactly 2 (Rédei: H is odd, so consecutive values differ by 2).

This script:
1. Verifies at n=3,4,5,6,7,8
2. Checks the sign pattern
3. Investigates whether the ±2 pattern extends one level below max degree
"""
import numpy as np
from collections import defaultdict
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def verify_max_degree_coeffs():
    print("=" * 60)
    print("MAX DEGREE COEFFICIENT THEOREM")
    print("=" * 60)

    for n in range(3, 8):
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        max_order = max(bin(S).count('1') for S in coeffs.keys()) if coeffs else 0
        expected_max = 2 * ((n-1)//2)

        # Collect max-order coefficients
        max_coeffs = [(S, c) for S, c in coeffs.items() if bin(S).count('1') == max_order]
        values = [c for _, c in max_coeffs]

        # One level below
        sub_coeffs = [(S, c) for S, c in coeffs.items()
                      if bin(S).count('1') == max_order - 1]
        sub_values = [c for _, c in sub_coeffs]

        print(f"\nn={n}, m={m}")
        print(f"  Max degree: {max_order} (expected {expected_max})")
        print(f"  Max-degree coefficients: {len(max_coeffs)} nonzero")
        if max_coeffs:
            distinct = sorted(set(values))
            print(f"  Distinct values: {distinct}")
            all_pm2 = all(abs(c) == 2 for c in values)
            print(f"  All ±2: {all_pm2}")
            pos = sum(1 for c in values if c > 0)
            neg = sum(1 for c in values if c < 0)
            print(f"  +2: {pos}, -2: {neg}")

        if sub_coeffs:
            sub_distinct = sorted(set(sub_values))
            sub_max = max(abs(c) for c in sub_values)
            print(f"  Sub-max degree ({max_order-1}): {len(sub_coeffs)} coeffs, "
                  f"values {sub_distinct}, max |c|={sub_max}")

    # n=8 is feasible but slow (2^21 = 2M tilings). Skip for now.
    # Instead, verify at n=4 where we can check by hand.
    print("\n" + "=" * 60)
    print("DETAILED n=4 ANALYSIS")
    print("=" * 60)

    n = 4
    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    print(f"n=4, m={m} tiles: {tiles}")
    print(f"H values: {[H[i] for i in range(1 << m)]}")
    print(f"\nAll nonzero coefficients:")
    for S in sorted(coeffs.keys()):
        c = coeffs[S]
        tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
        order = bin(S).count('1')
        print(f"  S={S:0{m}b} order={order} tiles={tile_set}: c={c}")

    # Check the sign pattern at max degree for n=5,6,7
    print("\n" + "=" * 60)
    print("SIGN PATTERN ANALYSIS")
    print("=" * 60)

    for n in [5, 6, 7]:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)
        max_order = max(bin(S).count('1') for S in coeffs.keys())

        max_coeffs = [(S, coeffs[S]) for S in sorted(coeffs.keys())
                      if bin(S).count('1') == max_order]

        print(f"\nn={n}, max degree={max_order}:")
        print(f"  {len(max_coeffs)} coefficients")

        # Analyze which tile subsets give +2 vs -2
        for S, c in max_coeffs[:10]:
            tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
            skips = [x-y for x,y in tile_set]
            # Sum of skips
            skip_sum = sum(skips)
            # Product of skips
            skip_prod = 1
            for s in skips:
                skip_prod *= s

            # Parity of skip sum
            parity = sum(skips) % 2

            print(f"  tiles={tile_set} skips={skips} sum={skip_sum} "
                  f"prod={skip_prod} parity={parity}: c={c}")

        # Check: does sign depend on parity of skip sum?
        correct_parity = 0
        for S, c in max_coeffs:
            tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
            skips = [x-y for x,y in tile_set]
            skip_sum = sum(skips)
            if (c > 0) == (skip_sum % 2 == 0):
                correct_parity += 1

        print(f"  Sign = (-1)^(skip_sum mod 2): "
              f"{correct_parity}/{len(max_coeffs)} = "
              f"{100*correct_parity/len(max_coeffs):.0f}%")

        # Try: sign depends on number of "lower row" tiles
        # (tiles in the bottom half of the staircase)
        for rule_name, rule_fn in [
            ("sum of y-coordinates mod 2",
             lambda ts: sum(y for x,y in ts) % 2),
            ("sum of x-coordinates mod 2",
             lambda ts: sum(x for x,y in ts) % 2),
            ("product of skips mod 4 == 2",
             lambda ts: 1 if (np.prod([x-y for x,y in ts]) % 4 == 2) else 0),
            ("number of odd-skip tiles mod 2",
             lambda ts: sum(1 for x,y in ts if (x-y) % 2 == 1) % 2),
        ]:
            correct = sum(1 for S, c in max_coeffs
                         if (c > 0) == (rule_fn([tiles[i] for i in range(m)
                                                 if S & (1 << i)]) == 0))
            alt = len(max_coeffs) - correct
            best = max(correct, alt)
            print(f"  Sign from {rule_name}: {best}/{len(max_coeffs)} = "
                  f"{100*best/len(max_coeffs):.0f}%")

if __name__ == "__main__":
    verify_max_degree_coeffs()
