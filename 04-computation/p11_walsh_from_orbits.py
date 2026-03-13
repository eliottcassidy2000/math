#!/usr/bin/env python3
"""
p11_walsh_from_orbits.py -- Compute Walsh coefficients at p=11 using orbit structure

Since we know there are exactly 4 extended orbits with known (alpha_j, H) values,
we can compute the Walsh decomposition without waiting for all 32 orientations.

We know which bits belong to which orbit (from orientation_symmetry_p11.py),
so we can assign the correct (alpha_1, alpha_2, alpha_3, H) to each bits value
and compute Walsh coefficients.

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def main():
    p = 11
    m = 5
    half = 32

    # From the computation + orbit analysis:
    # Extended orbits (bits -> orbit class):
    orbit_members = {
        'A': [0b00000, 0b00011, 0b00110, 0b01010, 0b10010,  # orbit 0 -> Interval class
              0b01101, 0b10101, 0b11001, 0b11100, 0b11111],
        'B': [0b00001, 0b01000, 0b01110, 0b10011, 0b10110,  # orbit 1 -> lowest H
              0b01001, 0b01100, 0b10001, 0b10111, 0b11110],
        'P': [0b00010, 0b11101],                              # orbit 2 -> Paley
        'C': [0b00100, 0b00111, 0b01011, 0b10000, 0b11010,  # orbit 3
              0b00101, 0b01111, 0b10100, 0b11000, 0b11011],
    }

    # Verify: all 32 bits covered
    all_bits = set()
    for members in orbit_members.values():
        all_bits.update(members)
    assert len(all_bits) == 32 and all_bits == set(range(32))

    # Alpha values per orbit (from computation):
    orbit_data = {
        'A': {'alpha_1': 18397, 'alpha_2': 11110, 'alpha_3': 1474, 'H': 93027,
               'c5': 484, 'c7': 3399, 'c9': 9350, 'c11': 5109},
        'B': {'alpha_1': 19629, 'alpha_2': 10912, 'alpha_3': 1188, 'H': 92411,
               'c5': 572, 'c7': 3729, 'c9': 10274, 'c11': 4999},
        'P': {'alpha_1': 21169, 'alpha_2': 10879, 'alpha_3': 1155, 'H': 95095,
               'c5': 594, 'c7': 3960, 'c9': 11055, 'c11': 5505},
        'C': {'alpha_1': 19541, 'alpha_2': 11220, 'alpha_3': 1188, 'H': 93467,
               'c5': 550, 'c7': 3586, 'c9': 10197, 'c11': 5153},
    }

    # Build per-bits data
    data = {}
    for label, members in orbit_members.items():
        for bits in members:
            data[bits] = orbit_data[label]

    # Verify OCF: H = 1 + 2*a1 + 4*a2 + 8*a3
    for label, d in orbit_data.items():
        H_check = 1 + 2*d['alpha_1'] + 4*d['alpha_2'] + 8*d['alpha_3']
        assert H_check == d['H'], f"OCF mismatch for {label}: {H_check} != {d['H']}"
    print("OCF verification passed for all orbits.")

    # Walsh decomposition
    print("\n" + "=" * 70)
    print("WALSH DECOMPOSITION AT p=11 (from orbit structure)")
    print("=" * 70)

    quantities = ['H', 'alpha_1', 'alpha_2', 'alpha_3', 'c5', 'c7', 'c9', 'c11']

    for name in quantities:
        vals = [data[bits][name] for bits in range(half)]
        mean = sum(vals) / half

        print(f"\n  {name}:")
        print(f"    Mean: {mean:.4f}")
        print(f"    Values unique: {sorted(set(vals))}")

        # Degree-1 Walsh
        has_deg1 = False
        for a in range(m):
            total = sum((1 if bits & (1 << a) else -1) * vals[bits]
                        for bits in range(half))
            h = total / half
            if abs(h) > 0.001:
                if not has_deg1:
                    print(f"    Degree-1:")
                    has_deg1 = True
                print(f"      h_hat[{a}] = {h:>12.4f}")

        # Degree-2 Walsh
        print(f"    Degree-2:")
        for a in range(m):
            for b in range(a + 1, m):
                total = sum((1 if bits & (1 << a) else -1) *
                            (1 if bits & (1 << b) else -1) * vals[bits]
                            for bits in range(half))
                h = total / half
                if abs(h) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    chi_ab = legendre(gap_a * gap_b, p)
                    sign_ok = (h > 0) == (chi_ab > 0) if h != 0 else True
                    print(f"      h_hat[{{{a},{b}}}] = {h:>12.4f}, "
                          f"gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}, "
                          f"{'OK' if sign_ok else 'FAIL'}")

        # Degree-3 Walsh
        has_deg3 = False
        for indices in combinations(range(m), 3):
            total = 0
            for bits in range(half):
                prod = 1
                for idx in indices:
                    prod *= (1 if bits & (1 << idx) else -1)
                total += vals[bits] * prod
            h = total / half
            if abs(h) > 0.001:
                if not has_deg3:
                    print(f"    Degree-3:")
                    has_deg3 = True
                print(f"      h_hat[{set(indices)}] = {h:>12.4f}")

        # Degree-4 Walsh
        has_deg4 = False
        for indices in combinations(range(m), 4):
            total = 0
            for bits in range(half):
                prod = 1
                for idx in indices:
                    prod *= (1 if bits & (1 << idx) else -1)
                total += vals[bits] * prod
            h = total / half
            if abs(h) > 0.001:
                if not has_deg4:
                    print(f"    Degree-4:")
                    has_deg4 = True
                gaps = tuple(i + 1 for i in indices)
                prod_gaps = 1
                for g in gaps:
                    prod_gaps *= g
                chi_prod = legendre(prod_gaps, p)
                sign_ok = (h > 0) == (chi_prod > 0) if h != 0 else True
                print(f"      h_hat[{set(indices)}] = {h:>12.4f}, "
                      f"gaps={gaps}, chi(prod)={chi_prod:+d}, "
                      f"{'OK' if sign_ok else 'FAIL'}")

        # Degree-5 Walsh (just the single one)
        total = 0
        for bits in range(half):
            prod = 1
            for idx in range(m):
                prod *= (1 if bits & (1 << idx) else -1)
            total += vals[bits] * prod
        h = total / half
        if abs(h) > 0.001:
            print(f"    Degree-5:")
            print(f"      h_hat[{{0,1,2,3,4}}] = {h:>12.4f}")

    # OCF decomposition at degree 2
    print(f"\n\n{'='*70}")
    print("OCF DECOMPOSITION: h_hat_H = 2*h_a1 + 4*h_a2 + 8*h_a3")
    print("=" * 70)

    for a in range(m):
        for b in range(a + 1, m):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)

            h_a1 = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) *
                        data[bits]['alpha_1'] for bits in range(half)) / half
            h_a2 = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) *
                        data[bits]['alpha_2'] for bits in range(half)) / half
            h_a3 = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) *
                        data[bits]['alpha_3'] for bits in range(half)) / half
            h_H = sum((1 if bits & (1<<a) else -1) * (1 if bits & (1<<b) else -1) *
                       data[bits]['H'] for bits in range(half)) / half

            recon = 2*h_a1 + 4*h_a2 + 8*h_a3

            print(f"\n  Pair ({a},{b}), gaps=({gap_a},{gap_b}), chi(ab)={chi_ab:+d}:")
            print(f"    h_hat_a1 = {h_a1:>12.4f} (sign {'+' if h_a1>0 else '-'})")
            print(f"    h_hat_a2 = {h_a2:>12.4f} (sign {'+' if h_a2>0 else '-'})")
            print(f"    h_hat_a3 = {h_a3:>12.4f} (sign {'+' if h_a3>0 else '-'})")
            print(f"    2*h_a1   = {2*h_a1:>12.4f}")
            print(f"    4*h_a2   = {4*h_a2:>12.4f}")
            print(f"    8*h_a3   = {8*h_a3:>12.4f}")
            print(f"    h_hat_H  = {h_H:>12.4f} (recon={recon:>12.4f})")

            total_contrib = abs(2*h_a1) + abs(4*h_a2) + abs(8*h_a3)
            if total_contrib > 0:
                print(f"    Contribution: a1={abs(2*h_a1)/total_contrib*100:.1f}% "
                      f"a2={abs(4*h_a2)/total_contrib*100:.1f}% "
                      f"a3={abs(8*h_a3)/total_contrib*100:.1f}%")

    # Variance decomposition
    print(f"\n\n{'='*70}")
    print("VARIANCE DECOMPOSITION")
    print("=" * 70)

    H_vals = [data[bits]['H'] for bits in range(half)]
    var_H = sum((h - sum(H_vals)/half)**2 for h in H_vals) / half

    # Walsh coefficients squared sum to variance
    total_walsh_var = 0
    for deg in range(1, m + 1):
        deg_var = 0
        for indices in combinations(range(m), deg):
            total = 0
            for bits in range(half):
                prod = 1
                for idx in indices:
                    prod *= (1 if bits & (1 << idx) else -1)
                total += H_vals[bits] * prod
            h = total / half
            deg_var += h**2
        print(f"  Degree-{deg} variance: {deg_var:>15.4f} ({deg_var/var_H*100:.2f}%)")
        total_walsh_var += deg_var

    print(f"  Total Walsh variance: {total_walsh_var:>15.4f}")
    print(f"  Direct variance:      {var_H:>15.4f}")
    print(f"  Match: {abs(total_walsh_var - var_H) < 0.01}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
