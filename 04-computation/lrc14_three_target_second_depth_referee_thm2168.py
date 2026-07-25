#!/usr/bin/env python3
"""Exact hostile controls for THM-2168.

This checks the load-bearing finite consequence rather than assuming its
classification: every pair of translated radius-1/14 character bands that
covers C_m x C_n (m|n, n<=18) is the complementary pair of fibres of one
order-two character.  It also prints an exact evaluation-kernel shear whose
scalar speed row is fixed while the transverse determinants drift.
"""

from fractions import Fraction
from math import gcd


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def mod_one(x):
    return x - (x.numerator // x.denominator)


def circle_norm(x):
    y = mod_one(x)
    return min(y, 1 - y)


def lcm(a, b):
    return a // gcd(a, b) * b


def character_order(m, n, p, q):
    return lcm(m // gcd(m, p), n // gcd(n, q))


def translated_value_masks(order):
    """All grid subsets cut out by a translated closed radius-1/14 arc."""
    critical = set()
    for j in range(order):
        v = Fraction(j, order)
        critical.add(mod_one(-v + Fraction(1, 14)))
        critical.add(mod_one(-v - Fraction(1, 14)))
    cuts = sorted(critical)
    probes = set(cuts)
    for i, left in enumerate(cuts):
        right = cuts[(i + 1) % len(cuts)]
        if i + 1 == len(cuts):
            right += 1
        probes.add(mod_one((left + right) / 2))
    masks = set()
    for theta in probes:
        mask = 0
        for j in range(order):
            if circle_norm(Fraction(j, order) + theta) <= Fraction(1, 14):
                mask |= 1 << j
        masks.add(mask)
    return masks


def finite_two_band_sweep():
    groups = 0
    characters = 0
    translated_states = 0
    half_states_total = 0
    covering_pairs = 0
    covering_groups = 0

    for n in range(1, 19):
        for m in range(1, n + 1):
            if n % m:
                continue
            size = m * n
            if size == 1:
                continue
            groups += 1
            elements = [(x, y) for x in range(m) for y in range(n)]
            half_states = set()

            for p in range(m):
                for q in range(n):
                    if p == 0 and q == 0:
                        continue
                    characters += 1
                    order = character_order(m, n, p, q)
                    value_index = []
                    for x, y in elements:
                        value = mod_one(
                            Fraction(p * x, m) + Fraction(q * y, n)
                        )
                        index = value * order
                        check(index.denominator == 1, "character order mismatch")
                        value_index.append(index.numerator % order)

                    for residue_mask in translated_value_masks(order):
                        translated_states += 1
                        group_mask = 0
                        for e, index in enumerate(value_index):
                            if (residue_mask >> index) & 1:
                                group_mask |= 1 << e
                        if group_mask.bit_count() * 2 == size:
                            half_states.add((tuple(value_index), order, group_mask))

            half_states_total += len(half_states)
            full = (1 << size) - 1
            local_pairs = 0
            states = sorted(half_states, key=lambda z: (z[1], z[0], z[2]))
            for i, (values_a, order_a, mask_a) in enumerate(states):
                for values_b, order_b, mask_b in states[i + 1 :]:
                    if (mask_a | mask_b) != full:
                        continue
                    local_pairs += 1
                    check(order_a == order_b == 2, "cover not order two")
                    check(values_a == values_b, "distinct order-two characters cover")
                    check((mask_a & mask_b) == 0, "covering half-fibres overlap")
                    check((mask_a ^ mask_b) == full, "half-fibres not complementary")
            if local_pairs:
                covering_groups += 1
                covering_pairs += local_pairs

    return {
        "groups": groups,
        "characters": characters,
        "translated_states": translated_states,
        "half_states": half_states_total,
        "covering_groups": covering_groups,
        "covering_pairs": covering_pairs,
    }


def shear_control():
    # Gamma=Z^2, h(x,y)=x, k=(0,1), g=(1,0).
    # The evaluated row (g,c,13u) is fixed at (1,2,39).
    scalar_row = (1, 2, 39)
    rows = []
    for n_shift, m_shift in ((0, 0), (1, 5), (12, -7), (1000, 2000)):
        c = (2, 13 + 13 * n_shift)
        u = (3, 1 + m_shift)
        evaluated = (1, c[0], 13 * u[0])
        check(evaluated == scalar_row, "shear changed scalar row")
        check((c[0] % 13, c[1] % 13) == (2, 0), "aligned type changed")
        blocker = (13 * u[0], 13 * u[1])
        check(blocker[0] % 13 == blocker[1] % 13 == 0, "blocker type changed")
        rows.append((n_shift, m_shift, c[1], u[1]))
    check(rows[-1][2] > 10000 and rows[-1][3] > 1000, "no determinant drift")
    return scalar_row, rows


def main():
    stats = finite_two_band_sweep()
    scalar_row, shear_rows = shear_control()
    print("THM-2168 EXACT HOSTILE CONTROLS")
    print("finite groups:", stats["groups"])
    print("nontrivial characters:", stats["characters"])
    print("translated band states:", stats["translated_states"])
    print("half-size states:", stats["half_states"])
    print("groups with two-band covers:", stats["covering_groups"])
    print("covering half-fibre pairs:", stats["covering_pairs"])
    print("all covers: same order-two character, complementary cosets = PASS")
    print("fixed evaluated row:", scalar_row)
    for row in shear_rows:
        print("shear N,M and determinants:", row)
    print("evaluation-shear blindness = PASS")


if __name__ == "__main__":
    main()
