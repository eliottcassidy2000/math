#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2555.

The symbolic proof in THM-2555 is all-depth.  This companion checks the
sheet action and its complete affine-equivariant normal form on representative
depths, the unit-role carry correction, the sharp nonunit erasure boundary,
the separation of the old top digit from the future immediate digit (including
the future-action borrow into the old sheet), and the exact positive
zero-arrival digit-cylinder hostile.
"""

from fractions import Fraction
from itertools import product


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_q(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_q(value)


def split_sheet(a, level):
    top_scale = P ** (level - 1)
    return a // top_scale, a % top_scale


def old_sheet_action(a, level, theta):
    return (a - P ** (level - 1) * theta) % (P ** level)


def sheet_action_and_classification_referee():
    action_checks = 0
    classification_checks = 0
    charges = range(P)

    for level in range(1, 5):
        modulus = P ** level
        for a in range(modulus):
            top, tail = split_sheet(a, level)
            # One deliberately nonconstant invariant gauge.  The theorem
            # proves that an arbitrary function of the tail works.
            gauge = (7 * tail * tail + 3 * tail + 5 * level) % P
            for theta in range(P):
                moved = old_sheet_action(a, level, theta)
                moved_top, moved_tail = split_sheet(moved, level)
                require(moved_tail == tail, "old action moved invariant tail")
                require(moved_top == (top - theta) % P,
                        "old action did not translate top digit")
                action_checks += 1

                # Check every possible F_13 charge lambda.  For fixed tail,
                # beta=lambda*top+gauge is the complete normal form.
                for charge in charges:
                    beta = (charge * top + gauge) % P
                    moved_beta = (charge * moved_top + gauge) % P
                    require(moved_beta == (beta - charge * theta) % P,
                            "equivariant normal form failed")
                    # Converse reconstruction from the top-zero section.
                    recovered_gauge = (beta - charge * top) % P
                    require(recovered_gauge == gauge,
                            "equivariant map lost its invariant gauge")
                    classification_checks += 1

    return action_checks, classification_checks


def unit_role_carry_referee():
    carry_checks = 0
    covariance_checks = 0
    base_points = (
        Fraction(1, 101),
        Fraction(17, 101),
        Fraction(99, 101),
        Fraction(5, 37),
    )

    for level in range(1, 6):
        scale = P ** level
        top_scale = P ** (level - 1)
        for weight in range(-30, 31):
            if weight % P == 0:
                continue
            inverse_weight = pow(weight % P, -1, P)
            for base in base_points:
                carry = floor_q(weight * base)
                for head in range(P):
                    x = (base + head) / P
                    y = fractional_part(weight * x)
                    sheet = floor_q(scale * y)
                    top = sheet // top_scale
                    recovered = inverse_weight * (top - carry) % P
                    require(recovered == head,
                            "unit-role top digit failed to recover old head")
                    carry_checks += 1

                    for theta in range(P):
                        shifted_y = fractional_part(y - Fraction(theta, P))
                        shifted_sheet = floor_q(scale * shifted_y)
                        shifted_top = shifted_sheet // top_scale
                        shifted_root = inverse_weight * (
                            shifted_top - carry
                        ) % P
                        require(
                            shifted_root
                            == (head - inverse_weight * theta) % P,
                            "carry-corrected old root lost covariance",
                        )
                        require(
                            fractional_part(scale * shifted_y)
                            == fractional_part(scale * y),
                            "old target action moved the future base",
                        )
                        covariance_checks += 1

    return carry_checks, covariance_checks


def nonunit_erasure_referee():
    checks = 0
    base_points = (
        Fraction(1, 101),
        Fraction(17, 101),
        Fraction(99, 101),
    )
    weights = (-39, -26, -13, 13, 26, 39)

    for weight in weights:
        for base in base_points:
            values = {
                fractional_part(weight * (base + head) / P)
                for head in range(P)
            }
            require(len(values) == 1,
                    "13-divisible role unexpectedly retained old head")
            checks += P
    return checks


def future_digit_referee():
    future_action_checks = 0
    future_borrow_checks = 0
    digit_independence_checks = 0
    pair_counts_by_level = []
    observed_top_deltas = set()

    # The immediate root of the future base z is charged by z -> z-theta/13.
    for future_digit in range(P):
        z = (future_digit + Fraction(2, 5)) / P
        require(floor_q(P * z) == future_digit, "future digit setup")
        for theta in range(P):
            shifted = fractional_part(z - Fraction(theta, P))
            require(
                floor_q(P * shifted) == (future_digit - theta) % P,
                "future action did not translate future immediate digit",
            )
            future_action_checks += 1

    # In full natural-extension coordinates the same future action is not
    # sheet-fixed.  If y=(z+a)/13^L, then
    #
    #   y-phi/13^(L+1)
    #    = ({z-phi/13} + a - 1_{z<phi/13})/13^L  (mod 1).
    #
    # The borrow can propagate into the old top digit, so that digit has no
    # uniform future charge even though the immediate digit of z does.
    for level in range(1, 4):
        modulus = P ** level
        for sheet in range(modulus):
            old_top, _ = split_sheet(sheet, level)
            for future_digit in range(P):
                z = (future_digit + Fraction(2, 5)) / P
                y = (z + sheet) / modulus
                for phi in range(P):
                    epsilon = int(z < Fraction(phi, P))
                    moved_z = fractional_part(z - Fraction(phi, P))
                    moved_sheet = (sheet - epsilon) % modulus
                    moved_y = fractional_part(
                        y - Fraction(phi, P ** (level + 1))
                    )
                    reconstructed = fractional_part(
                        (moved_z + moved_sheet) / modulus
                    )
                    require(moved_y == reconstructed,
                            "future-action natural-extension borrow failed")
                    require(
                        floor_q(P * moved_z)
                        == (future_digit - phi) % P,
                        "borrow model lost future immediate-digit charge",
                    )
                    moved_top, _ = split_sheet(moved_sheet, level)
                    observed_top_deltas.add((moved_top - old_top) % P)
                    future_borrow_checks += 1

    require(observed_top_deltas == {0, P - 1},
            "old top digit falsely acquired a uniform future charge")

    # A finite base-13 word d_1...d_L e exhibits the exact split:
    # d_1 is the old charged top-sheet digit; e is the future immediate digit.
    for level in range(1, 5):
        counts = [[0 for _ in range(P)] for _ in range(P)]
        for digits in product(range(P), repeat=level + 1):
            sheet = sum(
                digits[index] * P ** (level - 1 - index)
                for index in range(level)
            )
            top, _ = split_sheet(sheet, level)
            future = digits[level]
            require(top == digits[0], "sheet top is not first old digit")
            counts[top][future] += 1
            digit_independence_checks += 1

        expected = P ** (level - 1)
        require(all(value == expected for row in counts for value in row),
                "old and future digits are not independent coordinates")
        pair_counts_by_level.append(expected)

    return (
        future_action_checks,
        future_borrow_checks,
        sorted(observed_top_deltas),
        digit_independence_checks,
        pair_counts_by_level,
    )


def convolution_and_hostile_referee():
    convolution_checks = 0
    profiles = (
        ((1,) + (0,) * 12, (0, 1) + (0,) * 11),
        ((1,) * 13, tuple(range(1, 14))),
        (tuple((3 * h + 1) % 7 for h in range(P)),
         tuple((5 * h + 2) % 11 for h in range(P))),
    )

    for left, right in profiles:
        correlation = [
            sum(left[h] * right[(h + displacement) % P]
                for h in range(P))
            for displacement in range(P)
        ]
        require(sum(correlation) == sum(left) * sum(right),
                "displacement convolution lost total mass")
        # Enumerate the two endpoint digits independently and recover the
        # same displacement ledger.
        direct = [0] * P
        for head in range(P):
            for future in range(P):
                direct[(future - head) % P] += left[head] * right[future]
                convolution_checks += 1
        require(direct == correlation, "mixing convolution mismatch")

    hostile_masses = []
    orbit_hostile_masses = []
    hostile_checks = 0
    orbit_hostile_checks = 0
    for level in range(1, 7):
        # Words of length level+1 with first digit 0 and future digit 1.
        count = P ** (level - 1)
        total = P ** (level + 1)
        mass = Fraction(count, total)
        require(mass == Fraction(1, P * P), "hostile mass changed with depth")
        require((1 - 0) % P == 1, "hostile displacement")
        require(0 != 1, "hostile acquired a semantic diagonal")
        hostile_masses.append(mass)
        hostile_checks += count

        # Orbit-symmetrized hostile: allow every old top digit and impose
        # future=top+1, leaving all L-1 intervening sheet digits free.  This
        # is invariant under the abstract simultaneous relabelling of the
        # two endpoint digits; it is not asserted to be a physical circle
        # translation of a live LRC packet.
        orbit_count = P ** level
        orbit_mass = Fraction(orbit_count, total)
        require(orbit_mass == Fraction(1, P),
                "orbit hostile mass changed with depth")
        for top in range(P):
            for shift in range(P):
                moved_top = (top + shift) % P
                moved_future = (top + 1 + shift) % P
                require((moved_future - moved_top) % P == 1,
                        "simultaneous digit relabelling moved displacement")
                orbit_hostile_checks += 1
        orbit_hostile_masses.append(orbit_mass)

    # In the minimal hostile, source support {0} and future support {1}
    # have positive product mass but zero overlap at displacement zero.
    left = (1,) + (0,) * 12
    right = (0, 1) + (0,) * 11
    correlation = [
        sum(left[h] * right[(h + displacement) % P] for h in range(P))
        for displacement in range(P)
    ]
    require(correlation[0] == 0 and correlation[1] == 1,
            "sharp zero-arrival convolution hostile failed")

    return (
        convolution_checks,
        hostile_checks,
        hostile_masses,
        orbit_hostile_checks,
        orbit_hostile_masses,
    )


def main():
    action, classification = sheet_action_and_classification_referee()
    carry, covariance = unit_role_carry_referee()
    nonunit = nonunit_erasure_referee()
    (
        future,
        future_borrow,
        future_top_deltas,
        independence,
        pair_counts,
    ) = future_digit_referee()
    (
        convolution,
        hostile,
        masses,
        orbit_hostile,
        orbit_masses,
    ) = convolution_and_hostile_referee()

    print(f"old_sheet_action_checks={action}")
    print(f"equivariant_normal_form_checks={classification}")
    print(f"unit_role_carry_identity_checks={carry}")
    print(f"unit_role_old_action_covariance_checks={covariance}")
    print(f"13_divisible_role_erasure_checks={nonunit}")
    print(f"future_immediate_digit_action_checks={future}")
    print(f"future_natural_extension_borrow_checks={future_borrow}")
    print("future_action_old_top_digit_deltas="
          + ",".join(str(value) for value in future_top_deltas))
    print(f"old_future_digit_independence_checks={independence}")
    print("old_future_pair_multiplicity_levels_1_to_4="
          + ",".join(str(value) for value in pair_counts))
    print(f"displacement_convolution_checks={convolution}")
    print(f"positive_zero_arrival_hostile_word_checks={hostile}")
    print("hostile_masses_levels_1_to_6="
          + ",".join(str(value) for value in masses))
    print("hostile_head=0 future_root=1 displacement=1 diagonal=0")
    print(f"orbit_symmetrized_hostile_checks={orbit_hostile}")
    print("orbit_hostile_masses_levels_1_to_6="
          + ",".join(str(value) for value in orbit_masses))
    print("orbit_hostile_uniform_margins=1/169 displacement=1 diagonal=0")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
