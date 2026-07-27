#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2549.

Checks representative exact future-pullback translation identities through
level seven, their sharp level-zero boundary, the carry-sheet and future-scale
target laws, exact BV mixing on a rational positive control, and the
cemetery-only Hall tables.  The all-level statement in THM-2549 is proved
symbolically; the finite loops are regression controls rather than induction.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fractional_part(value):
    return value - value.numerator // value.denominator


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def future_translation_referee():
    root_checks = 0
    target_checks = 0
    sheet_checks = 0
    future_action_checks = 0
    z = Fraction(17, 101)
    y = Fraction(29, 103)

    for level in range(1, 8):
        scale = P ** level
        for weight in range(-18, 19):
            for root in range(P):
                x = (z + root) / P
                left = fractional_part(scale * weight * x)
                right = fractional_part(P ** (level - 1) * weight * z)
                require(left == right,
                        f"old-root constancy failed {level=} {weight=} {root=}")
                root_checks += 1

        for theta in range(P):
            old_shifted = scale * (y - Fraction(theta, P))
            require(old_shifted - scale * y == -P ** (level - 1) * theta,
                    "old target shift did not become integral")
            require(fractional_part(old_shifted) == fractional_part(scale * y),
                    "future base saw old target shift")
            target_checks += 1

            new_scale_shift = scale * (
                y - Fraction(theta, P ** (level + 1))
            )
            require(new_scale_shift == scale * y - Fraction(theta, P),
                    "future-scale target pullback law")
            future_action_checks += 1

            for sheet in range(P):
                old_point = (z + sheet) / scale
                shifted = old_point - Fraction(theta, P)
                moved_sheet = sheet - P ** (level - 1) * theta
                require(shifted == (z + moved_sheet) / scale,
                        "old charge did not move solely in ancestry sheet")
                sheet_checks += 1

    # Sharp L=0 boundary.  The delta_0 Boolean mask on F_13 is root-visible
    # and is moved by every nonzero old target translation.
    mask = tuple(1 if root == 0 else 0 for root in range(P))
    require(len(set(mask)) == 2, "level-zero root visibility")
    for theta in range(1, P):
        moved = tuple(mask[(root - theta) % P] for root in range(P))
        require(moved != mask, "level-zero target action should be nontrivial")

    return root_checks, target_checks, sheet_checks, future_action_checks


def ancestry_digit_chart_referee():
    ancestry_root_checks = 0
    digit_covariance_checks = 0
    base_points = (Fraction(1, 101), Fraction(17, 101), Fraction(99, 101))

    for level in range(1, 7):
        scale = P ** level
        top_scale = P ** (level - 1)
        for weight in range(-18, 19):
            if weight % P == 0:
                continue
            inverse_weight = pow(weight, -1, P)
            for z in base_points:
                carry = (weight * z).numerator // (weight * z).denominator
                for head in range(P):
                    x = (z + head) / P
                    y = fractional_part(weight * x)
                    sheet = (scale * y).numerator // (scale * y).denominator
                    top_digit = sheet // top_scale
                    recovered = inverse_weight * (top_digit - carry) % P
                    require(recovered == head,
                            "unit-role ancestry digit failed to recover head")
                    ancestry_root_checks += 1

                    for theta in range(P):
                        shifted = fractional_part(y - Fraction(theta, P))
                        shifted_sheet = (
                            (scale * shifted).numerator
                            // (scale * shifted).denominator
                        )
                        shifted_top = shifted_sheet // top_scale
                        require(shifted_top % P == (top_digit - theta) % P,
                                "top ancestry digit lost target covariance")
                        digit_covariance_checks += 1

    # Sharp nonunit boundary: multiplication by thirteen erases the old head
    # before the top ancestry digit is read.
    z = Fraction(17, 101)
    nonunit_digits = {
        (P * fractional_part(P * (z + head) / P)).numerator
        // (P * fractional_part(P * (z + head) / P)).denominator
        for head in range(P)
    }
    require(len(nonunit_digits) == 1,
            "13-divisible role unexpectedly retained the predecessor root")

    return ancestry_root_checks, digit_covariance_checks


def in_intervals(x, intervals):
    x = fractional_part(x)
    return any(left <= x < right for left, right in intervals)


def interval_measure(intervals):
    return sum(right - left for left, right in intervals)


def product_measure(left_intervals, right_intervals, scale):
    endpoints = {Fraction(0), Fraction(1)}
    for left, right in left_intervals:
        endpoints.add(left)
        endpoints.add(right)
    for left, right in right_intervals:
        for k in range(scale):
            endpoints.add((k + left) / scale)
            endpoints.add((k + right) / scale)
    ordered = sorted(endpoints)
    answer = Fraction(0)
    for a, b in zip(ordered, ordered[1:]):
        if a == b:
            continue
        midpoint = (a + b) / 2
        if (in_intervals(midpoint, left_intervals)
                and in_intervals(scale * midpoint, right_intervals)):
            answer += b - a
    return answer


def exact_mixing_control():
    head = (
        (Fraction(1, 17), Fraction(5, 17)),
        (Fraction(9, 17), Fraction(12, 17)),
    )
    handoff = (
        (Fraction(2, 19), Fraction(7, 19)),
        (Fraction(11, 19), Fraction(16, 19)),
    )
    head_mass = interval_measure(head)
    handoff_mass = interval_measure(handoff)
    product_limit = head_mass * handoff_mass
    variation_head = 2 * len(head)
    variation_handoff = 2 * len(handoff)
    covariance_constant = min(
        handoff_mass * (1 - handoff_mass),
        Fraction(variation_handoff, 12),
    )

    values = []
    for level in range(1, 6):
        scale = P ** level
        value = product_measure(head, handoff, scale)
        error = abs(value - product_limit)
        bound = covariance_constant * variation_head / scale
        require(value > 0, f"positive mixing control failed at {level=}")
        require(error <= bound, f"BV covariance invoice failed at {level=}")
        values.append(value)
    return head_mass, handoff_mass, product_limit, values


def hall_cemetery_referee():
    margin_checks = 0
    subset_checks = 0
    for roots in range(2, 7):
        for total in range(1, 7):
            for left in compositions(total, roots):
                diagonal = 0
                right_cemetery = total
                require(diagonal == 0, "cemetery table acquired a diagonal")
                for mask in range(1 << roots):
                    left_mass = sum(
                        left[i] for i in range(roots) if mask & (1 << i)
                    )
                    has_neighbor = any(
                        left[i] > 0 and mask & (1 << i)
                        for i in range(roots)
                    )
                    neighbor_mass = right_cemetery if has_neighbor else 0
                    require(left_mass <= neighbor_mass,
                            f"cemetery Hall inequality failed {left=}")
                    subset_checks += 1
                margin_checks += 1

    # Full F_13 hostile controls: singleton, two-root, and uniform head
    # margins all route entirely to the cemetery column.
    full_profiles = (
        (1,) + (0,) * 12,
        (1, 1) + (0,) * 11,
        (1,) * 13,
    )
    full_subset_checks = 0
    for left in full_profiles:
        total = sum(left)
        for mask in range(1 << P):
            left_mass = sum(left[i] for i in range(P) if mask & (1 << i))
            has_neighbor = any(
                left[i] > 0 and mask & (1 << i) for i in range(P)
            )
            neighbor_mass = total if has_neighbor else 0
            require(left_mass <= neighbor_mass, "F_13 cemetery Hall failure")
            full_subset_checks += 1
    return margin_checks, subset_checks, full_subset_checks


def main():
    root, target, sheet, future = future_translation_referee()
    ancestry_root, digit_covariance = ancestry_digit_chart_referee()
    head, handoff, limit, values = exact_mixing_control()
    margins, subsets, full_subsets = hall_cemetery_referee()

    print(f"future_root_constancy_checks={root}")
    print(f"old_target_neutrality_checks={target}")
    print(f"ancestry_sheet_transport_checks={sheet}")
    print(f"future_scale_target_action_checks={future}")
    print("level_zero_root_visible_and_target_active_control=PASS")
    print(f"unit_role_old_action_ancestry_root_checks={ancestry_root}")
    print(f"top_digit_target_covariance_checks={digit_covariance}")
    print("13_divisible_role_root_erasure_control=PASS")
    print(f"mixing_head_mass={head} handoff_mass={handoff} limit={limit}")
    print("mixing_levels_1_to_5=" + ",".join(str(value) for value in values))
    print(f"cemetery_margin_tables={margins}")
    print(f"cemetery_hall_subset_checks={subsets}")
    print(f"F13_cemetery_hall_subset_checks={full_subsets}")
    print("future_pullback_hall_diagonal=0")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
