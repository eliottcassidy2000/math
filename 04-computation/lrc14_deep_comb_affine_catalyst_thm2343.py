#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2343."""

from fractions import Fraction
from itertools import product


P = 13
GROUP = tuple(product(range(P), repeat=2))
ZERO = (0, 0)
GROUP_SIZE = P**2


def require(condition: bool, message: str) -> None:
    """Raise under both ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def sub(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def negate(point: tuple[int, int]) -> tuple[int, int]:
    return ((-point[0]) % P, (-point[1]) % P)


def dot(left: tuple[int, int], right: tuple[int, int]) -> int:
    return (left[0] * right[0] + left[1] * right[1]) % P


def translation_orthogonality_atlas() -> int:
    """Check that character modulation translates every target delta."""
    for difference in GROUP:
        histogram = [0] * P
        for character in GROUP:
            histogram[dot(character, difference)] += 1
        if difference == ZERO:
            require(
                histogram == [GROUP_SIZE] + [0] * (P - 1),
                "translated target delta lost its mass",
            )
        else:
            require(
                histogram == [P] * P,
                "translated target delta stopped cancelling",
            )

    for multiplier in range(1, P):
        phase = (0, multiplier)
        require(
            {add(endpoint_target, phase) for endpoint_target in GROUP}
            == set(GROUP),
            "deep phase stopped translating the target group bijectively",
        )
    return 12 * GROUP_SIZE**2


def affine_masks(phase: tuple[int, int]) -> tuple[set, set, set]:
    """Translate the pure-other, pure-deep, and fork target loci by -p."""
    pure_other = {
        sub((other, 0), phase)
        for other in range(1, P)
    }
    pure_deep = {
        sub((0, deep), phase)
        for deep in range(1, P)
    }
    fork = {
        sub((other, deep), phase)
        for other in range(1, P)
        for deep in range(1, P)
    }
    return pure_other, pure_deep, fork


def mask_atlas() -> tuple[int, int, int, int]:
    """Check the affine partition and its sharp excluded inverse point."""
    rows = 0
    zero_in_deep = 0
    for multiplier in range(1, P):
        phase = (0, multiplier)
        inverse_point = negate(phase)
        pure_other, pure_deep, fork = affine_masks(phase)
        require(
            (len(pure_other), len(pure_deep), len(fork)) == (12, 12, 144),
            "affine support-mask sizes changed",
        )
        require(
            pure_other.isdisjoint(pure_deep)
            and pure_other.isdisjoint(fork)
            and pure_deep.isdisjoint(fork),
            "affine support masks overlap",
        )
        require(
            pure_other | pure_deep | fork == set(GROUP) - {inverse_point},
            "affine masks stopped partitioning the complement of -p",
        )
        require(inverse_point not in pure_other | pure_deep | fork, "-p entered a nonzero target mask")
        require(ZERO not in pure_other and ZERO not in fork, "zero entered the wrong affine mask")
        require(ZERO in pure_deep, "zero endpoint target stopped catalysing the deep axis")
        zero_in_deep += 1
        rows += 1
    return rows, zero_in_deep, 12, 144


def shift_coefficients(
    endpoint: dict[tuple[int, int], Fraction],
    phase: tuple[int, int],
) -> dict[tuple[int, int], Fraction]:
    """Apply A_H(q)=A_K(q-p)."""
    return {
        target: endpoint.get(sub(target, phase), Fraction(0))
        for target in GROUP
    }


def hostile_and_positive_controls() -> tuple[int, int, Fraction, Fraction]:
    """Check the inverse-character hostile and constant-response catalyst."""
    phase = (0, 5)
    inverse_point = negate(phase)

    inverse_character_endpoint = {inverse_point: Fraction(7, 3)}
    zero_only_full = shift_coefficients(inverse_character_endpoint, phase)
    require(
        zero_only_full[ZERO] == Fraction(7, 3)
        and sum(value != 0 for value in zero_only_full.values()) == 1,
        "inverse-character boundary stopped landing at zero",
    )

    constant_endpoint_response = {ZERO: Fraction(11, 4)}
    deep_only_full = shift_coefficients(constant_endpoint_response, phase)
    require(
        deep_only_full[phase] == Fraction(11, 4)
        and sum(value != 0 for value in deep_only_full.values()) == 1,
        "constant endpoint response stopped landing on the deep axis",
    )

    generic_endpoint = {
        point: Fraction(
            1 + point[0] + 2 * point[1],
            1 + ((3 * point[0] + point[1]) % 5),
        )
        for point in GROUP
    }
    full = shift_coefficients(generic_endpoint, phase)
    full_nonzero_energy = sum(
        value * value
        for target, value in full.items()
        if target != ZERO
    )
    projected_endpoint_energy = sum(
        value * value
        for point, value in generic_endpoint.items()
        if point != inverse_point
    )
    require(
        full_nonzero_energy == projected_endpoint_energy > 0,
        "inverse-character projection energy changed",
    )
    return (
        sum(value != 0 for value in zero_only_full.values()),
        sum(value != 0 for value in deep_only_full.values()),
        full_nonzero_energy,
        generic_endpoint[inverse_point],
    )


translation_checks = translation_orthogonality_atlas()
mask_rows, zero_catalyst_rows, pure_size, fork_size = mask_atlas()
(
    inverse_hostile_support,
    constant_catalyst_support,
    projection_energy,
    excluded_coefficient,
) = hostile_and_positive_controls()

require(translation_checks == 12 * GROUP_SIZE**2, "translation check census changed")
require(mask_rows == 12 == zero_catalyst_rows, "deep-axis multiplier census changed")
require(inverse_hostile_support == 1, "inverse hostile support changed")
require(constant_catalyst_support == 1, "constant catalyst support changed")

print("theorem=THM-2343")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"target_group_size={GROUP_SIZE}")
print(f"deep_axis_nonzero_phases={mask_rows}")
print(f"translation_orthogonality_checks={translation_checks}")
print("full_target_translation=A_H(q)=A_K(q-p)")
print("zero_only_full_iff=endpoint_response_is_inverse_deep_character")
print("constant_endpoint_response=full_current_supported_at_nonzero_deep_axis")
print(f"inverse_character_hostile_support={inverse_hostile_support}")
print(f"constant_response_catalyst_support={constant_catalyst_support}")
print(f"affine_pure_other_mask_size={pure_size}")
print(f"affine_pure_deep_mask_size={pure_size}")
print(f"affine_fork_mask_size={fork_size}")
print(f"zero_endpoint_target_in_pure_deep_masks={zero_catalyst_rows}")
print("affine_masks_partition=G_minus_inverse_deep_point")
print(f"projection_control_energy={projection_energy}")
print(f"projection_control_excluded_coefficient={excluded_coefficient}")
print("endpoint_inverse_character_exclusion=OPEN")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
