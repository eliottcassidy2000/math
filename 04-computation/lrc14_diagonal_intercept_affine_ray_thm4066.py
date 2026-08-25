#!/usr/bin/env python3
"""Exact wall-cell certificate for THM-4066.

The primary path constructs every strict fully-spoiled phase component from
the literal labelled lift masks.  It does not import the affine-circuit
implementations of THM-4030/4032.
"""

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def masks(
    y: F,
    d: int,
    exceptions: tuple[int, ...],
    closed: bool = False,
) -> tuple[tuple[int, ...], ...]:
    compare = (lambda value: value <= F(1, 14)) if closed else (
        lambda value: value < F(1, 14)
    )
    return tuple(
        tuple(
            label
            for label in range(d)
            if compare(distance(speed, (y + label) / d))
        )
        for speed in exceptions
    )


def fully_spoiled(y: F, d: int, exceptions: tuple[int, ...]) -> bool:
    return set().union(*map(set, masks(y, d, exceptions))) == set(range(d))


def wall_components(
    d: int,
    exceptions: tuple[int, ...],
) -> tuple[tuple[F, F], ...]:
    """Return the strict fully-spoiled components in the gauge 0<=y<=1."""
    walls = {F(0), F(1)}
    for speed in exceptions:
        for label in range(d):
            for integer in range(-2, speed + 3):
                for sign in (-1, 1):
                    wall = F(d, speed) * (
                        F(integer) + sign * F(1, 14)
                    ) - label
                    if 0 <= wall <= 1:
                        walls.add(wall)
    ordered = sorted(walls)
    components: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        if not fully_spoiled((left + right) / 2, d, exceptions):
            continue
        if (
            components
            and components[-1][1] == left
            and fully_spoiled(left, d, exceptions)
        ):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return tuple(components)


def inverse_components(
    components: tuple[tuple[F, F], ...],
    dilation: int,
) -> tuple[tuple[F, F], ...]:
    return tuple(
        sorted(
            (
                (left + sheet) / dilation,
                (right + sheet) / dilation,
            )
            for left, right in components
            for sheet in range(dilation)
        )
    )


def check_template(
    d: int,
    primitive: tuple[int, ...],
    expected: tuple[tuple[F, F], ...],
    dilations: tuple[int, ...],
) -> int:
    require(wall_components(d, primitive) == expected, f"base window d={d}")
    phase_checks = 0
    for dilation in dilations:
        require(gcd(dilation, d) == 1, "unit dilation")
        scaled = tuple(dilation * speed for speed in primitive)
        require(
            wall_components(d, scaled) == inverse_components(expected, dilation),
            f"component pullback d={d},h={dilation}",
        )
        for denominator in range(1, 80):
            for numerator in range(denominator):
                y = F(numerator, denominator)
                base_masks = masks(dilation * y, d, primitive)
                scaled_masks = masks(y, d, scaled)
                label_image = tuple(
                    tuple(sorted((dilation * label) % d for label in mask))
                    for mask in scaled_masks
                )
                require(
                    label_image == base_masks,
                    f"label conjugacy d={d},h={dilation},y={y}",
                )
                phase_checks += 1
    return phase_checks


def clearance(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed, phase) for speed in speeds)


def main() -> None:
    d2_base = ((F(13, 35), F(8, 21)), (F(13, 21), F(22, 35)))
    d3_base = ((F(11, 56), F(3, 14)), (F(11, 14), F(45, 56)))
    d4_base = ((F(5, 49), F(1, 7)), (F(6, 7), F(44, 49)))

    phase_checks = 0
    phase_checks += check_template(2, (3, 5), d2_base, (1, 3, 5, 7, 9, 11))
    phase_checks += check_template(
        3,
        (1, 4, 5),
        d3_base,
        (1, 2, 4, 5, 7, 8, 10, 11),
    )
    phase_checks += check_template(4, (2, 7, 9), d4_base, (1, 3, 5, 7, 9, 11))

    require(
        masks(F(23, 112), 3, (1, 4, 5)) == ((0,), (2,), (1,)),
        "d3 first mask word",
    )
    require(
        masks(F(89, 112), 3, (1, 4, 5)) == ((2,), (0,), (1,)),
        "d3 second mask word",
    )
    require(
        masks(F(6, 49), 4, (2, 7, 9)) == ((0, 2), (1,), (3,)),
        "d4 first mask word",
    )
    require(
        masks(F(43, 49), 4, (2, 7, 9)) == ((1, 3), (2,), (0,)),
        "d4 second mask word",
    )

    endpoint_rows = (
        (2, (3, 5), (F(13, 35), F(8, 21), F(13, 21), F(22, 35))),
        (3, (1, 4, 5), (F(11, 56), F(3, 14), F(11, 14), F(45, 56))),
        (4, (2, 7, 9), (F(5, 49), F(1, 7), F(6, 7), F(44, 49))),
    )
    endpoint_checks = 0
    for d, primitive, endpoints in endpoint_rows:
        for endpoint in endpoints:
            require(
                not fully_spoiled(endpoint, d, primitive),
                "strict endpoint must leak",
            )
            closed_union = set().union(
                *map(set, masks(endpoint, d, primitive, closed=True))
            )
            require(closed_union == set(range(d)), "closed endpoint must cover")
            endpoint_checks += 1

    # If h is not a unit modulo d, label multiplication is not a permutation.
    # The true nonunit spoiled measure is far larger than the false pullback.
    nonunit = wall_components(3, (3, 12, 15))
    nonunit_measure = sum((right - left for left, right in nonunit), F(0))
    false_pullback_measure = sum(
        (right - left for left, right in inverse_components(d3_base, 3)),
        F(0),
    )
    require(nonunit_measure == F(51, 140), "nonunit hostile measure")
    require(false_pullback_measure == F(1, 28), "false pullback measure")
    require(nonunit_measure != false_pullback_measure, "unit hypothesis is load-bearing")

    h10 = tuple(range(1, 11))
    d3_row = tuple(3 * speed for speed in h10) + (5, 20, 25)
    d4_row = tuple(4 * speed for speed in h10) + (22, 77, 99)
    require(len(d3_row) == len(set(d3_row)) == 13, "d3 runner/disjointness gate")
    require(len(d4_row) == len(set(d4_row)) == 13, "d4 runner/disjointness gate")
    require(all(speed % 3 for speed in (5, 20, 25)), "d3 exception types")
    require(22 % 4 == 2 and 77 % 2 and 99 % 2, "d4 exception types")
    require(24 * 5 >= 11 * 10 and 5 * 5 < 11 * 10, "d3 strict cone gain")
    require(21 * 11 >= 22 * 10 and 27 * 11 < 44 * 10, "d4 strict cone gain")
    require(clearance(d3_row, F(1, 33)) == F(1, 11), "d3 explicit witness")
    require(clearance(d4_row, F(1, 44)) == F(1, 11), "d4 explicit witness")

    print("LRC14_DIAGONAL_INTERCEPT_AFFINE_RAY_THM4066")
    print("scope=LRC14_OPEN;exact_affine_ray_families_only")
    print("pullback=h*D_hdelta(y)=D_delta(hy);gcd(h,d)=1")
    print("d2_control_components", d2_base, "width", F(1, 105), "measure", F(2, 105))
    print("d3_components", d3_base, "width", F(1, 56), "measure", F(1, 28))
    print("d3_scaled=2h_components_width_1/(56h);closure=24h>=11M")
    print("d4_components", d4_base, "width", F(2, 49), "measure", F(4, 49))
    print("d4_scaled=2h_components_width_2/(49h);closure=21h>=22M")
    print("endpoint_checks", endpoint_checks, "strict_leaks_closed_covers")
    print("nonunit_h3_measure", nonunit_measure, "false_pullback", false_pullback_measure)
    print("phase_label_conjugacy_checks", phase_checks)
    print("d3_gain_factor=24/5;example_witness=1/33;clearance=1/11")
    print("d4_gain_factor=14/9;example_witness=1/44;clearance=1/11")
    print("checks", CHECKS)
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
