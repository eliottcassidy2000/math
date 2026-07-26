#!/usr/bin/env python3
"""Dependency-free exact companion for THM-2379.

The script checks the two strict-open thirteen-shift count identities on
their complete common rational-cell refinement.  It then exhausts every
anchored deep-safe / failed-safe-factor profile on that refinement,
including the complement sign, mixed-corner and energy bounds.  The
remaining checks are arithmetic: the 91-unit Abel multiplier gate, the
four-way blocker-status loss, and exact open neighborhoods in the
THM-2269 hostile speed family.

Every proof check uses ``require`` rather than ``assert``, so all checks
remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys
from typing import Iterable, Sequence


P = 13
NONZERO_COLOURS = (P - 1) ** 2
GUARD_AND_UNIT_SPEEDS = (
    ("H", 1, 2),
    ("q1", 2, 1),
    ("q2", 339, 1),
    ("q3", 677, 1),
    ("q4", 1015, 1),
    ("q5", 1353, 1),
    ("c1", 13, 1),
    ("c2", 169, 1),
    ("c3", 2197, 1),
)


def require(condition: bool, message: str) -> None:
    """Raise under both ordinary and optimized Python if a check fails."""
    if not condition:
        raise RuntimeError(message)


def fractional_part(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def circular_distance(value: Fraction) -> Fraction:
    reduced = fractional_part(value)
    return min(reduced, 1 - reduced)


def circle_distance(left: Fraction, right: Fraction) -> Fraction:
    difference = abs(fractional_part(left) - fractional_part(right))
    return min(difference, 1 - difference)


def danger(width: int, value: Fraction) -> int:
    """The strict-open indicator d_width(value)."""
    require(width in (1, 2), "only the unit and guard widths are in scope")
    return int(circular_distance(value) < Fraction(width, 14))


def shifted_danger(width: int, value: Fraction, shift: int) -> int:
    return danger(width, value - Fraction(shift, P))


def target_boundaries(width: int) -> set[Fraction]:
    """Boundaries of all value -> d_width(value-r/13)."""
    half_width = Fraction(width, 14)
    return {
        fractional_part(Fraction(shift, P) + sign * half_width)
        for shift in range(P)
        for sign in (-1, 1)
    }


def successor_boundaries(width: int) -> set[Fraction]:
    """Boundaries of value -> d_width(13 value)."""
    root_half_width = Fraction(width, 14 * P)
    return {
        fractional_part(Fraction(sheet, P) + sign * root_half_width)
        for sheet in range(P)
        for sign in (-1, 1)
    }


def cell_midpoints(boundaries: Iterable[Fraction]) -> tuple[Fraction, ...]:
    ordered = tuple(sorted(set(boundaries)))
    require(len(ordered) >= 2, "rational-cell boundary set is too small")
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        midpoints.append(fractional_part((left + right) / 2))
    return tuple(midpoints)


def mixed_corner(deep: Sequence[int], role: Sequence[int]) -> Fraction:
    """Sum of all normalized Fourier coefficients with both colours nonzero."""
    require(
        len(deep) == P and len(role) == P,
        "mixed profile stopped having thirteen shifts",
    )
    value_at_origin = deep[0] * role[0]
    origin_column = deep[0] * sum(role)
    origin_row = role[0] * sum(deep)
    total = sum(deep) * sum(role)
    return (
        value_at_origin
        - Fraction(origin_row + origin_column, P)
        + Fraction(total, P * P)
    )


def mixed_energy(deep: Sequence[int], role: Sequence[int]) -> Fraction:
    """Exact normalized Parseval energy in the fully mixed Fourier corner."""
    deep_total = sum(deep)
    role_total = sum(role)
    grand_mean = Fraction(deep_total * role_total, P * P)
    square_sum = Fraction(0)
    for row in range(P):
        row_mean = Fraction(deep[row] * role_total, P)
        for column in range(P):
            column_mean = Fraction(role[column] * deep_total, P)
            residual = (
                deep[row] * role[column]
                - row_mean
                - column_mean
                + grand_mean
            )
            square_sum += residual * residual
    return square_sum / (P * P)


def multiplier_is_live(width: int, multiplier: int) -> bool:
    """Whether sin(pi*width*multiplier/7) is nonzero, exactly."""
    require(multiplier != 0, "zero multiplier is outside a mixed colour")
    return (width * multiplier) % 7 != 0


def status_set(center: Fraction) -> frozenset[str]:
    return frozenset(
        name
        for name, speed, width in GUARD_AND_UNIT_SPEEDS
        if danger(width, speed * center)
    )


def boundary_radius(
    center: Fraction, speed: int, width: int
) -> Fraction:
    """Exact distance in x to the nearest strict-comb boundary."""
    image = fractional_part(speed * center)
    half_width = Fraction(width, 14)
    image_radius = min(
        circle_distance(image, half_width),
        circle_distance(image, 1 - half_width),
    )
    return image_radius / speed


def common_status_radius(
    center: Fraction,
) -> tuple[Fraction, tuple[str, ...]]:
    candidates = tuple(
        (boundary_radius(center, speed, width), name)
        for name, speed, width in GUARD_AND_UNIT_SPEEDS
    )
    radius = min(value for value, _ in candidates)
    limiting = tuple(name for value, name in candidates if value == radius)
    return radius, limiting


def check_hostile_neighborhood(
    label: str,
    center: Fraction,
    expected_status: frozenset[str],
    expected_radius: Fraction,
    expected_limiting: tuple[str, ...],
) -> None:
    require(status_set(center) == expected_status, f"{label}: center status changed")
    radius, limiting = common_status_radius(center)
    require(radius == expected_radius, f"{label}: open radius changed")
    require(limiting == expected_limiting, f"{label}: limiting comb changed")

    for sign in (-1, 1):
        require(
            status_set(center + sign * radius / 2) == expected_status,
            f"{label}: status changed inside the certified open radius",
        )

    boundary_hits = []
    for sign in (-1, 1):
        endpoint = center + sign * radius
        for name, speed, width in GUARD_AND_UNIT_SPEEDS:
            if name not in limiting:
                continue
            boundary_hits.append(
                circular_distance(speed * endpoint) == Fraction(width, 14)
            )
    require(any(boundary_hits), f"{label}: claimed radius reaches no comb boundary")


def main() -> None:
    # Keep the stored exact transcript byte-stable on Windows as well as
    # POSIX; theorem frontmatter hashes the LF working-tree bytes.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # The target-shift and successor boundaries are the same 26 rational
    # points for each width.  Their width-1/width-2 union is the complete
    # 52-cell refinement on which every indicator in both identities is
    # constant.
    boundary_sets = {}
    for width in (1, 2):
        target = target_boundaries(width)
        successor = successor_boundaries(width)
        require(len(target) == 26, f"L={width}: target boundary count changed")
        require(
            target == successor,
            f"L={width}: target and successor boundaries diverged",
        )
        boundary_sets[width] = target
    require(
        boundary_sets[1].isdisjoint(boundary_sets[2]),
        "unit and guard boundary refinements unexpectedly overlap",
    )
    samples = cell_midpoints(boundary_sets[1] | boundary_sets[2])
    require(len(samples) == 52, "combined rational refinement changed")

    for value in samples:
        for width in (1, 2):
            shifted_count = sum(
                shifted_danger(width, value, shift) for shift in range(P)
            )
            successor_count = 2 * width - danger(width, P * value)
            require(
                shifted_count == successor_count,
                f"L={width}: exact shift identity failed at {value}",
            )
            complement_count = sum(
                1 - shifted_danger(width, value, shift)
                for shift in range(P)
            )
            require(
                complement_count
                == P - 2 * width + danger(width, P * value),
                f"L={width}: complement shift identity failed at {value}",
            )

    # Exhaust every pair of common open cells satisfying the two anchor
    # hypotheses d_1(y)=0 and d_L(z)=1.
    profile_rows = {}
    distinct_profile_rows = {}
    corner_minima = {}
    energy_minima = {}
    for width in (1, 2):
        profile_count = 0
        distinct_profiles = set()
        corner_magnitudes = []
        energies = []
        uniform_corner = Fraction(P - 2 * width, P * P)
        uniform_energy = Fraction(
            (P - 2 * width) ** 2,
            P**4 * NONZERO_COLOURS,
        )
        for deep_value in samples:
            if danger(1, deep_value):
                continue
            deep = tuple(
                shifted_danger(1, deep_value, shift) for shift in range(P)
            )
            require(deep[0] == 0, "deep anchored-zero hypothesis changed")
            require(
                sum(deep) == 2 - danger(1, P * deep_value),
                "deep shift count changed inside corner exhaustion",
            )
            for role_value in samples:
                if not danger(width, role_value):
                    continue
                role = tuple(
                    shifted_danger(width, role_value, shift)
                    for shift in range(P)
                )
                complement = tuple(1 - bit for bit in role)
                require(role[0] == 1, "failed-safe-factor anchor changed")
                require(complement[0] == 0, "repair complement lost anchored zero")
                require(
                    all(
                        deep[row] * role[column]
                        + deep[row] * complement[column]
                        == deep[row]
                        for row in range(P)
                        for column in range(P)
                    ),
                    "danger/complement table stopped being shift-independent",
                )

                successor_deep = danger(1, P * deep_value)
                successor_role = danger(width, P * role_value)
                expected_magnitude = Fraction(
                    (2 - successor_deep)
                    * (P - 2 * width + successor_role),
                    P * P,
                )
                minus_corner = mixed_corner(deep, role)
                plus_corner = mixed_corner(deep, complement)
                require(
                    minus_corner == -expected_magnitude,
                    "negative anchored corner identity changed",
                )
                require(
                    plus_corner == expected_magnitude,
                    "complement anchored corner identity changed",
                )
                require(
                    plus_corner == -minus_corner,
                    "mixed complement sign algebra changed",
                )
                require(
                    all(deep[row] * complement[0] == 0 for row in range(P)),
                    "repair table stopped vanishing on its original role slice",
                )
                require(
                    expected_magnitude >= uniform_corner,
                    "uniform corner floor failed",
                )

                minus_energy = mixed_energy(deep, role)
                plus_energy = mixed_energy(deep, complement)
                require(
                    minus_energy == plus_energy,
                    "complement changed the fully mixed energy",
                )
                require(
                    minus_energy >= minus_corner**2 / NONZERO_COLOURS,
                    "144-colour Cauchy floor failed",
                )
                require(
                    minus_energy >= uniform_energy,
                    "uniform mixed-energy floor failed",
                )

                profile_count += 1
                distinct_profiles.add((deep, role))
                corner_magnitudes.append(expected_magnitude)
                energies.append(minus_energy)

        expected_profiles = 315 if width == 1 else 585
        require(
            profile_count == expected_profiles,
            f"L={width}: anchored profile count changed",
        )
        expected_distinct_profiles = 69 if width == 1 else 161
        require(
            len(distinct_profiles) == expected_distinct_profiles,
            f"L={width}: distinct anchored bit-profile count changed",
        )
        require(
            min(corner_magnitudes) == uniform_corner,
            f"L={width}: exact corner minimum changed",
        )
        profile_rows[width] = profile_count
        distinct_profile_rows[width] = len(distinct_profiles)
        corner_minima[width] = min(corner_magnitudes)
        energy_minima[width] = min(energies)

    unit_coefficient = Fraction(11, P * P * NONZERO_COLOURS)
    guard_coefficient = Fraction(9, P * P * NONZERO_COLOURS)
    require(unit_coefficient == Fraction(11, 24336), "unit coefficient changed")
    require(guard_coefficient == Fraction(1, 2704), "guard coefficient changed")
    unit_energy = Fraction(11**2, P**4 * NONZERO_COLOURS)
    guard_energy = Fraction(9**2, P**4 * NONZERO_COLOURS)

    # A representative bank with a symmetric lift-index range for every
    # nonzero mod-13
    # character.  The exact sine-zero rule, rather than floating point,
    # decides whether a probe multiplier is live.
    lift_range = tuple(range(-13, 14))
    raw_multipliers = tuple(
        residue + P * lift
        for residue in range(1, P)
        for lift in lift_range
    )
    require(len(raw_multipliers) == 324, "representative lift bank changed")
    live_multipliers = tuple(
        multiplier
        for multiplier in raw_multipliers
        if multiplier_is_live(1, multiplier)
    )
    killed_multipliers = tuple(
        multiplier
        for multiplier in raw_multipliers
        if not multiplier_is_live(1, multiplier)
    )
    require(
        len(live_multipliers) == 277 and len(killed_multipliers) == 47,
        "representative 7-divisibility counts changed",
    )
    require(
        all(gcd(abs(multiplier), 91) == 1 for multiplier in live_multipliers),
        "a live probe multiplier stopped being a 91-unit",
    )
    require(
        all(
            multiplier % 7 == 0 and multiplier % P != 0
            for multiplier in killed_multipliers
        ),
        "a killed probe multiplier has the wrong arithmetic type",
    )
    require(
        all(
            multiplier_is_live(1, multiplier)
            == multiplier_is_live(2, multiplier)
            for multiplier in raw_multipliers
        ),
        "unit and guard probe zero sets diverged",
    )
    live_pairs = 0
    for first in live_multipliers:
        for second in live_multipliers:
            require(
                gcd(abs(first * second), 91) == 1,
                "a live multiplier pair stopped being a 91-unit pair",
            )
            live_pairs += 1
    require(live_pairs == 76729, "representative live-pair count changed")

    # THM-2370 gives rho >= delta/n before any blocker-status split.
    # Demanding one literal status among the two nondeep blocker bits
    # incurs exactly the worst-case four-way pigeonhole loss.
    labelled_unit = unit_coefficient / 4
    labelled_guard = guard_coefficient / 4
    require(
        labelled_unit == Fraction(11, 97344),
        "four-way labelled unit constant changed",
    )
    require(
        labelled_guard == Fraction(1, 10816),
        "four-way labelled guard constant changed",
    )
    labelled_unit_energy = unit_energy / 16
    labelled_guard_energy = guard_energy / 16
    require(
        labelled_unit_energy
        == Fraction(121, 16 * P**4 * NONZERO_COLOURS),
        "four-way labelled unit energy changed",
    )
    require(
        labelled_guard_energy
        == Fraction(81, 16 * P**4 * NONZERO_COLOURS),
        "four-way labelled guard energy changed",
    )

    # Exact local hostiles in the THM-2269 speed family.  The first two
    # centers show that the c1/c2 blocker-status word may be empty in both
    # the guard-deletion and unit-deletion branches.  The latter two retain
    # the canonical branch center and its first image as positive controls.
    hostile_rows = (
        (
            "empty-guard",
            Fraction(1, 8),
            frozenset({"H"}),
            Fraction(1, 25256),
            ("q5",),
        ),
        (
            "empty-unit",
            Fraction(1, 2),
            frozenset({"q1"}),
            Fraction(3, 15379),
            ("c3",),
        ),
        (
            "thm2269-x0",
            Fraction(79, 338),
            frozenset({"c1"}),
            Fraction(64, 533533),
            ("q5",),
        ),
        (
            "thm2269-x1",
            Fraction(1, 26),
            frozenset({"H", "q2", "q3", "q4", "q5"}),
            Fraction(1, 41041),
            ("q5",),
        ),
    )
    for hostile in hostile_rows:
        check_hostile_neighborhood(*hostile)

    print("THM-2379 guard/unit factor-repair exact companion")
    print(
        "rational refinement: "
        f"L1 boundaries={len(boundary_sets[1])}; "
        f"L2 boundaries={len(boundary_sets[2])}; "
        f"combined cells={len(samples)}"
    )
    for width in (1, 2):
        print(
            f"L={width}: anchored cell pairs={profile_rows[width]}; "
            f"distinct bit-profile pairs={distinct_profile_rows[width]}; "
            f"corner_min={corner_minima[width]}; "
            f"actual_cell_energy_min={energy_minima[width]}"
        )
    print(
        "coefficient floors before status split: "
        f"unit={unit_coefficient}; guard={guard_coefficient}"
    )
    print(
        "energy floors before status split: "
        f"unit={unit_energy}; guard={guard_energy}"
    )
    print(
        "four-way labelled floors: "
        f"unit={labelled_unit}; guard={labelled_guard}"
    )
    print(
        "four-way labelled energy floors: "
        f"unit={labelled_unit_energy}; guard={labelled_guard_energy}"
    )
    print(
        "91-unit bank: "
        f"raw={len(raw_multipliers)}; live={len(live_multipliers)}; "
        f"septimal_zero={len(killed_multipliers)}; live_pairs={live_pairs}"
    )
    for label, center, expected, radius, limiting in hostile_rows:
        print(
            f"{label}: center={center}; danger={','.join(sorted(expected))}; "
            f"open_radius={radius}; first_boundary={','.join(limiting)}"
        )
    print("VERDICT: anchored failed-safe mass forces a mixed complement probe")


if __name__ == "__main__":
    main()
