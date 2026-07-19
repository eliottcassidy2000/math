#!/usr/bin/env python3
"""Exact audit for THM-1233, all-prefix component-span compactification.

The proof is analytic.  This dependency-free replay verifies every rational
constant and ceiling recurrence, and independently merges exact tooth
intervals on a finite bank to audit the suffix-component span lemmas.

Tournament audit: speed order is a transitive loss quotient.  The faithful
vertices are suffix-union components rooted at their slowest tooth, with
prefix-survivor components as occupants.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    """Optimization-safe certificate check (unlike Python's ``assert``)."""
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: F) -> int:
    return -((-x.numerator) // x.denominator)


def tooth_intervals(speed: int, left: F, right: F) -> list[tuple[F, F]]:
    """Open teeth of one comb which meet a finite real-lift window."""
    first = floor_fraction(speed * left) - 2
    last = ceil_fraction(speed * right) + 2
    result: list[tuple[F, F]] = []
    for index in range(first, last + 1):
        lo = F(14 * index - 1, 14 * speed)
        hi = F(14 * index + 1, 14 * speed)
        if lo < right and left < hi:
            result.append((max(lo, left), min(hi, right)))
    return result


def union_components(
    speeds: tuple[int, ...], left: F = F(-1), right: F = F(2)
) -> list[tuple[F, F]]:
    """Components of a finite window of a union of strict open teeth.

    Intervals which only abut are not merged; a later interval can merge them
    only if it has positive overlap with both, exactly as for the open union.
    """
    intervals = sorted(
        interval
        for speed in speeds
        for interval in tooth_intervals(speed, left, right)
    )
    components: list[tuple[F, F]] = []
    for lo, hi in intervals:
        if not components or lo >= components[-1][1]:
            components.append((lo, hi))
        else:
            old_lo, old_hi = components[-1]
            components[-1] = (old_lo, max(old_hi, hi))
    return components


def component_span_bank() -> tuple[int, int, int, int]:
    generic = one = two = three = 0
    for size in range(1, 7):
        for speeds in combinations(range(1, 11), size):
            slow = speeds[0]
            components = union_components(speeds)
            for lo, hi in components:
                span = hi - lo
                # Clipping at the audit-window boundary can only shorten.
                require(span < F(13, 7 * slow), "generic span bound failed")
                generic += 1
                if size == 1:
                    require(span <= F(1, 7 * slow), "one-comb span bound failed")
                    one += 1
                elif size == 2:
                    require(
                        span < F(1, 7 * speeds[0]) + F(2, 7 * speeds[1]),
                        "two-comb exact span bound failed",
                    )
                    require(span < F(3, 7 * speeds[0]), "two-comb coarse span failed")
                    two += 1
                elif size == 3:
                    require(
                        span
                        < F(1, 7 * speeds[0])
                        + F(2, 7 * speeds[1])
                        + F(4, 7 * speeds[2]),
                        "three-comb exact span bound failed",
                    )
                    require(span < F(1, speeds[0]), "three-comb coarse span failed")
                    three += 1
    return generic, one, two, three


def constant_ledger() -> tuple[
    dict[int, F], dict[int, F], dict[int, F], dict[int, F]
]:
    f_mass: dict[int, F] = {}
    normalized_length: dict[int, F] = {}
    physical_coefficient: dict[int, F] = {}
    ratio_coefficient: dict[int, F] = {}

    # Suffix-component span is span_factor[j]/d_(j+1).
    span_factor = {
        1: F(13, 7),
        2: F(13, 7),
        3: F(1),
        4: F(3, 7),
        5: F(1, 7),
    }
    expected_ratio = {
        1: F(91, 29),
        2: F(91, 22),
        3: F(49, 15),
        4: F(21, 8),
        5: F(7),
    }

    for j in range(1, 6):
        f_mass[j] = 1 - j * F(7, 36)
        normalized_length[j] = f_mass[j] / F(7, 6)
        physical_coefficient[j] = normalized_length[j] * F(6, 7)
        require(f_mass[j] == F(36 - 7 * j, 36), "f-mass ledger failed")
        require(
            normalized_length[j] == F(36 - 7 * j, 42),
            "normalized survivor ledger failed",
        )
        require(
            physical_coefficient[j] == F(36 - 7 * j, 49),
            "physical survivor ledger failed",
        )
        ratio_coefficient[j] = span_factor[j] / physical_coefficient[j]
        require(ratio_coefficient[j] == expected_ratio[j], "ratio ledger failed")

    return f_mass, normalized_length, physical_coefficient, ratio_coefficient


def ceiling_and_ratio_audit(
    ratio_coefficient: dict[int, F]
) -> tuple[dict[int, F], list[int], list[F], int]:
    component_adjacent_expected = {
        1: F(273, 29),
        2: F(455, 22),
        3: F(343, 15),
        4: F(189, 8),
        5: F(77),
    }
    component_adjacent = {
        j: ratio_coefficient[j] * (2 * j + 1) for j in range(1, 6)
    }
    require(component_adjacent == component_adjacent_expected, "adjacent ledger failed")
    adjacent = dict(component_adjacent)
    # THM-1232's functional gate d3/c<84/5 and d2/c>1 sharpen the
    # component-only d3/d2<455/22 row.
    adjacent[2] = F(84, 5)

    ceiling_checks = 0
    for c in range(1, 101):
        for d in range(c + 1, 20 * c + 1):
            require(
                ceil_fraction(F(6 * d + c, 7 * c)) <= ceil_fraction(F(d, c)),
                "ceiling monotonicity failed",
            )
            ceiling_checks += 1

    def tooth_count_cap(strict_upper: F) -> int:
        return ceil_fraction((6 * strict_upper + 1) / 7)

    # Direct recursive projective caps from the theorem and THM-1232's
    # independent functional gate at the third speed.
    x1 = F(13, 6)
    cap1 = tooth_count_cap(x1)
    require(cap1 == 2, "first tooth cap failed")
    x2 = ratio_coefficient[1] * (1 + cap1)
    require(x2 == F(273, 29), "second ratio cap failed")
    cap2 = tooth_count_cap(x2)
    require(cap2 == 9, "second tooth cap failed")

    component_x3 = ratio_coefficient[2] * (1 + 2 * cap2)
    require(component_x3 == F(1729, 22), "component third ratio cap failed")
    x3 = F(84, 5)
    require(x3 < component_x3, "functional third ratio is not sharper")
    cap3 = tooth_count_cap(x3)
    require(cap3 == 15, "third tooth cap failed")

    x4 = ratio_coefficient[3] * (1 + cap1 + cap2 + cap3)
    require(x4 == F(441, 5), "fourth ratio cap failed")
    cap4 = tooth_count_cap(x4)
    require(cap4 == 76, "fourth tooth cap failed")

    x5 = ratio_coefficient[4] * (1 + cap1 + cap2 + cap3 + cap4)
    require(x5 == F(2163, 8), "fifth ratio cap failed")
    cap5 = tooth_count_cap(x5)
    require(cap5 == 232, "fifth tooth cap failed")

    x6 = ratio_coefficient[5] * (1 + cap1 + cap2 + cap3 + cap4 + cap5)
    require(x6 == 2345, "sixth ratio cap failed")

    return (
        adjacent,
        [cap1, cap2, cap3, cap4, cap5],
        [x2, x3, x4, x5, x6],
        ceiling_checks,
    )


def tournament_audit() -> None:
    speeds = (7, 9, 13, 20, 31, 48)
    scores = tuple(sum(speed > other for other in speeds) for speed in speeds)
    require(scores == (0, 1, 2, 3, 4, 5), "tournament score audit failed")
    print("TOURNAMENT_LOSS_AUDIT")
    print("observable=log(d_j/d_i), gauge=positive sign")
    print("scores=0,1,2,3,4,5 cycles=0 SCCs=1,1,1,1,1,1")
    print("hamiltonian_paths=1 tie_path=7,9,13,20,31,48")
    print("faithful_vertices=suffix-union-components rooted at slow teeth")


def main() -> None:
    print("THM-1233 ALL-PREFIX COMPONENT-SPAN COMPACTIFICATION")
    f_mass, normalized, physical, ratio = constant_ledger()
    for j in range(1, 6):
        print(
            f"j={j} f_mass={f_mass[j]} normalized_survivor={normalized[j]} "
            f"physical_survivor={physical[j]}/c ratio_coefficient={ratio[j]}"
        )

    adjacent, caps, absolute_caps, ceiling_checks = ceiling_and_ratio_audit(ratio)
    print("ADJACENT_RATIO_CHAIN")
    print("d1/c<13/6")
    for j in range(1, 6):
        print(f"d{j + 1}/d{j}<{adjacent[j]}")

    print("DIRECT_PROJECTIVE_RECURSION")
    print("ceil_caps_d1_to_d5=" + ",".join(map(str, caps)))
    print("strict_ratio_caps_d2_to_d6=" + ",".join(map(str, absolute_caps)))
    print("d2/c<273/29 functional_d3/c<84/5")
    print("d4/c<441/5 d5/c<2163/8 d6/c<2345")
    print(f"ceiling_monotonicity_checks={ceiling_checks}")

    generic, one, two, three = component_span_bank()
    print("FINITE_EXACT_COMPONENT_TELEMETRY_NOT_PROOF")
    print(f"generic_component_checks={generic}")
    print(f"one_comb_checks={one} two_comb_checks={two} three_comb_checks={three}")
    print("generic_span<13/(7*s1), two_span<1/(7*s1)+2/(7*s2)")
    print("three_span<1/(7*s1)+2/(7*s2)+4/(7*s3)<1/s1")

    tournament_audit()
    print("PROVED=all-prefix survivor/component composition and projective compactness")
    print("OPEN=unbounded carrier-phase tangent stalk and compact-box exclusion")
    print("VERIFIED_EXACT")


if __name__ == "__main__":
    main()
