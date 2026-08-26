#!/usr/bin/env python3
"""Exact endpoint-toggle boundary controls for THM-4231."""

from fractions import Fraction
from math import gcd


POOL_DENOMINATOR = 18_241_159_416_480
EXTREMAL_BODY = (170, 176, 190, 193, 240, 252, 264, 286, 290)
OLD_COVER = (85, 88, 143, 168, 193, 240, 252, 264, 290)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def ceil_fraction(value: Fraction) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def toggle_geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    common = 1
    for speed in speeds:
        common = lcm(common, 14 * speed)

    delta_at: dict[int, int] = {}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            leave = (14 * tooth + 1) * unit
            enter = (14 * tooth + 13) * unit
            delta_at[leave] = delta_at.get(leave, 0) - 1
            delta_at[enter] = delta_at.get(enter, 0) + 1

    failed = len(speeds)
    previous = 0
    previous_safe = False
    mass = 0
    components = 0
    for tick in sorted(delta_at):
        safe = failed == 0
        if safe:
            mass += tick - previous
            if not previous_safe:
                components += 1
        previous_safe = safe
        failed += delta_at[tick]
        require(failed >= 0, "negative failure count")
        previous = tick
    safe = failed == 0
    if safe:
        mass += common - previous
        if not previous_safe:
            components += 1
    require(failed == len(speeds), "toggle sweep did not close")
    return common, mass, components


def direct_activation(mass: Fraction, components: int) -> int:
    surplus = 45 * mass - 4
    require(surplus > 0, "nonpositive base surplus")
    return ceil_fraction(Fraction(108 * components, 7) / surplus)


def analytic_gap(
    mass: Fraction, components: int, reciprocal_sum: Fraction
) -> Fraction:
    return (45 * mass - 4) / 63 - Fraction(6 * components, 49) * reciprocal_sum


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    extremal_common, extremal_ticks, extremal_components = toggle_geometry(
        EXTREMAL_BODY
    )
    extremal_mass = Fraction(extremal_ticks, extremal_common)
    require(
        extremal_mass == Fraction(4_579_301_272_924, POOL_DENOMINATOR)
        and extremal_components == 618
        and direct_activation(extremal_mass, extremal_components) == 1_307,
        "extremal base geometry changed",
    )

    cover_common, cover_ticks, cover_components = toggle_geometry(OLD_COVER)
    cover_mass = Fraction(cover_ticks, cover_common)
    require(
        cover_mass == Fraction(4_802_564_195_362, POOL_DENOMINATOR)
        and cover_components == 506
        and direct_activation(cover_mass, cover_components) == 995,
        "old-cover base geometry changed",
    )

    symmetric_1306 = analytic_gap(
        extremal_mass, extremal_components, Fraction(2, 1_306)
    )
    symmetric_1307 = analytic_gap(
        extremal_mass, extremal_components, Fraction(2, 1_307)
    )
    distinct_1306 = analytic_gap(
        extremal_mass,
        extremal_components,
        Fraction(1, 1_306) + Fraction(1, 1_307),
    )
    distinct_1307 = analytic_gap(
        extremal_mass,
        extremal_components,
        Fraction(1, 1_307) + Fraction(1, 1_308),
    )
    bridge_1306 = analytic_gap(
        extremal_mass,
        extremal_components,
        Fraction(1, 1_306) + Fraction(1, 1_308),
    )
    require(symmetric_1306 < 0 < symmetric_1307, "symmetric boundary changed")
    require(distinct_1306 < 0 < distinct_1307, "distinct-pair boundary changed")
    require(bridge_1306 > 0, "1306 cofinal bridge ceased to be positive")

    literal_rows = []
    for q, r in ((1_306, 1_307), (1_307, 1_308)):
        common, ticks, components = toggle_geometry(EXTREMAL_BODY + (q, r))
        mass = Fraction(ticks, common)
        delta = 63 * mass - 4
        require(delta > 0, "literal boundary control ceased to be safe")
        literal_rows.append((q, r, components, mass, delta))

    print("LRC14_COMPOSITE_DESCENT_LITERAL_CONTROLS_THM4231")
    print(
        f"EXTREMAL_BASE MASS {fraction_text(extremal_mass)} "
        f"COMPONENTS {extremal_components} DIRECT_K 1307"
    )
    print(
        f"OLD_COVER_BASE MASS {fraction_text(cover_mass)} "
        f"COMPONENTS {cover_components} DIRECT_K 995"
    )
    print(
        f"SYMMETRIC_GAPS Q1306 {fraction_text(symmetric_1306)} "
        f"Q1307 {fraction_text(symmetric_1307)} SIGNS NEGATIVE,POSITIVE"
    )
    print(
        f"DISTINCT_GAPS PAIR1306,1307 {fraction_text(distinct_1306)} "
        f"PAIR1307,1308 {fraction_text(distinct_1307)} SIGNS NEGATIVE,POSITIVE"
    )
    print(
        f"COFINAL_BRIDGE PAIR1306,1308 {fraction_text(bridge_1306)} POSITIVE YES"
    )
    for q, r, components, mass, delta in literal_rows:
        print(
            f"LITERAL_PAIR {q},{r} COMPONENTS {components} "
            f"MASS {fraction_text(mass)} DELTA63 {fraction_text(delta)} POSITIVE YES"
        )
    print("VERDICT CERTIFICATE_BOUNDARY_NOT_LITERAL_COUNTEREXAMPLE")


if __name__ == "__main__":
    main()
