#!/usr/bin/env python3
"""Canonical endpoint-toggle audit of THM-4231's hybrid Q=1290 quadrant."""

from fractions import Fraction
from hashlib import sha256
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

POOL_DENOMINATOR = 18_241_159_416_480
BASE_TICKS = 4_579_301_272_924
BASE_COMPONENTS = 618
B_STAR = (170, 176, 190, 193, 240, 252, 264, 286, 290)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def toggle_geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    denominator = 1
    for speed in speeds:
        denominator = lcm(denominator, 14 * speed)

    change: dict[int, int] = {}
    for speed in speeds:
        unit = denominator // (14 * speed)
        for tooth in range(speed):
            leave_danger = (14 * tooth + 1) * unit
            enter_danger = (14 * tooth + 13) * unit
            change[leave_danger] = change.get(leave_danger, 0) - 1
            change[enter_danger] = change.get(enter_danger, 0) + 1

    failed = len(speeds)
    previous_wall = 0
    previous_safe = False
    mass = 0
    components = 0
    for wall in sorted(change):
        safe = failed == 0
        if safe:
            mass += wall - previous_wall
            if not previous_safe:
                components += 1
        previous_safe = safe
        failed += change[wall]
        require(failed >= 0, "negative failure multiplicity")
        previous_wall = wall
    safe = failed == 0
    if safe:
        mass += denominator - previous_wall
        if not previous_safe:
            components += 1
    require(failed == len(speeds), "endpoint state did not close")
    return denominator, mass, components


def analytic_gap(q: int, r: int) -> Fraction:
    mass = Fraction(BASE_TICKS, POOL_DENOMINATOR)
    base_surplus = (45 * mass - 4) / 63
    error_scale = Fraction(6 * BASE_COMPONENTS, 49)
    return base_surplus - error_scale * (Fraction(1, q) + Fraction(1, r))


def main() -> None:
    r0_rows: list[str] = []
    pair_rows: list[str] = []
    minimum: tuple[Fraction, int, int, int, int, int] | None = None
    minimum_count = 0
    maximum_denominator = 0

    for q in range(1290, 1307):
        r0 = q + 1
        while analytic_gap(q, r0) < 0:
            r0 += 1
        require(r0 == 2614 - q, f"unexpected analytic boundary at q={q}")
        previous_gap = analytic_gap(q, r0 - 1)
        boundary_gap = analytic_gap(q, r0)
        require(previous_gap < 0 <= boundary_gap, f"bad r0 bracket at q={q}")
        r0_rows.append(
            f"{q},{r0},{fraction_text(previous_gap)},{fraction_text(boundary_gap)}"
        )

        for r in range(q + 1, r0):
            denominator, mass, components = toggle_geometry(B_STAR + (q, r))
            delta_ticks = 63 * mass - 4 * denominator
            require(delta_ticks > 0, f"literal failure at {(q, r)}")
            delta = Fraction(delta_ticks, denominator)
            pair_rows.append(f"{q},{r},{denominator},{mass},{components},{delta_ticks}")
            maximum_denominator = max(maximum_denominator, denominator)
            candidate = (delta, q, r, denominator, mass, components)
            if minimum is None or candidate[0] < minimum[0]:
                minimum = candidate
                minimum_count = 1
            elif candidate[0] == minimum[0]:
                minimum_count += 1

    require(len(r0_rows) == 17, "r0 ledger size")
    require(len(pair_rows) == 289, "literal pair ledger size")
    require(minimum is not None, "missing minimum")
    delta, q_min, r_min, denominator_min, mass_min, components_min = minimum
    require(
        (q_min, r_min, components_min, delta) ==
        (1300, 1305, 1044, Fraction(18_886_235_531, 2_585_198_330)),
        "minimum literal margin changed",
    )
    require(minimum_count == 1, "minimum is not unique")

    r0_payload = ("\n".join(r0_rows) + "\n").encode("ascii")
    pair_payload = ("\n".join(pair_rows) + "\n").encode("ascii")

    print("LRC14_HYBRID_Q1290_PRIMARY_TOGGLE")
    print(
        "DIRECT_CENSUS MAX_K 1307 MAX_COUNT 1 SECOND_K 1290 "
        "SECOND_COUNT 1 EXCEPTION B_STAR"
    )
    for row in r0_rows:
        q, r0, previous, boundary = row.split(",")
        print(f"R0 Q {q} FIRST_DIRECT_R {r0} PREVIOUS_GAP {previous} AT_GAP {boundary}")
    print(
        f"R0_LEDGER COUNT {len(r0_rows)} FORMULA 2614-Q SHA256 "
        f"{sha256(r0_payload).hexdigest()}"
    )
    print(
        f"LITERAL_LEDGER PAIRS {len(pair_rows)} STRICT_SAFE {len(pair_rows)} "
        f"NONPOSITIVE 0 SHA256 {sha256(pair_payload).hexdigest()}"
    )
    print(
        f"MIN_LITERAL Q {q_min} R {r_min} D {denominator_min} MASS {mass_min} "
        f"COMPONENTS {components_min} DELTA_TICKS {63 * mass_min - 4 * denominator_min} "
        f"DELTA {fraction_text(delta)} COUNT {minimum_count}"
    )
    print(f"MAX_DENOMINATOR {maximum_denominator}")
    print(
        "VERDICT ALL_DISTINCT_Q_R_GE_1290_SAFE_RELATIVE_TO_DIRECT_CENSUS "
        "Q1290_NOT_CLAIMED_LITERAL_MINIMAL LRC14_OPEN"
    )


if __name__ == "__main__":
    main()
