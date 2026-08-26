#!/usr/bin/env python3
"""Independent midpoint-wall audit of THM-4231's hybrid Q=1290 quadrant."""

from fractions import Fraction
from hashlib import sha256
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

POOL_D = 18_241_159_416_480
BODY_MASS = 4_579_301_272_924
BODY_COMPONENTS = 618
EXTREMAL_BODY = (170, 176, 190, 193, 240, 252, 264, 286, 290)


def common_multiple(values: tuple[int, ...]) -> int:
    answer = 1
    for value in values:
        answer = answer // gcd(answer, value) * value
    return answer


def text_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def midpoint_geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    denominator = common_multiple(tuple(14 * speed for speed in speeds))
    walls = {0, denominator}
    for speed in speeds:
        unit = denominator // (14 * speed)
        walls.update((14 * tooth + offset) * unit
                     for tooth in range(speed) for offset in (1, 13))
    ordered = sorted(walls)
    twice_period = 2 * denominator
    mass = 0
    components = 0
    was_safe = False
    for left, right in zip(ordered, ordered[1:]):
        twice_midpoint = left + right
        safe = True
        for speed in speeds:
            residue = (speed * twice_midpoint) % twice_period
            if 7 * residue < denominator or 7 * residue > 13 * denominator:
                safe = False
                break
        if safe:
            mass += right - left
            if not was_safe:
                components += 1
        was_safe = safe
    return denominator, mass, components


def direct_margin(q: int, r: int) -> Fraction:
    surplus = Fraction(45 * BODY_MASS - 4 * POOL_D, 63 * POOL_D)
    loss = Fraction(6 * BODY_COMPONENTS, 49) * (Fraction(1, q) + Fraction(1, r))
    return surplus - loss


def main() -> None:
    boundary_records: list[str] = []
    literal_records: list[str] = []
    least_delta: Fraction | None = None
    least_records: list[tuple[int, int, int, int, int]] = []
    largest_period = 0

    for q in range(1290, 1307):
        candidates = (r for r in range(q + 1, 2000) if direct_margin(q, r) >= 0)
        r0 = next(candidates)
        before = direct_margin(q, r0 - 1)
        at = direct_margin(q, r0)
        assert r0 == 2614 - q
        assert before < 0 <= at
        boundary_records.append(
            ",".join((str(q), str(r0), text_fraction(before), text_fraction(at)))
        )

        for r in range(q + 1, r0):
            denominator, mass, components = midpoint_geometry(
                tuple(sorted(EXTREMAL_BODY + (q, r)))
            )
            delta_ticks = 63 * mass - 4 * denominator
            assert delta_ticks > 0
            delta = Fraction(delta_ticks, denominator)
            literal_records.append(
                ",".join(map(str, (q, r, denominator, mass, components, delta_ticks)))
            )
            largest_period = max(largest_period, denominator)
            record = (q, r, denominator, mass, components)
            if least_delta is None or delta < least_delta:
                least_delta = delta
                least_records = [record]
            elif delta == least_delta:
                least_records.append(record)

    assert len(boundary_records) == 17
    assert len(literal_records) == sum(range(1, 34, 2)) == 289
    assert least_delta == Fraction(18_886_235_531, 2_585_198_330)
    assert least_records == [(1300, 1305, 91_205_797_082_400, 16_367_136_156_560, 1044)]

    boundary_bytes = ("\n".join(boundary_records) + "\n").encode("ascii")
    literal_bytes = ("\n".join(literal_records) + "\n").encode("ascii")

    print("LRC14_HYBRID_Q1290_INDEPENDENT_MIDPOINT")
    print(
        "DIRECT_CENSUS MAX_K 1307 MAX_COUNT 1 SECOND_K 1290 "
        "SECOND_COUNT 1 EXCEPTION B_STAR"
    )
    for row in boundary_records:
        q, r0, before, at = row.split(",")
        print(f"R0 Q {q} FIRST_DIRECT_R {r0} PREVIOUS_GAP {before} AT_GAP {at}")
    print(
        f"R0_LEDGER COUNT {len(boundary_records)} FORMULA 2614-Q SHA256 "
        f"{sha256(boundary_bytes).hexdigest()}"
    )
    print(
        f"LITERAL_LEDGER PAIRS {len(literal_records)} STRICT_SAFE {len(literal_records)} "
        f"NONPOSITIVE 0 SHA256 {sha256(literal_bytes).hexdigest()}"
    )
    q_min, r_min, denominator_min, mass_min, components_min = least_records[0]
    print(
        f"MIN_LITERAL Q {q_min} R {r_min} D {denominator_min} MASS {mass_min} "
        f"COMPONENTS {components_min} DELTA_TICKS {63 * mass_min - 4 * denominator_min} "
        f"DELTA {text_fraction(least_delta)} COUNT {len(least_records)}"
    )
    print(f"MAX_DENOMINATOR {largest_period}")
    print(
        "VERDICT ALL_DISTINCT_Q_R_GE_1290_SAFE_RELATIVE_TO_DIRECT_CENSUS "
        "Q1290_NOT_CLAIMED_LITERAL_MINIMAL LRC14_OPEN"
    )


if __name__ == "__main__":
    main()
