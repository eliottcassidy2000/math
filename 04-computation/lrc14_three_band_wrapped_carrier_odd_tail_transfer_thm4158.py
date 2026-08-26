#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4158's three-band carrier transfer."""

from __future__ import annotations

import hashlib
from fractions import Fraction
from math import comb, gcd


EXPECTED_SEMANTIC_SHA256 = "b1d8f0a1cad02303d3560e1deb471794d3b9fdcff0b1e1cffe12379674f6cdb6"
EXPECTED_SEMANTIC_FNV64 = "416565bf17102710"

EXPECTED_FINITE_COUNTS = {
    20: 75_582,
    40: 812_850_987,
    80: 3_595_550_244_611,
    120: 397_529_462_747_261,
    160: 10_616_582_432_233_990,
    200: 132_777_517_674_540_845,
}


DELTA = Fraction(1, 14)


def gap(speed: int, phase: Fraction) -> Fraction:
    value = speed * phase
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def carrier(first: int) -> tuple[Fraction, Fraction]:
    return Fraction(1, 14 * first), Fraction(8, 7 * (12 * first + 1))


def carrier_width(first: int) -> Fraction:
    left, right = carrier(first)
    return right - left


def blocks(first: int) -> tuple[tuple[int, int], ...]:
    """Closed-form common-safe bands; first>=2 has exactly three."""
    answer = []
    for tooth in range(3):
        lower = first * (14 * tooth + 1)
        upper = ((12 * first + 1) * (14 * tooth + 13)) // 16
        assert lower <= upper
        answer.append((lower, upper))
    return tuple(answer)


def pool(first: int) -> tuple[int, ...]:
    return tuple(
        speed
        for lower, upper in blocks(first)
        for speed in range(lower, upper + 1)
    )


def interval_safe_for_speed(first: int, speed: int) -> bool:
    """Direct quotient reconstruction, independent of the band formula."""
    lower_denominator = 14 * first
    upper_denominator = 7 * (12 * first + 1)
    lower_integer, lower_remainder = divmod(speed, lower_denominator)
    upper_integer, upper_remainder = divmod(8 * speed, upper_denominator)
    return (
        lower_integer == upper_integer
        and 14 * lower_remainder >= lower_denominator
        and 14 * lower_remainder <= 13 * lower_denominator
        and 14 * upper_remainder >= upper_denominator
        and 14 * upper_remainder <= 13 * upper_denominator
    )


def direct_pool(first: int) -> tuple[int, ...]:
    # A connected safe phase arc has length 6/7.  Above this exact bound no
    # speed can be safe throughout the carrier.
    last_possible = 12 * first * (12 * first + 1) // (4 * first - 1)
    return tuple(
        speed
        for speed in range(1, last_possible + 1)
        if interval_safe_for_speed(first, speed)
    )


def tails_have_safe_lift(
    first: int, tail_a: int, tail_b: int, right_override: Fraction | None = None
) -> bool:
    """Decide both physical sheets on every exact strict-wall cell."""
    left, right = carrier(first)
    if right_override is not None:
        right = right_override
    assert left <= right
    walls = {left, right}
    for speed in (tail_a, tail_b):
        for integer in range(-2, speed + 3):
            for sign in (-1, 1):
                lower_wall = 2 * Fraction(14 * integer + sign, 14 * speed)
                for wall in (lower_wall, lower_wall - 1):
                    if left <= wall <= right:
                        walls.add(wall)
    ordered = sorted(walls)
    probes = set(ordered)
    probes.update((u + v) / 2 for u, v in zip(ordered, ordered[1:]))
    return any(
        min(gap(tail_a, lift), gap(tail_b, lift)) >= DELTA
        for y in probes
        for lift in (y / 2, (1 + y) / 2)
    )


def containing_bad_tooth(speed: int, point: Fraction) -> tuple[Fraction, Fraction] | None:
    phase = speed * point
    integer = (2 * phase.numerator + phase.denominator) // (2 * phase.denominator)
    if abs(phase - integer) >= Fraction(1, 7):
        return None
    return (
        Fraction(7 * integer - 1, 7 * speed),
        Fraction(7 * integer + 1, 7 * speed),
    )


def pool_count_at_bound(first: int, bound: int) -> int:
    return sum(
        max(0, min(bound, upper) - lower + 1)
        for lower, upper in blocks(first)
        if lower <= bound
    )


def eleven_body_count(bound: int) -> int:
    answer = 0
    for first in range(2, bound + 1):
        available = pool_count_at_bound(first, bound) - 1
        if available >= 10:
            answer += comb(available, 10)
    return answer


def linear_power_integral(
    constant: Fraction, slope: Fraction, left: Fraction, right: Fraction
) -> Fraction:
    assert slope != 0
    return (
        (constant + slope * right) ** 11
        - (constant + slope * left) ** 11
    ) / (11 * slope)


def density_constant() -> Fraction:
    pieces = (
        (Fraction(0), Fraction(4, 123), Fraction(0), Fraction(63, 4)),
        (Fraction(4, 123), Fraction(1, 29), Fraction(1), Fraction(-15)),
        (Fraction(1, 29), Fraction(4, 81), Fraction(0), Fraction(14)),
        (Fraction(4, 81), Fraction(1, 15), Fraction(1), Fraction(-25, 4)),
        (Fraction(1, 15), Fraction(4, 39), Fraction(0), Fraction(35, 4)),
        (Fraction(4, 39), Fraction(1), Fraction(1), Fraction(-1)),
    )
    integral = sum(
        linear_power_integral(constant, slope, left, right)
        for left, right, constant, slope in pieces
    )
    return 11 * integral


def fnv64(data: bytes) -> str:
    value = 14695981039346656037
    for byte in data:
        value ^= byte
        value = value * 1099511628211 & ((1 << 64) - 1)
    return f"{value:016x}"


def main() -> None:
    # The common carrier has the exact anchored width inherited from THM-4151.
    for first in range(2, 1001):
        left, right = carrier(first)
        assert left < right
        assert right - left == Fraction(
            4 * first - 1, 14 * first * (12 * first + 1)
        )

    # Endpoint inequalities yield precisely three bands for every m>=2.
    # A fourth band is possible only at m=1; all later bands are empty.
    assert ((12 * 1 + 1) * 55) // 16 == 44
    assert 43 <= 44
    for first in range(2, 1001):
        assert direct_pool(first) == pool(first)
        for speed in pool(first):
            assert interval_safe_for_speed(first, speed)
        fourth_lower = 43 * first
        fourth_upper = ((12 * first + 1) * 55) // 16
        assert fourth_lower > fourth_upper
        fifth_lower = 57 * first
        fifth_upper = ((12 * first + 1) * 69) // 16
        assert fifth_lower > fifth_upper

    # Odd-wall parity: small tails clear at the anchored endpoint; large tails
    # either miss the bad carrier or have no more than G_m room to its right.
    endpoint_floor = Fraction(1, 2)
    large_tooth_checks = 0
    for first in range(2, 301):
        left, right = carrier(first)
        upper_lift = (1 + left) / 2
        for speed in range(1, 12 * first + 1, 2):
            value = gap(speed, upper_lift)
            assert value > DELTA
            endpoint_floor = min(endpoint_floor, value)

        first_large = 12 * first + 1
        exact_tooth = containing_bad_tooth(first_large, left)
        assert exact_tooth == (
            Fraction(6, 7 * first_large),
            Fraction(8, 7 * first_large),
        )
        tooth_left, tooth_right = exact_tooth
        assert left - tooth_left == Fraction(1, 14 * first * first_large)
        assert tooth_right - left == carrier_width(first)
        assert right == tooth_right

        shortened_right = (left + right) / 2
        assert not tails_have_safe_lift(first, 1, first_large, shortened_right)
        assert tails_have_safe_lift(first, 1, first_large)

        for speed in range(first_large, 60 * first + 2, 2):
            tooth = containing_bad_tooth(speed, left)
            if tooth is None:
                continue
            tooth_left, tooth_right = tooth
            assert left - tooth_left >= Fraction(1, 14 * first * speed)
            assert tooth_right - left <= Fraction(
                4 * first - 1, 14 * first * speed
            )
            assert tooth_right - left <= carrier_width(first)
            large_tooth_checks += 1

    # Direct strict-wall replay over all small tail pairs, without using the
    # parity proof.  The body is automatically safe on both lifts throughout.
    direct_tail_rows = 0
    for first in range(2, 11):
        odd_tails = range(1, 76, 2)
        for tail_b in odd_tails:
            for tail_a in range(1, tail_b, 2):
                assert tails_have_safe_lift(first, tail_a, tail_b)
                direct_tail_rows += 1
    assert direct_tail_rows == 6_327

    # The divisor-complete m=7 family lies outside both named first-window
    # gates and kills the old x=1/12 fixed clock.
    pool_seven = pool(7)
    assert blocks(7) == ((7, 69), (105, 143), (203, 217))
    assert len(pool_seven) == 117
    anchors = (7, 120, 126, 143)
    assert set(anchors).issubset(pool_seven)
    assert gcd(gcd(gcd(*anchors[:2]), anchors[2]), anchors[3]) == 1
    owner = {
        2: 120,
        3: 120,
        4: 120,
        5: 120,
        6: 120,
        7: 126,
        8: 120,
        9: 126,
        10: 120,
        11: 143,
        12: 120,
        13: 143,
        14: 126,
    }
    assert all(owner[modulus] in anchors for modulus in range(2, 15))
    assert all(owner[modulus] % modulus == 0 for modulus in range(2, 15))
    divisor_family_count = comb(len(pool_seven) - len(anchors), 7)
    assert divisor_family_count == 38_620_298_376
    assert 16 * 143 > 156 * 7 + 13
    assert 27 * (13 * 7 - 143) - 4 * 7 * 143 < 0
    assert gap(2 * 120, Fraction(1, 12)) == 0

    # Exact finite counts and the six-piece Riemann density constant.
    finite_counts = {
        bound: eleven_body_count(bound)
        for bound in (20, 40, 80, 120, 160, 200)
    }
    assert finite_counts == EXPECTED_FINITE_COUNTS
    density = density_constant()
    expected_density = Fraction(
        848953086913769850118498851618778832628468542103282298683365532079,
        2481088067163593416217816176836483026480276818419826456353950662656,
    )
    assert density == expected_density
    old_density = Fraction(35, 39) ** 10
    density_gain = density - old_density
    assert density_gain > 0

    ledger = (
        "carrier=[1/(14m),8/(7(12m+1))];width=(4m-1)/(14m(12m+1));"
        "bands=[m,floor(13(12m+1)/16)]|[15m,floor(27(12m+1)/16)]|"
        "[29m,floor(41(12m+1)/16)];m>=2;direct_pool_m2..1000=match;"
        f"endpoint_floor_m2..300={endpoint_floor};large_tooth_checks={large_tooth_checks};"
        "shortened_carrier_tails(1,12m+1)=unsafe;closed_endpoint=safe;"
        f"direct_tail_rows={direct_tail_rows};"
        "m7_pool=[7,69]|[105,143]|[203,217];pool_size=117;"
        "anchors=7,120,126,143;divisor_complete_2..14;gcd=1;"
        f"m7_family={divisor_family_count};fails_4148_4151;clock_1/12_killed;"
        "finite_counts="
        + ",".join(f"{bound}:{count}" for bound, count in finite_counts.items())
        + f";density={density};old_density={old_density};gain={density_gain}"
    )
    semantic = hashlib.sha256(ledger.encode()).hexdigest()
    semantic_fnv = fnv64(ledger.encode())
    if EXPECTED_SEMANTIC_SHA256:
        assert semantic == EXPECTED_SEMANTIC_SHA256
    if EXPECTED_SEMANTIC_FNV64:
        assert semantic_fnv == EXPECTED_SEMANTIC_FNV64

    print("THM4158_THREE_BAND_WRAPPED_CARRIER_ODD_TAIL_TRANSFER_20260825")
    print("status=PASS;scope=exact rational and finite combinatorial controls")
    print("carrier=[1/(14m),8/(7(12m+1))];m>=2")
    print("bands=[m,floor(13(12m+1)/16)]|[15m,floor(27(12m+1)/16)]|[29m,floor(41(12m+1)/16)]")
    print("direct_pool_reconstruction=m2..1000_match;fourth_band_only_m1")
    print(f"endpoint_test_floor_through_m300={endpoint_floor}")
    print(f"large_bad_tooth_checks={large_tooth_checks}")
    print("sharp_control=shortening_right_endpoint_fails_for_tails_(1,12m+1)")
    print(f"direct_strict_wall_tail_rows={direct_tail_rows}")
    print("m7_pool=[7,69]|[105,143]|[203,217];size=117")
    print("m7_anchors=7,120,126,143;primitive_divisor_complete_through_14")
    print(f"m7_family_count={divisor_family_count};outside_4148_4151={divisor_family_count}")
    print("fixed_clock_control=gap(240/12)=0")
    print("finite_family_counts=" + ",".join(f"{n}:{v}" for n, v in finite_counts.items()))
    print(f"asymptotic_density={density}")
    print(f"old_density={old_density};density_gain={density_gain}")
    print(f"semantic_sha256={semantic}")
    print(f"semantic_fnv64={semantic_fnv}")


if __name__ == "__main__":
    main()
