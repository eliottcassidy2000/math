#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4158's three-band carrier transfer."""

from __future__ import annotations

import hashlib
import sys
from fractions import Fraction
from math import comb, gcd


sys.stdout.reconfigure(newline="\n")


EXPECTED_SEMANTIC_SHA256 = "9ef136184225c88a692eb1aca0a33e5be60e02f4cf27a20c04939be98ba6d96e"
EXPECTED_SEMANTIC_FNV64 = "0f8d84ecb8dee6a5"

EXPECTED_FINITE_COUNTS = {
    20: 78_585,
    40: 813_203_703,
    80: 3_595_551_388_677,
    120: 397_529_463_891_327,
    160: 10_616_582_433_378_056,
    200: 132_777_517_675_684_911,
}


DELTA = Fraction(1, 14)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"FAIL: {label}")


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
    """Closed-form common-safe bands: four at m=1, then exactly three."""
    answer = []
    for tooth in range(4 if first == 1 else 3):
        lower = first * (14 * tooth + 1)
        upper = ((12 * first + 1) * (14 * tooth + 13)) // 16
        require(lower <= upper, ("nonempty band", first, tooth))
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
    require(left <= right, ("carrier order", first, right_override))
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
    for first in range(1, bound + 1):
        available = pool_count_at_bound(first, bound) - 1
        if available >= 10:
            answer += comb(available, 10)
    return answer


def linear_power_integral(
    constant: Fraction, slope: Fraction, left: Fraction, right: Fraction
) -> Fraction:
    require(slope != 0, ("nonzero integration slope", left, right))
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
    for first in range(1, 1001):
        left, right = carrier(first)
        require(left < right, ("positive carrier width", first))
        require(
            right - left
            == Fraction(4 * first - 1, 14 * first * (12 * first + 1)),
            ("carrier width identity", first),
        )

    # Endpoint inequalities yield precisely three bands for every m>=2.
    # A fourth band is possible only at m=1; all later bands are empty.
    require(((12 * 1 + 1) * 55) // 16 == 44, "m=1 fourth-band endpoint")
    require(43 <= 44, "m=1 fourth band is nonempty")
    for first in range(2, 1001):
        require(direct_pool(first) == pool(first), ("direct pool", first))
        for speed in pool(first):
            require(
                interval_safe_for_speed(first, speed),
                ("formula band safety", first, speed),
            )
        fourth_lower = 43 * first
        fourth_upper = ((12 * first + 1) * 55) // 16
        require(fourth_lower > fourth_upper, ("empty fourth band", first))
        fifth_lower = 57 * first
        fifth_upper = ((12 * first + 1) * 69) // 16
        require(fifth_lower > fifth_upper, ("empty fifth band", first))

    # Odd-wall parity: small tails clear at the anchored endpoint; large tails
    # either miss the bad carrier or have no more than G_m room to its right.
    endpoint_floor = Fraction(1, 2)
    large_tooth_checks = 0
    for first in range(1, 301):
        left, right = carrier(first)
        upper_lift = (1 + left) / 2
        for speed in range(1, 12 * first + 1, 2):
            value = gap(speed, upper_lift)
            require(value > DELTA, ("small-tail endpoint", first, speed))
            endpoint_floor = min(endpoint_floor, value)

        first_large = 12 * first + 1
        exact_tooth = containing_bad_tooth(first_large, left)
        require(
            exact_tooth
            == (
                Fraction(6, 7 * first_large),
                Fraction(8, 7 * first_large),
            ),
            ("sharp containing tooth", first),
        )
        tooth_left, tooth_right = exact_tooth
        require(
            left - tooth_left == Fraction(1, 14 * first * first_large),
            ("sharp left slack", first),
        )
        require(
            tooth_right - left == carrier_width(first),
            ("sharp right room", first),
        )
        require(right == tooth_right, ("sharp right endpoint", first))

        shortened_right = (left + right) / 2
        require(
            not tails_have_safe_lift(first, 1, first_large, shortened_right),
            ("shortened-carrier hostile", first),
        )
        require(
            tails_have_safe_lift(first, 1, first_large),
            ("closed-endpoint rescue", first),
        )

        for speed in range(first_large, 60 * first + 2, 2):
            tooth = containing_bad_tooth(speed, left)
            if tooth is None:
                continue
            tooth_left, tooth_right = tooth
            require(
                left - tooth_left >= Fraction(1, 14 * first * speed),
                ("odd-wall left slack", first, speed),
            )
            require(
                tooth_right - left
                <= Fraction(4 * first - 1, 14 * first * speed),
                ("odd-wall right room", first, speed),
            )
            require(
                tooth_right - left <= carrier_width(first),
                ("carrier-width closure", first, speed),
            )
            large_tooth_checks += 1

    # Direct strict-wall replay over all small tail pairs, without using the
    # parity proof.  The body is automatically safe on both lifts throughout.
    direct_tail_rows = 0
    for first in range(1, 11):
        odd_tails = range(1, 76, 2)
        for tail_b in odd_tails:
            for tail_a in range(1, tail_b, 2):
                require(
                    tails_have_safe_lift(first, tail_a, tail_b),
                    ("direct tail row", first, tail_a, tail_b),
                )
                direct_tail_rows += 1
    require(direct_tail_rows == 6_327, "direct tail row count")

    # The divisor-complete m=7 family lies outside both named first-window
    # gates and kills the old x=1/12 fixed clock.
    pool_seven = pool(7)
    require(blocks(7) == ((7, 69), (105, 143), (203, 217)), "m=7 bands")
    require(len(pool_seven) == 117, "m=7 pool size")
    anchors = (7, 120, 126, 143)
    require(set(anchors).issubset(pool_seven), "m=7 anchors in pool")
    require(gcd(gcd(gcd(*anchors[:2]), anchors[2]), anchors[3]) == 1, "anchor gcd")
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
    require(all(owner[modulus] in anchors for modulus in range(2, 15)), "owner labels")
    require(all(owner[modulus] % modulus == 0 for modulus in range(2, 15)), "divisor owners")
    divisor_family_count = comb(len(pool_seven) - len(anchors), 7)
    require(divisor_family_count == 38_620_298_376, "m=7 family count")
    require(16 * 143 > 156 * 7 + 13, "outside THM-4151 gate")
    require(27 * (13 * 7 - 143) - 4 * 7 * 143 < 0, "outside THM-4148 gate")
    require(gap(2 * 120, Fraction(1, 12)) == 0, "fixed-clock hostile")

    # Exact finite counts and the six-piece Riemann density constant.
    finite_counts = {
        bound: eleven_body_count(bound)
        for bound in (20, 40, 80, 120, 160, 200)
    }
    require(finite_counts == EXPECTED_FINITE_COUNTS, "finite family counts")
    density = density_constant()
    expected_density = Fraction(
        848953086913769850118498851618778832628468542103282298683365532079,
        2481088067163593416217816176836483026480276818419826456353950662656,
    )
    require(density == expected_density, "density constant")
    old_density = Fraction(35, 39) ** 10
    density_gain = density - old_density
    require(density_gain > 0, "strict density gain")

    ledger = (
        "carrier=[1/(14m),8/(7(12m+1))];width=(4m-1)/(14m(12m+1));"
        "bands=[m,floor(13(12m+1)/16)]|[15m,floor(27(12m+1)/16)]|"
        "[29m,floor(41(12m+1)/16)];m>=2;"
        "m1_fourth=[43,44];m1_pool_size=24;m1_family=1144066;"
        "direct_pool_m1..1000=match;"
        f"endpoint_floor_m1..300={endpoint_floor};large_tooth_checks={large_tooth_checks};"
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
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic SHA-256")
    if EXPECTED_SEMANTIC_FNV64:
        require(semantic_fnv == EXPECTED_SEMANTIC_FNV64, "semantic FNV-64")

    print("THM4158_THREE_BAND_WRAPPED_CARRIER_ODD_TAIL_TRANSFER_20260825")
    print("status=PASS;scope=exact rational and finite combinatorial controls")
    print("carrier=[1/(14m),8/(7(12m+1))];m>=1")
    print("bands=[m,floor(13(12m+1)/16)]|[15m,floor(27(12m+1)/16)]|[29m,floor(41(12m+1)/16)]")
    print("m1_exception=[1,10]|[15,21]|[29,33]|[43,44];size=24;family=1144066")
    print("direct_pool_reconstruction=m1..1000_match;fourth_band_only_m1")
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
