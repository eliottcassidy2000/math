#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4151's anchored first-window transfer."""

from __future__ import annotations

import hashlib
from fractions import Fraction
from math import comb, gcd


EXPECTED_SEMANTIC_SHA256 = "84adc34f0dc8ffb415455827a227f86cab1bea3169d73d3912721bd815e1b514"
EXPECTED_SEMANTIC_FNV64 = "ff11c60c57deec04"


def gap(speed: int, phase: Fraction) -> Fraction:
    value = speed * phase
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Fraction) -> Fraction:
    return min(gap(speed, phase) for speed in speeds)


def width(first: int, last: int) -> Fraction:
    return Fraction(13, 14 * last) - Fraction(1, 14 * first)


def anchored_gate(first: int) -> Fraction:
    return Fraction(4 * first - 1, 14 * first * (12 * first + 1))


def last_cap(first: int) -> int:
    return (156 * first + 13) // 16


def eleven_body_count(bound: int) -> int:
    answer = 0
    for first in range(1, bound + 1):
        available = min(bound, last_cap(first)) - first
        if available >= 10:
            answer += comb(available, 10)
    return answer


def translation_minimum(span: int) -> int:
    return max(1, (16 * span - 13 + 139) // 140)


def window_has_safe_lift(first: int, last: int, a: int, b: int) -> bool:
    """Decide a two-tail first window by every exact strict-wall cell."""
    left = Fraction(1, 14 * first)
    right = Fraction(13, 14 * last)
    walls = {left, right}
    for speed in (a, b):
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
        min(gap(a, lift), gap(b, lift)) >= Fraction(1, 14)
        for y in probes
        for lift in (y / 2, (1 + y) / 2)
    )


def physical_walls(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for integer in range(speed + 1):
            for sign in (-1, 1):
                wall = Fraction(14 * integer + sign, 14 * speed)
                if 0 <= wall <= 1:
                    walls.add(wall)
    return tuple(sorted(walls))


def safe_components(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return all exact closed safe components, including isolated walls."""
    delta = Fraction(1, 14)
    ordered = physical_walls(speeds)
    pieces: list[tuple[Fraction, Fraction]] = []
    for left, right in zip(ordered, ordered[1:]):
        middle = (left + right) / 2
        if clearance(speeds, middle) < delta:
            continue
        assert clearance(speeds, left) >= delta
        assert clearance(speeds, right) >= delta
        if pieces and pieces[-1][1] == left:
            pieces[-1] = (pieces[-1][0], right)
        else:
            pieces.append((left, right))
    for wall in ordered:
        if clearance(speeds, wall) < delta:
            continue
        if not any(left <= wall <= right for left, right in pieces):
            pieces.append((wall, wall))
    return tuple(sorted(pieces))


def fnv64(data: bytes) -> str:
    value = 14695981039346656037
    for byte in data:
        value ^= byte
        value = value * 1099511628211 & ((1 << 64) - 1)
    return f"{value:016x}"


def main() -> None:
    delta = Fraction(1, 14)

    # The width inequality is equivalent to the affine cap 16M<=156m+13.
    for first in range(1, 1001):
        cap = last_cap(first)
        assert 16 * cap <= 156 * first + 13 < 16 * (cap + 1)
        assert width(first, cap) >= anchored_gate(first)
        assert width(first, cap + 1) < anchored_gate(first)
        assert cap + 1 < 13 * first

    # At the lower endpoint, the upper lift strictly clears every odd tail
    # through 12m-1.  The first possible hostile is b=12m+1.
    endpoint_floor = Fraction(1, 2)
    for first in range(1, 501):
        y0 = Fraction(1, 14 * first)
        upper_lift = (1 + y0) / 2
        for speed in range(1, 12 * first + 1, 2):
            value = gap(speed, upper_lift)
            assert value > delta
            endpoint_floor = min(endpoint_floor, value)

        hostile = 12 * first + 1
        carrier_left = Fraction(6, 7 * hostile)
        carrier_right = Fraction(8, 7 * hostile)
        lattice_slack = y0 - carrier_left
        remaining_room = carrier_right - y0
        assert lattice_slack == Fraction(1, 14 * first * hostile)
        assert remaining_room == anchored_gate(first)

        # The first maximum outside the cap leaves its whole first window in
        # this open hostile carrier for the pair (1,12m+1).
        first_bad_last = last_cap(first) + 1
        y1 = Fraction(13, 14 * first_bad_last)
        assert carrier_left < y0 <= y1 < carrier_right < Fraction(1, 7)

    # The smallest exact hostile is rescued outside the first window, proving
    # that cap sharpness is mechanism sharpness rather than an unsafe row.
    hostile_left = Fraction(1, 14)
    hostile_right = Fraction(13, 154)
    assert gap(1, hostile_right / 2) == Fraction(13, 308) < delta
    assert gap(13, (1 + hostile_right) / 2) == Fraction(15, 308) < delta
    rescue_phase = Fraction(5, 24)
    rescue_row = tuple(2 * h for h in range(1, 12)) + (1, 13)
    rescue_clearance = clearance(rescue_row, rescue_phase)
    assert rescue_clearance == Fraction(1, 12)

    # Every finite additive pattern is admitted after an explicit translate.
    for span in range(0, 2001):
        first = translation_minimum(span)
        assert first + span <= last_cap(first)
        if first > 1:
            assert first - 1 + span > last_cap(first - 1)
    assert translation_minimum(46) == 6
    assert last_cap(6) == 59
    assert 52 <= last_cap(6)  # {6,...,52} is a 47-label certified block.
    max_block_at_six = last_cap(6) - 6 + 1
    assert max_block_at_six == 54

    # Direct, topology-free replay on every strict-wall cell for 1,355 small
    # primitive rows, plus the first carrier-sharp failure at each minimum.
    primitive_pairs = tuple(
        (p, q)
        for q in range(3, 52, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1
    )
    assert len(primitive_pairs) == 271
    direct_rows = 0
    for first in range(1, 6):
        for p, q in primitive_pairs:
            assert window_has_safe_lift(first, last_cap(first), p, q)
            direct_rows += 1
    assert direct_rows == 1_355
    for first in range(1, 26):
        assert not window_has_safe_lift(
            first, last_cap(first) + 1, 1, 12 * first + 1
        )

    # At m=2 the cap is globally sharp for the theorem's arbitrary-cardinality
    # scope: adding the first excluded maximum produces a literal unsafe row.
    hostile_body = tuple(2 * h for h in range(2, 22))
    hostile_body_components = safe_components(hostile_body)
    assert hostile_body_components == (
        (Fraction(1, 56), Fraction(13, 588)),
        (Fraction(281, 588), Fraction(27, 56)),
        (Fraction(29, 56), Fraction(307, 588)),
        (Fraction(575, 588), Fraction(55, 56)),
    )
    hostile_row = hostile_body + (1, 25)
    hostile_walls = physical_walls(hostile_row)
    hostile_midpoints = tuple(
        (left + right) / 2 for left, right in zip(hostile_walls, hostile_walls[1:])
    )
    assert len(hostile_walls) == 946
    assert max(clearance(hostile_row, wall) for wall in hostile_walls) == Fraction(12, 175)
    assert max(clearance(hostile_row, point) for point in hostile_midpoints) == Fraction(1, 16)
    assert safe_components(hostile_row) == ()

    # Exact finite controls for the positive-density eleven-body family.
    expected_counts = {
        20: 75_582,
        40: 792_864_735,
        80: 3_548_681_310_136,
        120: 392_890_789_426_111,
        160: 10_500_809_430_042_208,
        200: 131_378_242_150_108_190,
    }
    for bound, expected in expected_counts.items():
        assert eleven_body_count(bound) == expected

    density = Fraction(35, 39) ** 10
    assert density == Fraction(2_758_547_353_515_625, 8_140_406_085_191_601)
    assert last_cap(1) == 10 and last_cap(2) == 20
    # THM-4148 adds only the minimum-one body {1,...,11} to the combined gate.
    assert eleven_body_count(20) + 1 == 75_583

    ledger = (
        "gate=(4m-1)/(14m(12m+1));cap=floor((156m+13)/16);"
        "endpoint=odd<=12m;hostile_b=12m+1;slack=1/(14mb);"
        "mechanism_hostile=(m,M,a,b)=(1,11,1,13);rescue=5/24,1/12;"
        "translate=ceil((16D-13)/140);span46_min=6;cap6=59;block6_52=47;"
        "direct_primitive_rows=1355;direct_first_failures=25;"
        "actual_cap_hostile=2{2..21}+{1,25}:unsafe;"
        "hostile_wall_probes=946;wallmax=12/175;midmax=1/16;"
        "anchored_counts20,40,80,120,160,200="
        "75582,792864735,3548681310136,392890789426111,"
        "10500809430042208,131378242150108190;"
        "density=(35/39)^10=2758547353515625/8140406085191601;"
        "combined_add_one_for_N>=11"
    )
    semantic = hashlib.sha256(ledger.encode()).hexdigest()
    semantic_fnv = fnv64(ledger.encode())
    if EXPECTED_SEMANTIC_SHA256:
        assert semantic == EXPECTED_SEMANTIC_SHA256
    if EXPECTED_SEMANTIC_FNV64:
        assert semantic_fnv == EXPECTED_SEMANTIC_FNV64

    print("THM4151_SCALE_SENSITIVE_ANCHORED_FIRST_WINDOW_TRANSFER_20260825")
    print("status=PASS;scope=exact rational and finite combinatorial controls")
    print("gate=(4m-1)/(14m(12m+1));integral_gate=16M<=156m+13")
    print(f"endpoint_test_floor_through_m500={endpoint_floor}")
    print("sharp_carrier=b=12m+1;left_slack=1/(14mb);remaining_room=gate")
    print("mechanism_hostile=(m=1,M=11,tails=(1,13),window=(1/14,13/154))")
    print(f"hostile_rescue=(phase=5/24,clearance={rescue_clearance})")
    print("translation_span_gate=140m>=16D-13")
    print("span46_first_translate=6;cap_at_6=59;block_6_52_labels=47;max_block_at_6=54")
    print("direct_wall_cells=1355_primitive_passes;25_first_failures_confirmed")
    print("actual_cap_hostile=2{2..21}+{1,25};safe_components=0")
    print("actual_cap_hostile_probes=946_walls;wallmax=12/175;midmax=1/16")
    print("anchored_family_counts=" + ",".join(f"{n}:{v}" for n, v in expected_counts.items()))
    print("combined_family_rule=anchored+1_for_N>=11;combined_count_at_20=75583")
    print(f"asymptotic_density={density}")
    print(f"semantic_sha256={semantic}")
    print(f"semantic_fnv64={semantic_fnv}")


if __name__ == "__main__":
    main()
