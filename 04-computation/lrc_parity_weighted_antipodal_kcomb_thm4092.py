#!/usr/bin/env python3
"""Fraction-exact primary referee for THM-4092.

The universe is the four displayed AP cores, their stated parity rows, a
finite threshold bank, and explicit low-height hostile covers.  The proof
carrier is the literal open antipodal danger comb, not sampled time.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations


DELTA = F(1, 14)
CORE_ROWS = {
    8: (F(11, 49), F(13, 56), 98),
    7: (F(4, 35), F(13, 98), 98),
    6: (F(4, 35), F(1, 7), 70),
    5: (F(4, 35), F(1, 7), 70),
}
THRESHOLDS = {
    (8, 0): 85,
    (8, 1): 106,
    (8, 2): 150,
    (8, 3): 281,
    (7, 0): 63,
    (7, 1): 90,
    (7, 2): 172,
    (6, 0): 76,
    (6, 1): 146,
    (5, 0): 181,
}
HOSTILES = {
    (8, 0): (10, 22, 26),
    (8, 1): (9, 10, 26),
    (8, 2): (9, 10, 11),
    (8, 3): (9, 11, 13),
    (7, 0): (8, 10, 12, 26),
    (7, 1): (8, 9, 10, 12),
    (7, 2): (8, 9, 10, 11),
    (6, 0): (8, 10, 14, 22, 26),
    (6, 1): (7, 8, 10, 12, 26),
    (5, 0): (6, 8, 10, 14, 22, 26),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parity_weight(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_safe(speed: int, theta: F) -> bool:
    return all(
        circle_norm(speed * (theta + F(half, 2))) >= DELTA
        for half in (0, 1)
    )


def local_teeth(speed: int, left: F, right: F) -> list[tuple[F, F, int]]:
    """All lifted open teeth meeting [left,right]."""
    weight = parity_weight(speed)
    count = weight * speed
    radius = DELTA / speed
    lower = count * (left - radius)
    upper = count * (right + radius)
    first = lower.numerator // lower.denominator - 2
    last = upper.numerator // upper.denominator + 3
    return [
        (F(index, count) - radius, F(index, count) + radius, speed)
        for index in range(first, last + 1)
        if F(index, count) + radius > left
        and F(index, count) - radius < right
    ]


def union_measure(intervals: list[tuple[F, F, int]], left: F, right: F) -> F:
    clipped = sorted((max(a, left), min(b, right)) for a, b, _ in intervals if b > left and a < right)
    merged: list[tuple[F, F]] = []
    for start, end in clipped:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return sum((end - start for start, end in merged), F(0))


def covers_open(intervals: list[tuple[F, F, int]], left: F, right: F) -> bool:
    """Whether a finite union of literal open intervals contains [left,right]."""
    active_right: F | None = None
    for start, end, _ in sorted(intervals):
        if active_right is None:
            if start < left:
                active_right = end
            continue
        if start < active_right:
            active_right = max(active_right, end)
        elif active_right <= right:
            return False
        if active_right > right:
            return True
    return active_right is not None and active_right > right


def sharp_window_bound(speed: int, length: F) -> F:
    weight = parity_weight(speed)
    return F(weight, 7) * length + F(7 - weight, 49 * speed)


def norm_extrema(speed: int, left: F, right: F) -> tuple[F, F]:
    points = {left, right}
    low_index = (speed * left).numerator // (speed * left).denominator - 2
    high_index = (speed * right).numerator // (speed * right).denominator + 3
    for index in range(low_index, high_index + 1):
        for point in (F(index, speed), F(2 * index + 1, 2 * speed)):
            if left <= point <= right:
                points.add(point)
    values = [circle_norm(speed * point) for point in points]
    return min(values), max(values)


def criterion(core_size: int, speeds: tuple[int, ...]) -> tuple[bool, F, F, int]:
    left, right, _ = CORE_ROWS[core_size]
    length = right - left
    total_weight = sum(parity_weight(speed) for speed in speeds)
    boundary = sum(
        (F(7 - parity_weight(speed), speed) for speed in speeds),
        F(0),
    )
    budget = 7 * length * (7 - total_weight)
    return total_weight < 7 and boundary < budget, boundary, budget, total_weight


def sample_speeds(core_size: int, odd_count: int) -> tuple[int, ...]:
    count = 11 - core_size
    floor = THRESHOLDS[(core_size, odd_count)]
    for speeds in combinations(range(floor, floor + 2 * count + 5), count):
        if sum(speed % 2 for speed in speeds) == odd_count:
            return speeds
    raise RuntimeError("parity sample not found")


def even_presentation_denominator(value: F) -> int:
    return value.denominator if value.denominator % 2 == 0 else 2 * value.denominator


def arrangement_witness(core_size: int, speeds: tuple[int, ...]) -> tuple[F, int, str]:
    left, right, q_bound = CORE_ROWS[core_size]
    candidates: list[tuple[F, int, str]] = [
        (left, even_presentation_denominator(left), "core_left"),
        (right, even_presentation_denominator(right), "core_right"),
    ]
    require(max(candidates[0][1], candidates[1][1]) <= q_bound, "core clock bound changed")
    for speed in speeds:
        for start, end, _ in local_teeth(speed, left, right):
            if left <= start <= right:
                candidates.append((start, 14 * speed, f"{speed}-left"))
            if left <= end <= right:
                candidates.append((end, 14 * speed, f"{speed}-right"))
    body = tuple(range(1, core_size + 1)) + speeds
    for theta, clock, source in candidates:
        if all(literal_safe(speed, theta) for speed in body):
            return theta, clock, source
    raise RuntimeError(f"no arrangement witness for {core_size},{speeds}")


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    # Exact core controls and the displayed interval atlas.
    core_gates = 0
    for core_size, (left, right, _) in CORE_ROWS.items():
        require(0 < left < right < 1, "invalid core interval")
        for speed in range(1, core_size + 1):
            minimum, maximum = norm_extrema(speed, left, right)
            require(minimum >= DELTA, f"lower core gate failed m={core_size},v={speed}")
            if speed % 2:
                require(maximum <= F(3, 7), f"upper odd gate failed m={core_size},v={speed}")
            core_gates += 1

    require(CORE_ROWS[8][1] - CORE_ROWS[8][0] == F(3, 392), "AP8 length changed")
    require(CORE_ROWS[7][1] - CORE_ROWS[7][0] == F(9, 490), "AP7 length changed")
    require(CORE_ROWS[6][1] - CORE_ROWS[6][0] == F(1, 35), "AP6 length changed")
    require(CORE_ROWS[5][1] - CORE_ROWS[5][0] == F(1, 35), "AP5 length changed")

    # Hostile interval endpoints on many rational grids test the sharp window law.
    fragmentation_gates = 0
    equality_gates = 0
    for speed in range(1, 81):
        weight = parity_weight(speed)
        period = F(1, weight * speed)
        duty = F(weight, 7)
        equality_length = duty * period
        equality_measure = union_measure(
            local_teeth(speed, -equality_length / 2, equality_length / 2),
            -equality_length / 2,
            equality_length / 2,
        )
        require(
            equality_measure == sharp_window_bound(speed, equality_length),
            f"sharp equality control failed v={speed}",
        )
        equality_gates += 1
        for denominator in range(5, 30):
            for numerator in range(denominator):
                left = F(numerator, denominator)
                length = F((3 * numerator + 2) % denominator + 1, denominator * 3)
                right = left + length
                actual = union_measure(local_teeth(speed, left, right), left, right)
                require(actual <= sharp_window_bound(speed, length), f"window bound failed v={speed}")
                fragmentation_gates += 1

    # Arithmetic of every uniform threshold row.
    threshold_gates = 0
    sample_rows: list[tuple[int, int, tuple[int, ...], str, int]] = []
    owner_gates = 0
    dilation_gates = 0
    source_histogram: dict[str, int] = {}
    for (core_size, odd_count), floor in sorted(THRESHOLDS.items()):
        count = 11 - core_size
        total_weight = count + odd_count
        left, right, _ = CORE_ROWS[core_size]
        length = right - left
        rational_threshold = F(7 * count - total_weight, 1) / (7 * length * (7 - total_weight))
        require(floor - 1 <= rational_threshold < floor, "integer threshold rounding failed")
        speeds = sample_speeds(core_size, odd_count)
        valid, boundary, budget, checked_weight = criterion(core_size, speeds)
        require(valid and checked_weight == total_weight and boundary < budget, "sample criterion failed")
        all_teeth = sum((local_teeth(speed, left, right) for speed in speeds), [])
        actual_union = union_measure(all_teeth, left, right)
        require(actual_union < length, "certified row covered its core interval")
        theta, clock, source = arrangement_witness(core_size, speeds)
        require(clock % 2 == 0 and clock <= 14 * max(speeds), "clock bound failed")
        scaled = theta * clock
        require(scaled.denominator == 1, "clock phase is not integral")
        label = scaled.numerator % clock
        opposite = (label + clock // 2) % clock
        body = tuple(range(1, core_size + 1)) + speeds
        for tail_class in range(1, 2 * clock, 2):
            require(
                not (
                    7 * mod_norm(tail_class * label, clock) < clock
                    and 7 * mod_norm(tail_class * opposite, clock) < clock
                ),
                "odd owner class survived complementary labels",
            )
            owner_gates += 1
        for dilation in range(1, 9):
            lifted = theta / dilation
            for speed in body:
                require(literal_safe(dilation * speed, lifted), "dilation transport failed")
                dilation_gates += 1
        source_histogram[source] = source_histogram.get(source, 0) + 1
        sample_rows.append((core_size, odd_count, speeds, source, clock))
        threshold_gates += 1

    # Low-height controls show that the selected intervals can genuinely be covered.
    hostile_gates = 0
    for (core_size, odd_count), speeds in HOSTILES.items():
        left, right, _ = CORE_ROWS[core_size]
        require(sum(speed % 2 for speed in speeds) == odd_count, "hostile parity changed")
        intervals = sum((local_teeth(speed, left, right) for speed in speeds), [])
        require(covers_open(intervals, left, right), f"hostile stopped covering {core_size},{odd_count}")
        hostile_gates += 1
    ceiling = (8, 9, 11, 13)
    require(sum(parity_weight(speed) for speed in ceiling) == 7, "weighted ceiling changed")
    left7, right7, _ = CORE_ROWS[7]
    require(
        covers_open(sum((local_teeth(speed, left7, right7) for speed in ceiling), []), left7, right7),
        "weight-seven hostile stopped covering",
    )

    print("THM-4092 primary Fraction-exact audit: PASS")
    print(f"core gates={core_gates}; fragmentation gates={fragmentation_gates}; sharp equalities={equality_gates}")
    print(f"threshold rows={threshold_gates}; owner gates={owner_gates}; dilation gates={dilation_gates}")
    print(f"low-height hostile covers={hostile_gates}; weight-seven hostile={ceiling}")
    print("uniform threshold table (core,odd,min_outlier):")
    print("  " + ", ".join(f"AP{m}/o{o}:{b}" for (m, o), b in sorted(THRESHOLDS.items())))
    print("sample certified rows:")
    for row in sample_rows:
        print(f"  AP{row[0]}/o{row[1]} speeds={row[2]} source={row[3]} N={row[4]}")
    print(f"arrangement source histogram={source_histogram}")


if __name__ == "__main__":
    main()
