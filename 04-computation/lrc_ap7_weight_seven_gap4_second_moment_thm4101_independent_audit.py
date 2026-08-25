#!/usr/bin/env python3
"""Independent original-theta wall audit of THM-4101.

Unlike the primary probe, this script does not merge comb intervals or use a
periodic-window helper.  It constructs every literal theta wall, classifies
walls and open cells directly, and integrates predicates cell by cell.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction
DELTA = Q(1, 14)
J0 = Q(4, 35)
J1 = Q(13, 98)
CUT = Q(1, 8)
JLEN = J1 - J0


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise RuntimeError(label)


def norm(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def safe_at(speed: int, theta: Q) -> bool:
    return (
        norm(speed * theta) >= DELTA
        and norm(speed * (theta + Q(1, 2))) >= DELTA
    )


def danger_at(speed: int, theta: Q) -> bool:
    return not safe_at(speed, theta)


def count_roots(speed: int) -> int:
    return speed if speed % 2 == 0 else 2 * speed


def wall_rows(speed: int) -> list[tuple[Q, int, str]]:
    roots = count_roots(speed)
    radius = DELTA / speed
    low = roots * (J0 - radius)
    high = roots * (J1 + radius)
    first = low.numerator // low.denominator - 2
    last = high.numerator // high.denominator + 3
    rows: list[tuple[Q, int, str]] = []
    for index in range(first, last + 1):
        for value, side in (
            (Q(index, roots) - radius, "left"),
            (Q(index, roots) + radius, "right"),
        ):
            if J0 <= value <= J1:
                rows.append((value, 14 * speed, f"{speed}-{side}"))
    return rows


def predicate_measure(walls: set[Q], predicate) -> Q:
    ordered = sorted(walls)
    total = Q(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if predicate(midpoint):
            total += right - left
    return total


def pair_measure(first: int, second: int) -> Q:
    walls = {J0, J1}
    walls.update(value for value, _, _ in wall_rows(first))
    walls.update(value for value, _, _ in wall_rows(second))
    return predicate_measure(
        walls,
        lambda theta: danger_at(first, theta) and danger_at(second, theta),
    )


def arrangement(speeds: tuple[int, ...]) -> tuple[Q, Q, int, int]:
    walls = {J0, J1}
    for speed in speeds:
        walls.update(value for value, _, _ in wall_rows(speed))
    ordered = sorted(walls)
    union = Q(0)
    pair_sum = Q(0)
    safe_cells = 0
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        multiplicity = sum(danger_at(speed, midpoint) for speed in speeds)
        if multiplicity:
            union += right - left
        else:
            safe_cells += 1
        pair_sum += Q(multiplicity * (multiplicity - 1), 2) * (right - left)
    uncovered_walls = sum(
        not any(danger_at(speed, theta) for speed in speeds)
        for theta in ordered
    )
    return union, pair_sum, safe_cells, uncovered_walls


def in_open_arc(residue: Q, start: Q, end: Q) -> bool:
    residue %= 1
    if start < 0:
        return residue > 1 + start or residue < end
    return start < residue < end


def frozen_at(gap_start: int, theta: Q) -> bool:
    phase = (2 * gap_start * theta) % 1
    if J0 <= theta <= CUT:
        return in_open_arc(phase, -Q(2, 35), Q(1, 7))
    if CUT <= theta <= J1:
        return in_open_arc(phase, -Q(1, 7), Q(4, 49))
    return False


def frozen_walls(gap_start: int) -> set[Q]:
    multiplier = 2 * gap_start
    walls = {J0, CUT, J1}
    for left, right, arc_left, arc_right in (
        (J0, CUT, -Q(2, 35), Q(1, 7)),
        (CUT, J1, -Q(1, 7), Q(4, 49)),
    ):
        low = multiplier * left - arc_right
        high = multiplier * right - arc_left
        first = low.numerator // low.denominator - 2
        last = high.numerator // high.denominator + 3
        for index in range(first, last + 1):
            for endpoint in (arc_left, arc_right):
                theta = (index + endpoint) / multiplier
                if left <= theta <= right:
                    walls.add(theta)
    return walls


def main() -> None:
    rho_left = Q(1, 5)
    rho_right = Q(11, 49)
    pair_main = rho_left * (CUT - J0) + rho_right * (J1 - CUT)
    pair_error = (
        rho_left * (1 - rho_left) + rho_right * (1 - rho_right)
    ) / 2
    gate("J length", JLEN == Q(9, 490))
    gate("pair main", pair_main == Q(927, 240100))
    gate("pair error", pair_error == Q(10027, 60025))

    # Literal original-theta core check on every core wall and open cell.
    core_walls = {J0, J1}
    for speed in range(1, 8):
        core_walls.update(value for value, _, _ in wall_rows(speed))
    core_points = sorted(core_walls)
    core_tests = core_points + [
        (left + right) / 2 for left, right in zip(core_points, core_points[1:])
    ]
    gate(
        "AP7 original-theta interval",
        all(safe_at(speed, theta) for theta in core_tests for speed in range(1, 8)),
    )

    # Independent pointwise classification.
    pointwise = []
    for multiplicity in range(5):
        left = Q(int(multiplicity > 0))
        right = Q(multiplicity) - Q(multiplicity * (multiplicity - 1), 4)
        gate(f"moment r={multiplicity}", left <= right)
        pointwise.append((multiplicity, left, right))

    # Exact wall/cell containment and overlap replay for 497 odd pairs.
    containment_tests = 0
    overlap_tests = 0
    minimum_slack: tuple[Q, int] | None = None
    for gap_start in range(9, 1002, 2):
        walls = frozen_walls(gap_start)
        walls.update(value for value, _, _ in wall_rows(gap_start))
        walls.update(value for value, _, _ in wall_rows(gap_start + 4))
        ordered = sorted(walls)
        tests = ordered + [
            (left + right) / 2 for left, right in zip(ordered, ordered[1:])
        ]
        for theta in tests:
            if frozen_at(gap_start, theta):
                gate(
                    f"frozen containment u={gap_start}",
                    danger_at(gap_start, theta) and danger_at(gap_start + 4, theta),
                )
            containment_tests += 1
        actual = pair_measure(gap_start, gap_start + 4)
        lower = pair_main - pair_error / gap_start
        gate(f"pair lower u={gap_start}", actual >= lower)
        slack = actual - lower
        if minimum_slack is None or slack < minimum_slack[0]:
            minimum_slack = (slack, gap_start)
        overlap_tests += 1

    # Recompute the two monotonicity boundary cases without importing the
    # primary tax routine.
    half_main = pair_main / 2
    case_265 = (
        Q(6, 49 * 264)
        + Q(5, 49) * (Q(1, 265) + Q(1, 267) + Q(1, 269))
        + pair_error / (2 * 265)
    )
    case_267 = (
        Q(6, 49 * 264)
        + Q(5, 49) * (Q(1, 265) + Q(1, 267) + Q(1, 271))
        + pair_error / (2 * 267)
    )
    margin_265 = half_main - case_265
    margin_267 = half_main - case_267
    gate("uniform case u=265", margin_265 > 0)
    gate("uniform case u>=267", margin_267 > 0)

    hostile = (8, 9, 11, 13)
    hostile_union, hostile_pairs, hostile_safe_cells, hostile_uncovered_walls = arrangement(hostile)
    gate("hostile cells", hostile_safe_cells == 0)
    gate("hostile walls", hostile_uncovered_walls == 0)
    gate("hostile measure", hostile_union == JLEN)

    positive = (264, 265, 267, 269)
    positive_union, positive_pairs, positive_safe_cells, positive_uncovered_walls = arrangement(positive)
    gate("positive survivor cells", positive_safe_cells > 0)
    gate("positive survivor measure", positive_union < JLEN)

    # An arrangement endpoint supplies the adaptive even clock independently.
    candidates = [(J0, 70, "core-left"), (J1, 98, "core-right")]
    for speed in positive:
        candidates.extend(wall_rows(speed))
    body = tuple(range(1, 8)) + positive
    witness = next(
        (theta, clock, source)
        for theta, clock, source in sorted(candidates)
        if all(safe_at(speed, theta) for speed in body)
    )
    gate("witness clock parity", witness[1] % 2 == 0)
    gate("witness label integrality", (witness[0] * witness[1]).denominator == 1)
    gate("witness clock bound", witness[1] <= 14 * max(positive))

    print("THM-4101 independent literal-theta audit: PASS")
    print(f"core_wall_cell_tests={len(core_tests)*7}")
    print(f"pointwise={pointwise}")
    print(
        f"gap4_overlap_rows={overlap_tests} containment_tests={containment_tests} "
        f"minimum_slack={minimum_slack}"
    )
    print(f"constants=P:{pair_main} E:{pair_error} L:{JLEN}")
    print(f"uniform_margins=u265:{margin_265} u267:{margin_267}")
    print(
        f"hostile={hostile} union={hostile_union} pair_sum={hostile_pairs} "
        f"safe_cells={hostile_safe_cells} uncovered_walls={hostile_uncovered_walls}"
    )
    print(
        f"positive={positive} union={positive_union} survivor={JLEN-positive_union} "
        f"pair_sum={positive_pairs} safe_cells={positive_safe_cells} "
        f"uncovered_walls={positive_uncovered_walls}"
    )
    print(f"arrangement_witness={witness} max_clock={14*max(positive)}")


if __name__ == "__main__":
    main()
