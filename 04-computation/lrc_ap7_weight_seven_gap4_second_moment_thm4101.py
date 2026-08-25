#!/usr/bin/env python3
"""Fraction-exact primary referee for THM-4101.

This is not a claim of LRC(14).  Its analytic carrier is
the AP7 interval J=[4/35,13/98], the four-comb pointwise second-moment
inequality, and a forced overlap between odd danger combs U_u and U_(u+4).
The finite universe below is hostile/positive controls plus all odd
u<=2001 for the pair-overlap replay.  The uniform >=264 conclusion itself
is symbolic: only two exact monotonicity boundary cases are needed.
"""

from __future__ import annotations

from fractions import Fraction as F


DELTA = F(1, 14)
LEFT = F(4, 35)
TURN = F(1, 8)
RIGHT = F(13, 98)
LENGTH = RIGHT - LEFT

# Frozen phase arcs on the two sides of 1/8.
RHO_MINUS = F(1, 5)
RHO_PLUS = F(11, 49)
PAIR_MAIN = RHO_MINUS * (TURN - LEFT) + RHO_PLUS * (RIGHT - TURN)
PAIR_ERROR = (
    RHO_MINUS * (1 - RHO_MINUS) + RHO_PLUS * (1 - RHO_PLUS)
) / 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def parity_weight(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def literal_safe(speed: int, theta: F) -> bool:
    return all(
        circle_norm(speed * (theta + F(half, 2))) >= DELTA
        for half in (0, 1)
    )


def local_teeth(speed: int) -> list[tuple[F, F]]:
    count = parity_weight(speed) * speed
    radius = DELTA / speed
    low = count * (LEFT - radius)
    high = count * (RIGHT + radius)
    first = low.numerator // low.denominator - 2
    last = high.numerator // high.denominator + 3
    return sorted(
        (max(LEFT, F(index, count) - radius),
         min(RIGHT, F(index, count) + radius))
        for index in range(first, last + 1)
        if F(index, count) + radius > LEFT
        and F(index, count) - radius < RIGHT
    )


def merge_open(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    merged: list[tuple[F, F]] = []
    for start, end in sorted(intervals):
        # Literal teeth are open.  A touching endpoint remains a safe seam.
        if merged and start < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def union_measure(speeds: tuple[int, ...]) -> F:
    return sum(
        (end - start for start, end in merge_open([
            tooth for speed in speeds for tooth in local_teeth(speed)
        ])),
        F(0),
    )


def intersection_measure(first: int, second: int) -> F:
    left_teeth = local_teeth(first)
    right_teeth = local_teeth(second)
    left_index = 0
    right_index = 0
    total = F(0)
    while left_index < len(left_teeth) and right_index < len(right_teeth):
        left_piece = left_teeth[left_index]
        right_piece = right_teeth[right_index]
        total += max(
            F(0),
            min(left_piece[1], right_piece[1])
            - max(left_piece[0], right_piece[0]),
        )
        if left_piece[1] < right_piece[1]:
            left_index += 1
        else:
            right_index += 1
    return total


def periodic_arc_measure(
    multiplier: int,
    arc_left: F,
    arc_right: F,
    interval_left: F,
    interval_right: F,
) -> F:
    """Measure where multiplier*theta mod 1 lies in one lifted open arc."""
    low = multiplier * interval_left - arc_right
    high = multiplier * interval_right - arc_left
    first = low.numerator // low.denominator - 2
    last = high.numerator // high.denominator + 3
    pieces = []
    for index in range(first, last + 1):
        start = F(index, multiplier) + arc_left / multiplier
        end = F(index, multiplier) + arc_right / multiplier
        if end > interval_left and start < interval_right:
            pieces.append((max(start, interval_left), min(end, interval_right)))
    return sum((end - start for start, end in pieces), F(0))


def first_moment_tax(speeds: tuple[int, ...]) -> F:
    return sum(
        (F(7 - parity_weight(speed), 49 * speed) for speed in speeds),
        F(0),
    )


def analytic_union_upper(speeds: tuple[int, ...], gap_start: int) -> F:
    return (
        LENGTH
        + first_moment_tax(speeds)
        - (PAIR_MAIN - PAIR_ERROR / gap_start) / 2
    )


def open_covers(speeds: tuple[int, ...]) -> bool:
    walls = {LEFT, RIGHT}
    for speed in speeds:
        for start, end in local_teeth(speed):
            walls.add(start)
            walls.add(end)
    ordered = sorted(walls)
    tests = ordered + [
        (start + end) / 2 for start, end in zip(ordered, ordered[1:])
    ]
    return all(
        any(not literal_safe(speed, theta) for speed in speeds)
        for theta in tests
    )


def norm_extrema(speed: int) -> tuple[F, F]:
    points = {LEFT, RIGHT}
    low = (speed * LEFT).numerator // (speed * LEFT).denominator - 2
    high = (speed * RIGHT).numerator // (speed * RIGHT).denominator + 3
    for index in range(low, high + 1):
        for point in (F(index, speed), F(2 * index + 1, 2 * speed)):
            if LEFT <= point <= RIGHT:
                points.add(point)
    values = [circle_norm(speed * point) for point in points]
    return min(values), max(values)


def main() -> None:
    require(LENGTH == F(9, 490), "AP7 interval length changed")
    require(TURN - LEFT == F(3, 280), "left split length changed")
    require(RIGHT - TURN == F(3, 392), "right split length changed")
    require(PAIR_MAIN == F(927, 240100), "pair main constant changed")
    require(PAIR_ERROR == F(10027, 60025), "pair error constant changed")

    # The chosen interval is genuinely AP7 antipodal-safe, not sampled-safe.
    for speed in range(1, 8):
        minimum, maximum = norm_extrema(speed)
        require(minimum >= DELTA, f"AP7 lower gate failed at {speed}")
        if speed % 2:
            require(maximum <= F(3, 7), f"AP7 odd upper gate failed at {speed}")

    # Pointwise, for r in {0,...,4}, 1_(r>0)<=r-C(r,2)/2.
    moment_values = []
    for multiplicity in range(5):
        indicator = int(multiplicity > 0)
        upper = F(multiplicity) - F(multiplicity * (multiplicity - 1), 4)
        require(F(indicator) <= upper, f"pointwise moment gate failed at {multiplicity}")
        moment_values.append((multiplicity, indicator, upper))

    # Lower-discrepancy replay.  For a period-1/n arc of duty rho on a
    # window of length ell, the exact minimum is at least
    # rho*ell-rho*(1-rho)/n.  The two frozen arcs are contained in
    # U_u cap U_(u+4) because 8*theta-1 ranges through
    # [-3/35,0] and [0,3/49] on the two subwindows.
    overlap_rows = 0
    worst_pair_slack: tuple[F, int] | None = None
    for gap_start in range(9, 2002, 2):
        multiplier = 2 * gap_start
        frozen = periodic_arc_measure(
            multiplier, -F(2, 35), F(1, 7), LEFT, TURN
        ) + periodic_arc_measure(
            multiplier, -F(1, 7), F(4, 49), TURN, RIGHT
        )
        lower = PAIR_MAIN - PAIR_ERROR / gap_start
        actual = intersection_measure(gap_start, gap_start + 4)
        require(frozen >= lower, f"periodic lower discrepancy failed at u={gap_start}")
        require(actual >= frozen, f"frozen pair containment failed at u={gap_start}")
        slack = actual - lower
        if worst_pair_slack is None or slack < worst_pair_slack[0]:
            worst_pair_slack = (slack, gap_start)
        overlap_rows += 1

    # Uniform minimum 264: e>=264, while odd speeds are >=265.
    # If u=265, distinctness forces the third odd speed w>=267.
    case_u265 = (
        F(6, 49 * 264)
        + F(5, 49) * (F(1, 265) + F(1, 267) + F(1, 269))
        + PAIR_ERROR / (2 * 265)
    )
    # If u>=267, monotonicity puts u at 267 and permits w=265.
    case_u267 = (
        F(6, 49 * 264)
        + F(5, 49) * (F(1, 265) + F(1, 267) + F(1, 271))
        + PAIR_ERROR / (2 * 267)
    )
    half_main = PAIR_MAIN / 2
    margin_u265 = half_main - case_u265
    margin_u267 = half_main - case_u267
    require(margin_u265 == F(489784241, 100536614409000), "u=265 margin changed")
    require(margin_u267 == F(203219143, 20256819706200), "u=267 margin changed")
    require(margin_u265 > 0 and margin_u267 > 0, "uniform >=264 gate failed")

    # The adjacent lower floor 263 is not certified by this frozen-arc bound.
    near_miss = (264, 263, 265, 267)
    near_miss_upper = analytic_union_upper(near_miss, 263)
    require(near_miss_upper >= LENGTH, "the recorded 263 near-miss disappeared")
    require(union_measure(near_miss) < LENGTH, "263 control unexpectedly covered J")

    hostile = (8, 9, 11, 13)
    require(open_covers(hostile), "THM-4092 weight-seven hostile stopped covering J")
    require(union_measure(hostile) == LENGTH, "hostile cover measure changed")
    require(analytic_union_upper(hostile, 9) > LENGTH, "hostile was falsely certified")

    boundary_positive = (264, 265, 267, 269)
    require(analytic_union_upper(boundary_positive, 265) < LENGTH, "boundary bank not certified")
    require(union_measure(boundary_positive) < LENGTH, "boundary bank covered J")

    separated_positive = (264, 265, 269, 1001)
    require(analytic_union_upper(separated_positive, 265) < LENGTH, "separated bank not certified")
    require(union_measure(separated_positive) < LENGTH, "separated bank covered J")

    print("THM-4101 AP7 weight-seven gap-four primary audit: PASS")
    print(f"J=[{LEFT},{RIGHT}] L={LENGTH} split={TURN}")
    print(f"pair_lower=P-E/u P={PAIR_MAIN} E={PAIR_ERROR}")
    print("pointwise_moment=" + ",".join(
        f"r{row[0]}:{row[1]}<={row[2]}" for row in moment_values
    ))
    print(f"pair_overlap_rows={overlap_rows} worst_exact_slack={worst_pair_slack}")
    print(f"uniform_floor=264 margin_u265={margin_u265} margin_u267={margin_u267}")
    print(
        f"hostile={hostile} union={union_measure(hostile)} "
        f"analytic_upper={analytic_union_upper(hostile,9)} covers_open=1"
    )
    print(
        f"boundary_positive={boundary_positive} union={union_measure(boundary_positive)} "
        f"survivor_measure={LENGTH-union_measure(boundary_positive)} "
        f"analytic_slack={LENGTH-analytic_union_upper(boundary_positive,265)}"
    )
    print(
        f"separated_positive={separated_positive} union={union_measure(separated_positive)} "
        f"analytic_slack={LENGTH-analytic_union_upper(separated_positive,265)}"
    )
    print(
        f"floor263_near_miss={near_miss} actual_survivor={LENGTH-union_measure(near_miss)} "
        f"analytic_deficit={LENGTH-near_miss_upper} status=criterion-only-not-counterexample"
    )


if __name__ == "__main__":
    main()
