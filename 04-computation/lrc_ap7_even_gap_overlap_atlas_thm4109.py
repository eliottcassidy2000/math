#!/usr/bin/env python3
"""Exact certificate for THM-4109's AP7 odd-pair gap atlas.

The script classifies phase resonances,
proves exact pair-overlap quasi-polynomials for the three whole-J resonant
gaps, certifies their uniform exact-second-moment height floors, and runs the
height-263 hostile screen for every gap that could asymptotically improve
THM-4101's original height 264.
"""

from __future__ import annotations

from fractions import Fraction as F


LEFT = F(4, 35)
RIGHT = F(13, 98)
LENGTH = RIGHT - LEFT
DELTA = F(1, 14)
ARC_RADIUS = F(2, 7)
RESIDUE_PERIOD = 980


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def triangle(phase: F) -> F:
    residue = phase % 1
    return max(F(0), ARC_RADIUS - min(residue, 1 - residue))


def averaged_overlap(gap: int) -> F:
    """Q_g=int_J max(0,2/7-||2g theta||) dtheta, exactly."""
    start = 2 * gap * LEFT
    end = 2 * gap * RIGHT
    walls = {start, end}
    first = start.numerator // start.denominator - 2
    last = end.numerator // end.denominator + 3
    for integer in range(first, last + 1):
        for wall in (F(integer) - ARC_RADIUS, F(integer), F(integer) + ARC_RADIUS):
            if start < wall < end:
                walls.add(wall)
    ordered = sorted(walls)
    integral = sum(
        ((right - left) * (triangle(left) + triangle(right)) / 2
         for left, right in zip(ordered, ordered[1:])),
        F(0),
    )
    return integral / (2 * gap)


def local_odd_teeth(speed: int) -> list[tuple[F, F]]:
    roots = 2 * speed
    radius = DELTA / speed
    low = roots * (LEFT - radius)
    high = roots * (RIGHT + radius)
    first = low.numerator // low.denominator - 1
    last = high.numerator // high.denominator + 2
    return [
        (max(LEFT, F(index, roots) - radius),
         min(RIGHT, F(index, roots) + radius))
        for index in range(first, last + 1)
        if F(index, roots) + radius > LEFT
        and F(index, roots) - radius < RIGHT
    ]


def exact_pair_overlap(start: int, gap: int) -> F:
    left_teeth = local_odd_teeth(start)
    right_teeth = local_odd_teeth(start + gap)
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


def whole_interval_resonance(gap: int) -> int | None:
    """Integer k if 2gJ lies in (k-2/7,k+2/7) and contains k."""
    start = 2 * gap * LEFT
    end = 2 * gap * RIGHT
    first = start.numerator // start.denominator - 1
    last = end.numerator // end.denominator + 2
    for integer in range(first, last + 1):
        if (
            start < integer < end
            and F(integer) - start < ARC_RADIUS
            and end - F(integer) < ARC_RADIUS
        ):
            return integer
    return None


def pair_formula_at(
    start: int,
    gap: int,
    mean: F,
    alpha: F,
    beta: F,
) -> F:
    return mean + (alpha * start + beta) / (start * (start + gap))


def residue_forms(gap: int, mean: F) -> dict[int, tuple[F, F]]:
    """Certify u(u+g)(M_g(u)-Q_g)=alpha_r*u+beta_r mod 980.

    For whole-J single-resonance gaps all floor/ceiling boundaries have
    denominators dividing 980.  On a fixed odd residue they are affine, and
    the arithmetic-progression sum of pair-component lengths has the stated
    rational form.  Four exact samples separately check both second
    differences before the form is used.
    """
    forms: dict[int, tuple[F, F]] = {}
    for residue in range(1, RESIDUE_PERIOD, 2):
        starts = [residue + index * RESIDUE_PERIOD for index in range(4)]
        transformed = [
            start * (start + gap) * (exact_pair_overlap(start, gap) - mean)
            for start in starts
        ]
        require(
            transformed[2] - 2 * transformed[1] + transformed[0] == 0,
            f"first affine check failed g={gap}, residue={residue}",
        )
        require(
            transformed[3] - 2 * transformed[2] + transformed[1] == 0,
            f"second affine check failed g={gap}, residue={residue}",
        )
        alpha = (transformed[1] - transformed[0]) / RESIDUE_PERIOD
        beta = transformed[0] - alpha * residue
        require(
            pair_formula_at(starts[3], gap, mean, alpha, beta)
            == exact_pair_overlap(starts[3], gap),
            "held-out formula check failed",
        )
        forms[residue] = (alpha, beta)
    return forms


def exact_uniform_floor(
    gap: int,
    mean: F,
    forms: dict[int, tuple[F, F]],
) -> tuple[int, tuple[F, int, int, str], tuple[F, int, int, str], int]:
    """Smallest B certified by the exact selected-pair second moment.

    For fixed B the worst even outlier is the least even e>=B.  If u is the
    least odd o>=B, distinctness makes the third odd the least allowed odd;
    otherwise that third odd is o.  On each residue ray the margin numerator
    is quadratic, and the exact derivative check below makes the first point
    decisive for the whole infinite ray.
    """
    pair_tax = F(5, 49)

    def pair_value(start: int) -> F:
        alpha, beta = forms[start % RESIDUE_PERIOD]
        return pair_formula_at(start, gap, mean, alpha, beta)

    def certify(height: int) -> tuple[bool, tuple[F, int, int, str], int]:
        least_odd = height if height % 2 else height + 1
        least_even = height if height % 2 == 0 else height + 1
        third = least_odd + 2
        if third == least_odd + gap:
            third += 2
        special_margin = (
            pair_value(least_odd) / 2
            - pair_tax * (F(1, least_odd) + F(1, least_odd + gap) + F(1, third))
            - F(6, 49 * least_even)
        )
        worst = (special_margin, least_odd, third, "special")
        failures = 0
        first_ray_speed = least_odd + 2
        constant_tax = F(5, 49 * least_odd) + F(6, 49 * least_even)
        quadratic_lead = mean / 2 - constant_tax
        for residue, (alpha, beta) in forms.items():
            start = residue
            while start < first_ray_speed:
                start += RESIDUE_PERIOD
            numerator = (
                quadratic_lead * start * (start + gap)
                + (alpha * start + beta) / 2
                - pair_tax * (2 * start + gap)
            )
            derivative = (
                quadratic_lead * (2 * start + gap)
                + alpha / 2
                - 2 * pair_tax
            )
            margin = numerator / (start * (start + gap))
            if margin < worst[0]:
                worst = (margin, start, least_odd, "ray")
            if numerator <= 0 or derivative <= 0:
                failures += 1
        return special_margin > 0 and failures == 0, worst, failures

    floor = 1
    while True:
        success, worst, failures = certify(floor)
        if success:
            previous_success, previous_worst, _ = certify(floor - 1)
            require(not previous_success, "reported floor is not first")
            return floor, worst, previous_worst, failures
        floor += 1


def exact_margin(bank: tuple[int, int, int, int], start: int, gap: int) -> F:
    even = next(speed for speed in bank if speed % 2 == 0)
    third = next(
        speed for speed in bank
        if speed % 2 and speed not in (start, start + gap)
    )
    return (
        exact_pair_overlap(start, gap) / 2
        - F(5, 49) * (F(1, start) + F(1, start + gap) + F(1, third))
        - F(6, 49 * even)
    )


def main() -> None:
    require(LENGTH == F(9, 490), "AP7 length changed")

    zero_gaps = []
    internal = []
    means: dict[int, F] = {}
    for gap in range(2, 202, 2):
        mean = averaged_overlap(gap)
        means[gap] = mean
        if mean == 0:
            zero_gaps.append(gap)
        resonance = whole_interval_resonance(gap)
        if resonance is not None:
            internal.append((gap, resonance))
    require(zero_gaps == [2, 6, 10], "zero-overlap gap classification changed")
    require(internal == [(4, 1), (8, 2), (12, 3)], "whole-J resonance classification changed")

    general_bound_gates = 0
    for gap in range(2, 202, 2):
        lower_error = F(5 * gap, 4) * LENGTH + F(2, 7)
        upper_error = F(5 * gap, 4) * LENGTH + 1
        for start in range(9, 100, 2):
            actual = exact_pair_overlap(start, gap)
            require(
                means[gap] - lower_error / start <= actual
                <= means[gap] + upper_error / start,
                f"moving-target bound failed g={gap},u={start}",
            )
            general_bound_gates += 1

    # The internal-resonance exact quasi-polynomial floors.
    internal_rows = []
    expected_floors = {4: 197, 8: 232, 12: 268}
    expected_means = {
        4: F(2187, 480200),
        8: F(927, 240100),
        12: F(1521, 480200),
    }
    for gap, resonance in internal:
        mean = means[gap]
        require(mean == expected_means[gap], f"overlap mean changed g={gap}")
        forms = residue_forms(gap, mean)
        floor, worst, previous_worst, failures = exact_uniform_floor(gap, mean, forms)
        require(floor == expected_floors[gap], f"uniform floor changed g={gap}")
        require(failures == 0 and worst[0] > 0 and previous_worst[0] <= 0, "floor signs changed")
        internal_rows.append((gap, resonance, mean, floor, worst, previous_worst))

    # Height-263 asymptotic screen.  If Q_g/2 is no larger than the fixed
    # third-odd/even tax, u->infinity is already hostile.  The general
    # moving-target period bound proves M_g(u)->Q_g.
    height = 263
    asymptotic_tax = F(5, 49 * 263) + F(6, 49 * 264)
    candidates = [
        gap for gap in range(2, 202, 2)
        if means[gap] / 2 > asymptotic_tax
    ]
    require(
        all(means[gap] / 2 != asymptotic_tax for gap in range(2, 202, 2)),
        "an equality case needs separate asymptotic treatment",
    )
    require(candidates == [4, 8, 12, 16, 20, 34, 38, 42, 46], "candidate screen changed")

    # For g>=202, one extra whole triangle bounds any residual phase window.
    large_gap_upper = 4 * LENGTH / 49 + F(2, 49 * 202)
    require(large_gap_upper < 2 * asymptotic_tax, "large-gap screen lost strictness")

    hostile_specs = {
        12: (263, 265),
        16: (263, 265),
        20: (263, 265),
        34: (263, 265),
        38: (267, 263),
        42: (267, 263),
        46: (267, 263),
    }
    hostile_rows = []
    for gap, (start, third) in hostile_specs.items():
        bank = tuple(sorted((264, third, start, start + gap)))
        margin = exact_margin(bank, start, gap)
        require(min(bank) >= 263 and margin < 0, f"height-263 hostile failed g={gap}")
        hostile_rows.append((gap, bank, exact_pair_overlap(start, gap), margin))

    print("THM-4109 AP7 even-gap overlap atlas exact certificate: PASS")
    print(f"J=[{LEFT},{RIGHT}] L={LENGTH}")
    print(f"zero_pair_overlap_gaps={zero_gaps}")
    print(f"whole_J_internal_resonances={internal}")
    print(f"general_bound_exact_gates={general_bound_gates}")
    print("internal exact-pair uniform floors:")
    for gap, resonance, mean, floor, worst, previous in internal_rows:
        print(
            f"  g={gap} k={resonance} Q={mean} floor={floor} "
            f"worst={worst} previous_floor_hostile={previous}"
        )
    print(f"height263_asymptotic_tax={asymptotic_tax}")
    print(f"height263_only_possible_gaps={candidates}")
    print(f"g>=202_Q_upper={large_gap_upper}<2tax={2*asymptotic_tax}")
    print("height263 exact criterion hostiles outside g=4,8:")
    for row in hostile_rows:
        print(f"  g={row[0]} bank={row[1]} pair={row[2]} margin={row[3]}")
    print("general moving-target lower bound:")
    print("  M_g(u)>=Q_g-(5*g*L/4+2/7)/u")
    print("general convergence upper bound:")
    print("  M_g(u)<=Q_g+(5*g*L/4+1)/u")
    print("conclusion=gap4 improves 264 to 197; gap8 broadens at 232; no other even gap beats 264")


if __name__ == "__main__":
    main()
