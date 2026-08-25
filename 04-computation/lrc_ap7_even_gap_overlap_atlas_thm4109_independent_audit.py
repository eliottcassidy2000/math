#!/usr/bin/env python3
"""Independent exact referee for THM-4109's AP7 even-gap overlap atlas.

This implementation does not merge the two original combs.  It has two
separate paths:

* a literal wall-cell integrator in the original theta coordinate; and
* for the three single-packet resonances, a direct sum of the two endpoint
  families on either side of theta=1/8.

The periodic triangle integral is evaluated by an exact antiderivative, not
by the wall quadrature used in the primary atlas.
"""

from __future__ import annotations

from fractions import Fraction as F


LEFT = F(4, 35)
MID = F(1, 8)
RIGHT = F(13, 98)
LENGTH = RIGHT - LEFT
PERIOD = 980
ODD_TAX = F(5, 49)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -floor_fraction(-value)


def circle_distance(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def odd_hit(speed: int, theta: F) -> bool:
    return circle_distance(2 * speed * theta) < F(1, 7)


def even_hit(speed: int, theta: F) -> bool:
    return circle_distance(speed * theta) < F(1, 14)


def odd_walls(speed: int) -> set[F]:
    phase_left = 2 * speed * LEFT
    phase_right = 2 * speed * RIGHT
    walls: set[F] = set()
    for index in range(floor_fraction(phase_left) - 2,
                       ceil_fraction(phase_right) + 3):
        for sign in (-1, 1):
            wall = F(7 * index + sign, 14 * speed)
            if LEFT < wall < RIGHT:
                walls.add(wall)
    return walls


def even_walls(speed: int) -> set[F]:
    phase_left = speed * LEFT
    phase_right = speed * RIGHT
    walls: set[F] = set()
    for index in range(floor_fraction(phase_left) - 2,
                       ceil_fraction(phase_right) + 3):
        for sign in (-1, 1):
            wall = F(14 * index + sign, 14 * speed)
            if LEFT < wall < RIGHT:
                walls.add(wall)
    return walls


def wall_pair_measure(start: int, gap: int) -> F:
    walls = {LEFT, RIGHT}
    walls.update(odd_walls(start))
    walls.update(odd_walls(start + gap))
    ordered = sorted(walls)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if odd_hit(start, midpoint) and odd_hit(start + gap, midpoint):
            total += right - left
    return total


def bank_union_measure(bank: tuple[int, int, int, int]) -> F:
    walls = {LEFT, RIGHT}
    for speed in bank:
        walls.update(even_walls(speed) if speed % 2 == 0 else odd_walls(speed))
    ordered = sorted(walls)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if any(
            even_hit(speed, midpoint) if speed % 2 == 0
            else odd_hit(speed, midpoint)
            for speed in bank
        ):
            total += right - left
    return total


def triangle_primitive(value: F) -> F:
    """Integral of max(0,2/7-||s||) from 0 to value."""
    integer = floor_fraction(value)
    residue = value - integer
    whole = F(4 * integer, 49)
    if residue <= F(2, 7):
        partial = F(2, 7) * residue - residue * residue / 2
    elif residue <= F(5, 7):
        partial = F(2, 49)
    else:
        partial = F(2, 49) + (residue - F(5, 7)) ** 2 / 2
    return whole + partial


def averaged_overlap(gap: int) -> F:
    return (
        triangle_primitive(2 * gap * RIGHT)
        - triangle_primitive(2 * gap * LEFT)
    ) / (2 * gap)


def packet_component_measure(start: int, gap: int) -> F:
    """Sum the explicit packet endpoints, independently of comb merging."""
    require(gap in (4, 8, 12), "single-packet formula only")
    n = 2 * start
    m = 2 * (start + gap)
    resonance = gap // 4
    first = floor_fraction(n * LEFT) - 3
    last = ceil_fraction(n * RIGHT) + 3
    total = F(0)
    for index in range(first, last + 1):
        # Before 1/8: the n-tooth indexed by p meets the m-tooth p+k.
        raw_left = F(index, 1) + resonance - F(1, 7)
        raw_right = F(index, 1) + F(1, 7)
        left = max(LEFT, raw_left / m)
        right = min(MID, raw_right / n)
        if left < right:
            total += right - left

        # After 1/8 the other two endpoints are active.
        raw_left = F(index, 1) - F(1, 7)
        raw_right = F(index, 1) + resonance + F(1, 7)
        left = max(MID, raw_left / n)
        right = min(RIGHT, raw_right / m)
        if left < right:
            total += right - left
    return total


def residue_forms(gap: int, mean: F) -> tuple[dict[int, tuple[F, F]], int]:
    """Derive at remote samples and audit both backward and forward."""
    forms: dict[int, tuple[F, F]] = {}
    gates = 0
    for residue in range(1, PERIOD, 2):
        derive_left = residue + 2 * PERIOD
        derive_right = residue + 3 * PERIOD

        def transformed(start: int) -> F:
            return (
                start * (start + gap)
                * (packet_component_measure(start, gap) - mean)
            )

        left_value = transformed(derive_left)
        right_value = transformed(derive_right)
        alpha = (right_value - left_value) / PERIOD
        beta = left_value - alpha * derive_left
        for start in (
            residue,
            residue + PERIOD,
            derive_left,
            derive_right,
            residue + 4 * PERIOD,
            residue + 7 * PERIOD,
        ):
            require(
                transformed(start) == alpha * start + beta,
                f"packet affine law failed g={gap},r={residue},u={start}",
            )
            gates += 1

        # Literal original-theta cells independently check a low and remote
        # representative of every residue class.
        for start in (residue, derive_left):
            require(
                wall_pair_measure(start, gap)
                == packet_component_measure(start, gap),
                f"wall/packet disagreement g={gap},r={residue},u={start}",
            )
            gates += 1
        forms[residue] = (alpha, beta)
    return forms, gates


def formula_pair(
    start: int,
    gap: int,
    mean: F,
    forms: dict[int, tuple[F, F]],
) -> F:
    alpha, beta = forms[start % PERIOD]
    return mean + (alpha * start + beta) / (start * (start + gap))


def certificate_margin(pair: F, start: int, gap: int, third: int, even: int) -> F:
    return (
        pair / 2
        - ODD_TAX * (F(1, start) + F(1, start + gap) + F(1, third))
        - F(6, 49 * even)
    )


def audit_uniform_height(
    gap: int,
    height: int,
    mean: F,
    forms: dict[int, tuple[F, F]],
) -> tuple[F, int, int, int]:
    """Certify every distinct bank at a fixed lower speed height."""
    least_odd = height if height % 2 else height + 1
    least_even = height if height % 2 == 0 else height + 1

    third = least_odd + 2
    if third == least_odd + gap:
        third += 2
    special = certificate_margin(
        formula_pair(least_odd, gap, mean, forms),
        least_odd,
        gap,
        third,
        least_even,
    )
    require(special > 0, f"special placement failed g={gap},B={height}")

    # If the selected pair does not start at the least odd speed, the third
    # odd tax is maximal at that least speed.  Each residue ray then has a
    # quadratic numerator.  Positive leading coefficient, initial value and
    # initial derivative make it positive on the entire continuous ray, and
    # hence on the required 980-step integer subray.
    ray_floor = least_odd + 2
    constant_tax = ODD_TAX / least_odd + F(6, 49 * least_even)
    lead = mean / 2 - constant_tax
    require(lead > 0, f"ray leading coefficient failed g={gap},B={height}")
    least_gate = special
    least_start = least_odd
    for residue, (alpha, beta) in forms.items():
        start = residue
        while start < ray_floor:
            start += PERIOD
        numerator = (
            lead * start * (start + gap)
            + (alpha * start + beta) / 2
            - ODD_TAX * (2 * start + gap)
        )
        derivative = (
            lead * (2 * start + gap)
            + alpha / 2
            - 2 * ODD_TAX
        )
        require(
            numerator > 0 and derivative > 0,
            f"ray gate failed g={gap},B={height},r={residue},u={start}",
        )
        boundary_margin = numerator / (start * (start + gap))
        if boundary_margin < least_gate:
            least_gate = boundary_margin
            least_start = start
    return least_gate, least_start, least_odd, least_even


def first_height_hostile(
    gap: int,
    height: int,
    limit: int = 2001,
) -> tuple[tuple[int, int, int, int], int, F, F]:
    """Find a literal exact failure of the selected-pair margin."""
    least_even = height if height % 2 == 0 else height + 1
    least_odd = height if height % 2 else height + 1
    for start in range(least_odd, limit + 1, 2):
        third = least_odd
        if third in (start, start + gap):
            third += 2
            if third in (start, start + gap):
                third += 2
        pair = wall_pair_measure(start, gap)
        margin = certificate_margin(pair, start, gap, third, least_even)
        if margin <= 0:
            bank = tuple(sorted((least_even, third, start, start + gap)))
            require(len(set(bank)) == 4 and min(bank) >= height, "bad hostile bank")
            return bank, start, pair, margin
    raise RuntimeError(f"no hostile found g={gap},B={height},limit={limit}")


def main() -> None:
    require(LENGTH == F(9, 490), "AP7 interval changed")
    require(
        [triangle_primitive(F(k)) for k in range(5)]
        == [F(4 * k, 49) for k in range(5)],
        "triangle antiderivative period gate failed",
    )
    expected_means = {
        4: F(2187, 480200),
        8: F(927, 240100),
        12: F(1521, 480200),
    }
    expected_heights = {4: 197, 8: 232, 12: 268}
    expected_previous = {
        4: ((196, 197, 203, 207), F(-28297, 11357603964)),
        8: ((231, 232, 233, 241), F(-500165, 73729113612)),
        12: ((267, 268, 269, 279), F(-54, 1604047)),
    }

    rows = []
    audit_gates = 0
    for gap in (4, 8, 12):
        mean = averaged_overlap(gap)
        require(mean == expected_means[gap], f"mean changed g={gap}")
        forms, gates = residue_forms(gap, mean)
        audit_gates += gates
        height = expected_heights[gap]
        least_gate, least_start, least_odd, least_even = audit_uniform_height(
            gap, height, mean, forms
        )
        hostile_bank, hostile_start, hostile_pair, hostile_margin = (
            first_height_hostile(gap, height - 1)
        )
        require(
            (hostile_bank, hostile_margin) == expected_previous[gap],
            f"previous-height hostile changed g={gap}",
        )
        rows.append((
            gap, mean, height, least_gate, least_start, least_odd, least_even,
            hostile_bank, hostile_start, hostile_pair, hostile_margin,
        ))

    # Independent asymptotic screen via the closed triangle antiderivative.
    asymptotic_tax = F(5, 49 * 263) + F(6, 49 * 264)
    candidates = [
        gap for gap in range(2, 202, 2)
        if averaged_overlap(gap) / 2 > asymptotic_tax
    ]
    require(
        candidates == [4, 8, 12, 16, 20, 34, 38, 42, 46],
        "height-263 asymptotic screen changed",
    )
    require(
        F(4, 49) * LENGTH + F(2, 49 * 202) < 2 * asymptotic_tax,
        "large-gap screen failed",
    )
    screen_hostiles = []
    for gap in candidates:
        if gap in (4, 8):
            continue
        bank, start, pair, margin = first_height_hostile(gap, 263)
        require(margin < 0, f"non-strict screen hostile g={gap}")
        screen_hostiles.append((gap, bank, start, pair, margin))

    canonical_hostile = (8, 9, 11, 13)
    canonical_union = bank_union_measure(canonical_hostile)
    require(canonical_union == LENGTH, "canonical AP7 hostile lost full measure")
    canonical_pair_rows = []
    for start, gap, third in ((9, 2, 13), (11, 2, 9), (9, 4, 11)):
        pair = wall_pair_measure(start, gap)
        margin = certificate_margin(pair, start, gap, third, 8)
        require(margin < 0, "canonical hostile unexpectedly certified")
        canonical_pair_rows.append((start, gap, pair, margin))

    print("AP7 even-gap overlap atlas independent exact audit: PASS")
    print(f"J=[{LEFT},{RIGHT}] L={LENGTH}")
    print(f"endpoint/quasipolynomial_gates={audit_gates}")
    print("single-packet uniform certificate thresholds:")
    for row in rows:
        print(
            f"  g={row[0]} Q={row[1]} first_height={row[2]} "
            f"least_boundary_gate={row[3]} at_u={row[4]} "
            f"least_odd={row[5]} least_even={row[6]} "
            f"previous_hostile={row[7]} selected_u={row[8]} "
            f"pair={row[9]} margin={row[10]}"
        )
    print(f"height263_asymptotic_candidates={candidates}")
    print("height263 independently searched certificate hostiles:")
    for row in screen_hostiles:
        print(
            f"  g={row[0]} bank={row[1]} selected_u={row[2]} "
            f"pair={row[3]} margin={row[4]}"
        )
    print(
        f"canonical_hostile={canonical_hostile} "
        f"union_measure={canonical_union} selected_pair_rows={canonical_pair_rows}"
    )
    print("scope=thresholds are sharp for the selected-pair second-moment certificate, not actual absorption")


if __name__ == "__main__":
    main()
