#!/usr/bin/env python3
"""Independent exact hostile for the sharp arbitrary-language boundary of THM-3499.

The program does not import either THM-3499 companion.  It checks the binary
shortlex endpoints and half-block counts directly, encloses the two level
harmonic limits by exact rational logarithm intervals, certifies separated
complete-stage normalized intervals for the alternating construction, and
types the depth-two affine F_2^2/K4/matching/tournament bookkeeping.
"""

from __future__ import annotations

import hashlib
import itertools
from decimal import Decimal, localcontext
from fractions import Fraction


EXPECTED_LEDGER_SHA256 = "e0f67a9bdda9e6e158969c54c023b390c17fed8c4ce68b0ddeffee4936824513"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def log_ratio_interval(
    numerator: int, denominator: int, terms: int = 32
) -> tuple[Fraction, Fraction]:
    """Exact enclosure of log(numerator/denominator) for numerator>denominator.

    With z=(a-b)/(a+b), use
      log(a/b)=2*sum_{j>=0} z^(2j+1)/(2j+1).
    The omitted positive tail is bounded geometrically after replacing every
    remaining denominator by the first omitted denominator.
    """

    require(numerator > denominator > 0, "log-ratio domain failure")
    require(terms >= 1, "log-ratio truncation must be positive")
    z = Fraction(numerator - denominator, numerator + denominator)
    partial = sum(
        (2 * z ** (2 * index + 1) / (2 * index + 1) for index in range(terms)),
        Fraction(0),
    )
    first_power = 2 * terms + 1
    tail = 2 * z**first_power / (first_power * (1 - z * z))
    return partial, partial + tail


def decimal_text(value: Fraction, places: int = 15) -> str:
    with localcontext() as context:
        context.prec = places + 20
        decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        return f"{decimal_value:.{places}f}"


def exact_harmonic_interval(start: int, stop: int) -> Fraction:
    require(1 <= start <= stop, "invalid harmonic interval")
    return sum((Fraction(1, value) for value in range(start, stop + 1)), Fraction(0))


def normalized_interval(
    numerator_interval: tuple[Fraction, Fraction],
    denominator_interval: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    numerator_low, numerator_high = numerator_interval
    denominator_low, denominator_high = denominator_interval
    require(numerator_low > 0 and denominator_low > 0, "nonpositive ratio interval")
    return numerator_low / denominator_high, numerator_high / denominator_low


def maximum_endpoint_distance(
    first: tuple[Fraction, Fraction], second: tuple[Fraction, Fraction]
) -> Fraction:
    return max(
        abs(first[0] - second[0]),
        abs(first[0] - second[1]),
        abs(first[1] - second[0]),
        abs(first[1] - second[1]),
    )


def canonical_edge(left: tuple[int, int], right: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return (left, right) if left < right else (right, left)


def main() -> None:
    log_three_halves = log_ratio_interval(3, 2)
    log_four_thirds = log_ratio_interval(4, 3)
    log_two = log_ratio_interval(2, 1)
    left_target = normalized_interval(log_three_halves, log_two)
    right_target = normalized_interval(log_four_thirds, log_two)

    require(left_target[0] > right_target[1], "the two target coefficients overlap")
    require(
        left_target[0] + right_target[0] <= 1 <= left_target[1] + right_target[1],
        "complementary target intervals do not enclose one",
    )

    ledger: list[str] = []
    for level in range(1, 10):
        level_start = 1 << level
        half_size = 1 << (level - 1)
        first_start = level_start
        first_stop = level_start + half_size - 1
        last_start = first_stop + 1
        last_stop = (1 << (level + 1)) - 1

        require(level_start == 1 << level, "binary shortlex level start failure")
        require(last_stop - level_start + 1 == 1 << level, "binary level size failure")
        require(first_stop - first_start + 1 == half_size, "first-half count failure")
        require(last_stop - last_start + 1 == half_size, "last-half count failure")
        require(first_stop + 1 == last_start, "half-block seam failure")

        first_mass = exact_harmonic_interval(first_start, first_stop)
        last_mass = exact_harmonic_interval(last_start, last_stop)
        first_error_bound = Fraction(1, 3 * (1 << level))
        last_error_bound = Fraction(1, 6 * (1 << level))
        require(first_mass >= log_three_halves[0], "first mass below its integral")
        require(
            first_mass <= log_three_halves[1] + first_error_bound,
            "first-mass rectangle error failure",
        )
        require(last_mass >= log_four_thirds[0], "last mass below its integral")
        require(
            last_mass <= log_four_thirds[1] + last_error_bound,
            "last-mass rectangle error failure",
        )
        ledger.append(
            f"level={level};interval={level_start}:{last_stop};"
            f"first={first_start}:{first_stop};last={last_start}:{last_stop};"
            f"masses={first_mass},{last_mass}"
        )

    # Stages start at k=0 and cover every positive level.  Even stages take
    # the first half, odd stages the last half.  For the completed stage with
    # T levels, N=2^(T+1)-1, hence T*log(2)<log(N)<(T+1)*log(2).
    first_level_count = 0
    last_level_count = 0
    total_levels = 0
    stage_intervals: dict[int, tuple[Fraction, Fraction]] = {}
    stage_rows: list[tuple[int, int, int, str, tuple[Fraction, Fraction]]] = []
    for stage in range(9):
        stage_length = 1 << (1 << stage)
        if stage % 2 == 0:
            first_level_count += stage_length
            choice = "first"
        else:
            last_level_count += stage_length
            choice = "last"
        total_levels += stage_length
        require(first_level_count + last_level_count == total_levels, "stage count failure")

        numerator_low = (
            first_level_count * log_three_halves[0]
            + last_level_count * log_four_thirds[0]
        )
        # Every chosen level is a left Riemann sum.  The sum of all possible
        # first-half errors over positive levels is at most 1/3; this also
        # bounds any mixture of first- and last-half choices.
        numerator_high = (
            first_level_count * log_three_halves[1]
            + last_level_count * log_four_thirds[1]
            + Fraction(1, 3)
        )
        denominator_low = total_levels * log_two[0]
        denominator_high = (total_levels + 1) * log_two[1]
        ratio_interval = normalized_interval(
            (numerator_low, numerator_high), (denominator_low, denominator_high)
        )
        stage_intervals[stage] = ratio_interval
        stage_rows.append((stage, stage_length, total_levels, choice, ratio_interval))
        ledger.append(
            f"stage={stage};length={stage_length};levels={total_levels};"
            f"counts={first_level_count},{last_level_count};ratio={ratio_interval[0]},{ratio_interval[1]}"
        )

    require(stage_intervals[8][0] > stage_intervals[7][1], "late stage intervals are not disjoint")
    require(
        maximum_endpoint_distance(stage_intervals[7], right_target) < Fraction(1, 10**15),
        "odd-stage interval is not close to the right target",
    )
    require(
        maximum_endpoint_distance(stage_intervals[8], left_target) < Fraction(1, 10**30),
        "even-stage interval is not close to the left target",
    )

    # Depth two is F_2^2.  Lexicographic comparison orients K4 transitively.
    # The three nonzero covectors split the six K4 edges into three matchings.
    vertices = tuple(itertools.product((0, 1), repeat=2))
    edges = tuple(itertools.combinations(vertices, 2))
    arcs = tuple((left, right) for left, right in edges if left < right)
    require(len(vertices) == 4, "depth-two vertex count failure")
    require(len(edges) == 6 and len(arcs) == 6, "K4 edge/arc count failure")
    require(len(tuple(itertools.combinations(edges, 2))) == 15, "six-object comparison count failure")

    forms = ((1, 0), (0, 1), (1, 1))
    matchings: dict[tuple[int, int], tuple[tuple[tuple[int, int], tuple[int, int]], ...]] = {}
    covered_edges: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for form in forms:
        fibres = []
        for value in (0, 1):
            fibre = tuple(
                vertex
                for vertex in vertices
                if (form[0] * vertex[0] + form[1] * vertex[1]) % 2 == value
            )
            require(len(fibre) == 2, "nonzero form does not give a 2+2 partition")
            fibres.append(fibre)
        matching = tuple(canonical_edge(*fibre) for fibre in fibres)
        require(len(set(sum((list(edge) for edge in matching), []))) == 4, "not a perfect matching")
        matchings[form] = matching
        covered_edges.extend(matching)

    require(len(set(covered_edges)) == 6, "three matchings do not partition E(K4)")
    require(set(covered_edges) == set(edges), "matching edge partition mismatch")
    for form in forms:
        ledger.append(f"form={form};matching={matchings[form]}")

    digest = hashlib.sha256("\n".join(ledger).encode("ascii")).hexdigest()
    require(
        digest == EXPECTED_LEDGER_SHA256,
        f"exact semantic ledger changed: {digest}",
    )

    print("THM-3510 independent exact hostile")
    print("binary levels n=1..9: endpoints, seams, equal half-counts, and rational masses exact")
    print(
        "normalized fixed targets: "
        f"first=[{decimal_text(left_target[0])},{decimal_text(left_target[1])}] "
        f"last=[{decimal_text(right_target[0])},{decimal_text(right_target[1])}]"
    )
    for stage, stage_length, levels, choice, interval in stage_rows[3:9]:
        print(
            f"stage k={stage} length={stage_length} levels={levels} choice={choice} "
            f"normalized=[{decimal_text(interval[0])},{decimal_text(interval[1])}]"
        )
    print("depth two: 4 word vertices, 6 K4 comparisons/arcs, 3 covector matchings")
    print("six edge-objects would require 15 pairwise comparisons, not 6")
    print(f"exact ledger sha256={digest}")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
