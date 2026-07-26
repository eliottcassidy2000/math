#!/usr/bin/env python3
"""Exact companion for THM-2441.

The theorem compares the seven-colour inverse-branch ancestry profiles at
two clocks whose reduced clock parameters agree modulo ``7 D_0``.  The
deterministic core exhausts small grid sources and all single intervals on
three 13-divisible grids.  A separate power-clock core visits every residue
class in the exact ``ord_{7 D_0}(13)`` cycle.  Random fixtures are secondary.

Every calculation uses ``Fraction`` arithmetic.  All load-bearing checks use
``require`` so that ``python -O`` runs the same audit.
"""

from fractions import Fraction as F
from math import floor, gcd
from random import Random


P = 7
Interval = tuple[F, F]
IntervalSet = tuple[Interval, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_fraction(value: F) -> int:
    return -floor(-value)


def validate(intervals: IntervalSet, name: str) -> None:
    previous = F(0)
    for index, (left, right) in enumerate(intervals):
        require(F(0) <= left < right <= F(1), f"{name}: bad interval")
        if index:
            require(previous <= left, f"{name}: overlapping intervals")
        previous = right


def contains(intervals: IntervalSet, point: F) -> bool:
    return any(left <= point < right for left, right in intervals)


def measure(intervals: IntervalSet) -> F:
    return sum((right - left for left, right in intervals), F(0))


def count_residue_range(left: int, right: int, residue: int) -> int:
    if left > right:
        return 0
    first = left + (residue - left) % P
    if first > right:
        return 0
    return 1 + (right - first) // P


def arithmetic_profile(
    source: IntervalSet,
    clock: int,
    terminal: F,
) -> tuple[int, ...]:
    """Evaluate ancestry by integer-range residue counts."""

    require(clock >= 1, "positive clock")
    require(F(0) <= terminal < F(1), "terminal coordinate")
    counts = [0] * P
    for left, right in source:
        first = max(0, ceil_fraction(clock * left - terminal))
        last = min(clock - 1, ceil_fraction(clock * right - terminal) - 1)
        for residue in range(P):
            counts[residue] += count_residue_range(first, last, residue)
    return tuple(counts)


def brute_profile(
    source: IntervalSet,
    clock: int,
    terminal: F,
) -> tuple[int, ...]:
    """Definition-level evaluator used at modest clocks."""

    counts = [0] * P
    for prefix in range(clock):
        physical = (terminal + prefix) / clock
        if contains(source, physical):
            counts[prefix % P] += 1
    return tuple(counts)


def endpoint_word(
    source: IntervalSet,
    clock: int,
) -> tuple[tuple[F, int, int], ...]:
    """Signed endpoint multiset, including events at the terminal cut."""

    word: list[tuple[F, int, int]] = []
    for left, right in source:
        for endpoint, sign in ((left, 1), (right, -1)):
            scaled = clock * endpoint
            phase = scaled - floor(scaled)
            word.append((phase, floor(scaled) % P, sign))
    return tuple(sorted(word))


def partition_points(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, ...]:
    points = {F(0), F(1)}
    points.update(left for left, _ in target)
    points.update(right for _, right in target)
    points.update(phase for phase, _, _ in endpoint_word(source, clock))
    return tuple(sorted(points))


def energy(profile: tuple[int, ...]) -> int:
    return sum(
        (profile[(residue + 1) % P] - profile[residue]) ** 2
        for residue in range(P)
    )


def gamma(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> F:
    total = F(0)
    points = partition_points(source, target, clock)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if contains(target, midpoint):
            total += (right - left) * energy(
                arithmetic_profile(source, clock, midpoint)
            )
    return total


def grid_mask_source(denominator: int, mask: int) -> IntervalSet:
    """Turn a union of denominator-grid cells into maximal intervals."""

    intervals: list[Interval] = []
    index = 0
    while index < denominator:
        if not (mask >> index) & 1:
            index += 1
            continue
        start = index
        while index < denominator and (mask >> index) & 1:
            index += 1
        intervals.append((F(start, denominator), F(index, denominator)))
    return tuple(intervals)


def grid_cell_source(denominator: int, selected: list[bool]) -> IntervalSet:
    mask = sum((1 << index) for index, value in enumerate(selected) if value)
    return grid_mask_source(denominator, mask)


def multiplicative_order_13(modulus: int) -> int:
    require(gcd(13, modulus) == 1, "13 must be a unit")
    value = 13 % modulus
    order = 1
    while value != 1:
        value = value * 13 % modulus
        order += 1
        require(order <= modulus, "order search failed")
    return order


class Counters:
    pair_fixtures = 0
    chamber_checks = 0
    brute_profile_checks = 0
    endpoint_word_checks = 0
    gamma_checks = 0


def compare_congruent_clocks(
    source: IntervalSet,
    target: IntervalSet,
    *,
    K: int,
    D0: int,
    N: int,
    t: int,
    numerator_mass: int,
    brute: bool,
) -> None:
    """Audit the pointwise uniform-drift identity for one congruent pair."""

    require(K >= 0 and D0 >= 1 and gcd(D0, 13) == 1, "denominator split")
    require(N >= 1 and t >= 1, "clock parameters")
    D = (13**K) * D0
    require(measure(source) * D == numerator_mass, "numerator mass")
    next_N = N + P * D0 * t
    clock = (13**K) * N
    next_clock = (13**K) * next_N
    increment = t * numerator_mass

    word = endpoint_word(source, clock)
    next_word = endpoint_word(source, next_clock)
    require(word == next_word, "signed endpoint word changed")
    Counters.endpoint_word_checks += 1

    points = tuple(
        sorted(
            {
                *partition_points(source, target, clock),
                *partition_points(source, target, next_clock),
            }
        )
    )
    require(
        partition_points(source, target, clock)
        == partition_points(source, target, next_clock),
        "relative chamber partition changed",
    )
    probes = {F(0)}
    probes.update(point for point in points if point < 1)
    probes.update((left + right) / 2 for left, right in zip(points, points[1:]))

    source_mass = measure(source)
    for terminal in sorted(probes):
        profile = arithmetic_profile(source, clock, terminal)
        next_profile = arithmetic_profile(source, next_clock, terminal)
        require(
            next_profile
            == tuple(value + increment for value in profile),
            "pointwise uniform drift failed",
        )

        centered = tuple(
            F(value) - F(clock) * source_mass / P for value in profile
        )
        next_centered = tuple(
            F(value) - F(next_clock) * source_mass / P
            for value in next_profile
        )
        require(centered == next_centered, "global centering changed")

        point_mean = sum(profile, 0) / F(P)
        next_point_mean = sum(next_profile, 0) / F(P)
        require(
            tuple(F(value) - point_mean for value in profile)
            == tuple(F(value) - next_point_mean for value in next_profile),
            "pointwise projection changed",
        )

        if brute:
            require(
                profile == brute_profile(source, clock, terminal),
                "first arithmetic/brute mismatch",
            )
            require(
                next_profile == brute_profile(source, next_clock, terminal),
                "second arithmetic/brute mismatch",
            )
            Counters.brute_profile_checks += 2
        Counters.chamber_checks += 1

    require(
        gamma(source, target, clock) == gamma(source, target, next_clock),
        "target-restricted Gamma changed",
    )
    Counters.gamma_checks += 1
    Counters.pair_fixtures += 1


def deterministic_grid_core(target: IntervalSet) -> tuple[int, int]:
    """Exhaust all grid-cell sources through denominator eight."""

    sources = 0
    pairs = 0
    for D0 in range(1, 9):
        for mask in range(1 << D0):
            source = grid_mask_source(D0, mask)
            numerator_mass = mask.bit_count()
            sources += 1
            for N in range(1, 5):
                for t in (1, 2):
                    compare_congruent_clocks(
                        source,
                        target,
                        K=0,
                        D0=D0,
                        N=N,
                        t=t,
                        numerator_mass=numerator_mass,
                        brute=True,
                    )
                    pairs += 1
    return sources, pairs


def deterministic_thirteen_grid_core(target: IntervalSet) -> tuple[int, int]:
    """Exhaust every single interval on grids 13, 26, and 39."""

    intervals = 0
    pairs = 0
    for D0 in (1, 2, 3):
        D = 13 * D0
        for left in range(D):
            for right in range(left + 1, D + 1):
                source = ((F(left, D), F(right, D)),)
                intervals += 1
                for N in (1, 2):
                    compare_congruent_clocks(
                        source,
                        target,
                        K=1,
                        D0=D0,
                        N=N,
                        t=1,
                        numerator_mass=right - left,
                        brute=True,
                    )
                    pairs += 1
    return intervals, pairs


def deterministic_power_core(target: IntervalSet) -> tuple[int, int, int]:
    """Visit every phase class of every small exact power-clock cycle."""

    source_count = 0
    residue_classes = 0
    maximum_period = 0
    for K in (0, 1):
        for D0 in range(1, 13):
            D = (13**K) * D0
            selected = [
                (index % 5 in (0, 2)) or (index % 11 == 7)
                for index in range(D)
            ]
            source = grid_cell_source(D, selected)
            numerator_mass = sum(selected)
            period = multiplicative_order_13(P * D0)
            maximum_period = max(maximum_period, period)
            require(period % 2 == 0, "physical reflection not periodic")
            for d in range(period):
                N = 13**d
                next_N = 13 ** (d + period)
                difference = next_N - N
                require(difference % (P * D0) == 0, "power congruence")
                t = difference // (P * D0)
                compare_congruent_clocks(
                    source,
                    target,
                    K=K,
                    D0=D0,
                    N=N,
                    t=t,
                    numerator_mass=numerator_mass,
                    brute=False,
                )
                residue_classes += 1
            source_count += 1
    return source_count, residue_classes, maximum_period


def secondary_random_audit(target: IntervalSet) -> int:
    rng = Random(2_441)
    fixtures = 0
    for _ in range(96):
        K = rng.randrange(2)
        D0 = rng.randrange(1, 11)
        D = (13**K) * D0
        selected = [rng.randrange(2) == 1 for _ in range(D)]
        source = grid_cell_source(D, selected)
        compare_congruent_clocks(
            source,
            target,
            K=K,
            D0=D0,
            N=rng.randrange(1, 21),
            t=rng.randrange(1, 4),
            numerator_mass=sum(selected),
            brute=False,
        )
        fixtures += 1
    return fixtures


def hostile_and_positive_controls() -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    F,
    tuple[F, ...],
]:
    full = ((F(0), F(1)),)
    hostile_target = ((F(1, 3), F(2, 3)),)
    small_profile = arithmetic_profile(full, 13, F(1, 2))
    large_profile = arithmetic_profile(full, 13**3, F(1, 2))
    require(small_profile == (2, 2, 2, 2, 2, 2, 1), "small hostile")
    require(
        large_profile == (314, 314, 314, 314, 314, 314, 313),
        "large hostile",
    )
    require(small_profile != large_profile, "raw word should not be periodic")
    require(
        tuple(F(value) - F(13, 7) for value in small_profile)
        == tuple(F(value) - F(13**3, 7) for value in large_profile),
        "hostile centered profiles",
    )
    hostile_gamma = gamma(full, hostile_target, 13)
    require(hostile_gamma == F(2, 3), "hostile Gamma value")
    require(
        gamma(full, hostile_target, 13**3) == hostile_gamma,
        "hostile Gamma period",
    )

    alternating_source = (
        (F(1, 39), F(8, 39)),
        (F(13, 39), F(17, 39)),
    )
    alternating_target = ((F(1, 5), F(4, 5)),)
    alternating = tuple(
        gamma(alternating_source, alternating_target, 13 ** (1 + d))
        for d in range(6)
    )
    require(
        alternating
        == (
            F(12, 5),
            F(16, 15),
            F(12, 5),
            F(16, 15),
            F(12, 5),
            F(16, 15),
        ),
        "period-two positive control",
    )
    require(multiplicative_order_13(21) == 2, "order modulo 21")
    return small_profile, large_profile, hostile_gamma, alternating


def main() -> None:
    target = (
        (F(1, 11), F(4, 11)),
        (F(7, 11), F(10, 11)),
    )
    validate(target, "target")

    grid_sources, grid_pairs = deterministic_grid_core(target)
    interval_sources, interval_pairs = deterministic_thirteen_grid_core(target)
    power_sources, power_classes, maximum_period = deterministic_power_core(target)
    random_fixtures = secondary_random_audit(target)
    small, large, hostile_gamma, alternating = hostile_and_positive_controls()

    print("THM-2441 exact ancestry event-period companion")
    print(
        "deterministic_grid_core:",
        f"sources={grid_sources}",
        f"congruent_pairs={grid_pairs}",
    )
    print(
        "deterministic_13_grid_core:",
        f"single_intervals={interval_sources}",
        f"congruent_pairs={interval_pairs}",
    )
    print(
        "deterministic_power_core:",
        f"sources={power_sources}",
        f"residue_classes={power_classes}",
        f"maximum_period={maximum_period}",
    )
    print("secondary_random_fixtures:", random_fixtures)
    print(
        "identity_checks:",
        f"pairs={Counters.pair_fixtures}",
        f"chambers={Counters.chamber_checks}",
        f"brute_profiles={Counters.brute_profile_checks}",
        f"endpoint_words={Counters.endpoint_word_checks}",
        f"gamma={Counters.gamma_checks}",
    )
    print("raw_hostile_profiles:", small, large)
    print("raw_hostile_gamma:", hostile_gamma)
    print("period_two_positive_gamma:", alternating)
    print(
        "VERIFIED: congruent reduced clocks differ by one exact uniform "
        "ancestry vector; centered event words and target Gamma are periodic."
    )


if __name__ == "__main__":
    main()
