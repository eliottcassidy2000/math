#!/usr/bin/env python3
"""Dependency-free exact hostile audit for THM-2421.

The two evaluators below are intentionally different.

* ``brute_profile`` tests every inverse branch ``(y+n)/R``.
* ``event_sweep`` starts from an interval-counting formula and updates the
  seven ancestry counts only at signed endpoint events.

All geometry and all integrals use ``Fraction`` arithmetic.  In particular,
the large-clock tests do not enumerate the ``R=13^k`` inverse branches.
"""

from fractions import Fraction as F
from math import floor
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


def interval_measure(intervals: IntervalSet) -> F:
    return sum((right - left for left, right in intervals), F(0))


def count_residue_range(left: int, right: int, residue: int) -> int:
    """Count n in the inclusive integer interval [left,right], n=residue mod 7."""

    if left > right:
        return 0
    first = left + (residue - left) % P
    if first > right:
        return 0
    return 1 + (right - first) // P


def arithmetic_profile(intervals: IntervalSet, clock: int, terminal: F) -> tuple[int, ...]:
    """Count inverse branches by integer-range arithmetic, without an R-loop."""

    require(clock >= 1, "clock")
    require(F(0) < terminal < F(1), "terminal must lie in an open chamber")
    counts = [0] * P
    for left, right in intervals:
        first = max(0, ceil_fraction(clock * left - terminal))
        last = min(clock - 1, ceil_fraction(clock * right - terminal) - 1)
        for residue in range(P):
            counts[residue] += count_residue_range(first, last, residue)
    return tuple(counts)


def brute_profile(intervals: IntervalSet, clock: int, terminal: F) -> tuple[int, ...]:
    """Independent O(R) definition-level evaluator."""

    counts = [0] * P
    for prefix in range(clock):
        physical = (terminal + prefix) / clock
        if contains(intervals, physical):
            counts[prefix % P] += 1
    return tuple(counts)


def signed_endpoint_events(
    intervals: IntervalSet,
    clock: int,
) -> dict[F, tuple[int, ...]]:
    """Interior event y -> aggregate signed jump in the seven ancestry bins."""

    raw: dict[F, list[int]] = {}
    for left, right in intervals:
        for endpoint, sign in ((left, 1), (right, -1)):
            scaled = clock * endpoint
            phase = scaled - floor(scaled)
            # Grid endpoints occur at the cut 0=1 of the terminal coordinate.
            # They affect only null boundary values, not any open chamber.
            if phase == 0:
                continue
            jump = raw.setdefault(phase, [0] * P)
            jump[floor(scaled) % P] += sign
    return {phase: tuple(jump) for phase, jump in raw.items()}


def energy(profile: tuple[int, ...]) -> int:
    return sum(
        (profile[(residue + 1) % P] - profile[residue]) ** 2
        for residue in range(P)
    )


def flat(profile: tuple[int, ...]) -> bool:
    return all(value == profile[0] for value in profile)


def partition_points(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, ...]:
    events = signed_endpoint_events(source, clock)
    points = {F(0), F(1), *events}
    for left, right in target:
        points.add(left)
        points.add(right)
    return tuple(sorted(points))


def brute_partition_points(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, ...]:
    """Definition-level breakpoint search over every endpoint/branch pair."""

    points = {F(0), F(1)}
    for left, right in target:
        points.add(left)
        points.add(right)
    for left, right in source:
        for endpoint in (left, right):
            for prefix in range(clock):
                phase = clock * endpoint - prefix
                if F(0) < phase < F(1):
                    points.add(phase)
    return tuple(sorted(points))


def event_sweep(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, F, int, int, tuple[tuple[int, ...], ...]]:
    """Return Gamma, nonflat target mass, chamber count, jump checks, profiles."""

    validate(source, "source")
    validate(target, "target")
    points = partition_points(source, target, clock)
    events = signed_endpoint_events(source, clock)
    require(len(points) >= 2, "empty partition")

    profiles: list[tuple[int, ...]] = []
    gamma = F(0)
    nonflat_mass = F(0)
    jump_checks = 0
    current: tuple[int, ...] | None = None

    for index, (left, right) in enumerate(zip(points, points[1:])):
        require(left < right, "degenerate chamber")
        midpoint = (left + right) / 2
        direct = arithmetic_profile(source, clock, midpoint)

        if index == 0:
            current = direct
        else:
            jump = events.get(left, (0,) * P)
            current = tuple(
                value + delta for value, delta in zip(current, jump)
            )
            if any(jump):
                jump_checks += 1
            require(current == direct, "signed endpoint jump formula failed")

        profiles.append(current)
        if contains(target, midpoint):
            length = right - left
            gamma += length * energy(current)
            if not flat(current):
                nonflat_mass += length

    require((gamma > 0) == (nonflat_mass > 0), "detector equivalence failed")
    return gamma, nonflat_mass, len(profiles), jump_checks, tuple(profiles)


def brute_integral(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, F, int]:
    """Definition-level chamber integration; used only at modest clocks."""

    points = brute_partition_points(source, target, clock)
    gamma = F(0)
    nonflat_mass = F(0)
    profile_checks = 0
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        brute = brute_profile(source, clock, midpoint)
        arithmetic = arithmetic_profile(source, clock, midpoint)
        require(brute == arithmetic, "brute/arithmetic profile mismatch")
        profile_checks += 1
        if contains(target, midpoint):
            length = right - left
            gamma += length * energy(brute)
            if not flat(brute):
                nonflat_mass += length
    return gamma, nonflat_mass, profile_checks


def integrated_profile(
    source: IntervalSet,
    target: IntervalSet,
    clock: int,
) -> tuple[F, ...]:
    """Compute integral_Q N_c(y)dy exactly on the event partition."""

    totals = [F(0)] * P
    points = partition_points(source, target, clock)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if contains(target, midpoint):
            profile = arithmetic_profile(source, clock, midpoint)
            for residue, value in enumerate(profile):
                totals[residue] += (right - left) * value
    return tuple(totals)


def full_flat_certificate(source: IntervalSet, clock: int) -> bool:
    """Endpoint rigidity certificate for Gamma_[0,1)=0."""

    events = signed_endpoint_events(source, clock)
    points = tuple(sorted({F(0), F(1), *events}))
    midpoint = (points[0] + points[1]) / 2
    return flat(arithmetic_profile(source, clock, midpoint)) and all(
        flat(jump) for jump in events.values()
    )


def random_interval_set(rng: Random, denominator: int, pairs: int) -> IntervalSet:
    cuts = sorted(rng.sample(range(1, denominator), 2 * pairs))
    return tuple(
        (F(cuts[2 * index], denominator), F(cuts[2 * index + 1], denominator))
        for index in range(pairs)
    )


def danger_intervals(speed: int) -> IntervalSet:
    """Exact D_speed={y: ||speed*y||<1/14}, split at the circle cut."""

    radius = F(1, 14 * speed)
    pieces: list[Interval] = [(F(0), radius), (F(1) - radius, F(1))]
    for center_index in range(1, speed):
        center = F(center_index, speed)
        pieces.append((center - radius, center + radius))
    return tuple(sorted(pieces))


def main() -> None:
    # Random exact hostile audit: every chamber is evaluated by the O(R)
    # definition, the integer-range formula, and the signed-event sweep.
    rng = Random(2_421)
    fixture_count = 0
    chamber_checks = 0
    jump_checks = 0
    flat_certificate_checks = 0
    full = ((F(0), F(1)),)
    for clock in (2, 3, 6, 7, 13, 17, 29, 169):
        for trial in range(24):
            source = random_interval_set(rng, 31 + 2 * (trial % 4), 3)
            target = random_interval_set(rng, 37 + 2 * (trial % 5), 2)
            swept = event_sweep(source, target, clock)
            brute = brute_integral(source, target, clock)
            require(swept[:2] == brute[:2], "integral evaluator mismatch")
            require(swept[2] == brute[2], "chamber-count mismatch")
            full_sweep = event_sweep(source, full, clock)
            require(
                (full_sweep[0] == 0) == full_flat_certificate(source, clock),
                "full-target endpoint rigidity failed",
            )
            fixture_count += 1
            chamber_checks += brute[2]
            jump_checks += swept[3]
            flat_certificate_checks += 1

    # Boundary and coincident-event controls.  At y=1/2 the first source
    # interval exits ancestry bin 0 exactly when the second enters bin 1.
    coincident_source = (
        (F(0), F(1, 26)),
        (F(3, 26), F(2, 13)),
    )
    coincident_events = signed_endpoint_events(coincident_source, 13)
    require(coincident_events == {F(1, 2): (-1, 1, 0, 0, 0, 0, 0)}, "coincident event")
    coincident = event_sweep(coincident_source, full, 13)
    coincident_brute = brute_integral(coincident_source, full, 13)
    require(coincident[:2] == coincident_brute[:2] == (F(2), F(1)), "coincident integral")

    # A source crossing the circle cut checks that grid events may be omitted
    # from the open-chamber sweep while exact initialization retains them.
    wrap_source = ((F(0), F(1, 10)), (F(9, 10), F(1)))
    wrap_target = ((F(0), F(2, 5)), (F(3, 5), F(1)))
    wrap_checks = 0
    for clock in (6, 7, 13, 29, 169):
        swept = event_sweep(wrap_source, wrap_target, clock)
        brute = brute_integral(wrap_source, wrap_target, clock)
        require(swept[:3] == brute[:2] + (brute[2],), "wrap control")
        wrap_checks += brute[2]

    # THM-2418's fixed real/even BV-two rank-one hostile is exactly flat at
    # every 13-adic clock.  The k=8 check has R=815,730,721 but only one
    # chamber and no prefix loop.
    rank_one_source = ((F(3, 13), F(10, 13)),)
    rank_one_clocks = []
    for depth in range(1, 9):
        clock = 13**depth
        gamma, nonflat_mass, chambers, jumps, profiles = event_sweep(
            rank_one_source,
            full,
            clock,
        )
        expected = (13 ** (depth - 1),) * P
        require(gamma == nonflat_mass == 0, "rank-one detector hostile")
        require(chambers == 1 and jumps == 0, "rank-one event shape")
        require(profiles == (expected,), "rank-one ancestry profile")
        rank_one_clocks.append(clock)

    # Complete one-interval classification.  The profile is flat exactly
    # when the R-scaled interval length is an integral multiple of seven.
    single_interval_checks = 0
    for clock in (2, 6, 7, 13, 14, 29):
        for left_numerator in range(29):
            for right_numerator in range(left_numerator + 1, 30):
                source = (
                    (F(left_numerator, 29), F(right_numerator, 29)),
                )
                scaled_length = clock * (source[0][1] - source[0][0])
                predicted_flat = (
                    scaled_length.denominator == 1
                    and scaled_length.numerator % P == 0
                )
                actual_flat = event_sweep(source, full, clock)[0] == 0
                require(actual_flat == predicted_flat, "single interval classifier")
                single_interval_checks += 1

    # The sevenfold threshold is sharp: simultaneous +1 (and then -1)
    # events in all residues change the common ancestry level but keep the
    # profile flat on every chamber.
    seven_sheet_source = tuple(
        (F(4 * prefix + 1, 52), F(4 * prefix + 3, 52))
        for prefix in range(P)
    )
    seven_sheet_events = signed_endpoint_events(seven_sheet_source, 13)
    require(
        seven_sheet_events
        == {F(1, 4): (1,) * P, F(3, 4): (-1,) * P},
        "seven-sheet event packets",
    )
    seven_sheet = event_sweep(seven_sheet_source, full, 13)
    require(seven_sheet[:2] == (F(0), F(0)), "seven-sheet flat hostile")
    require(
        seven_sheet[4] == ((0,) * P, (1,) * P, (0,) * P),
        "seven-sheet profiles",
    )

    # Positive and false-sufficient-test controls at a large clock.  Prefix
    # cylinders have grid endpoints, so their profiles are constant in y.
    large_clock = 13**6
    one_cylinder = ((F(0), F(1, large_clock)),)
    one = event_sweep(one_cylinder, full, large_clock)
    require(one[:2] == (F(2), F(1)), "one-cylinder positive control")
    require(one[4] == ((1, 0, 0, 0, 0, 0, 0),), "one-cylinder profile")

    # Seven retained cells, but residue counts (2,1,1,1,1,1,0).  Hence total
    # ancestry is 0 modulo 7 while the full vector is nonflat and Gamma=6.
    divisible_nonflat = (
        (F(0), F(6, large_clock)),
        (F(7, large_clock), F(8, large_clock)),
    )
    divisible = event_sweep(divisible_nonflat, full, large_clock)
    require(divisible[:2] == (F(6), F(1)), "divisible-total control")
    require(divisible[4] == ((2, 1, 1, 1, 1, 1, 0),), "divisible profile")
    require(sum(divisible[4][0]) % P == 0, "total-mod-seven hostile")

    # A stronger cancellation hostile for the fixed-pair classifier of
    # THM-2418.  Every integrated carry mass is 1/2, although the pointwise
    # ancestry profile is (1,1,1,0,0,0,0) on the first half and its
    # complementary rotation on the second half.
    balanced_rotation_source = tuple(
        (
            (F(prefix, 13), F(2 * prefix + 1, 26))
            if prefix < 3
            else (F(2 * prefix + 1, 26), F(prefix + 1, 13))
        )
        for prefix in range(P)
    )
    balanced = event_sweep(balanced_rotation_source, full, 13)
    balanced_integral = integrated_profile(balanced_rotation_source, full, 13)
    require(balanced[:2] == (F(2), F(1)), "balanced-rotation energy")
    require(balanced_integral == (F(1, 2),) * P, "balanced integrated profile")

    # Exact terminal restriction.  D_7 has mass 1/7; it leaves the flat
    # hostile flat and scales the one-cylinder energy from 2 to 2/7.
    terminal_d7 = danger_intervals(7)
    require(interval_measure(terminal_d7) == F(1, 7), "D7 measure")
    d7_flat = event_sweep(rank_one_source, terminal_d7, 13**5)
    d7_positive = event_sweep(one_cylinder, terminal_d7, large_clock)
    require(d7_flat[:2] == (F(0), F(0)), "D7 flat control")
    require(d7_positive[:2] == (F(2, 7), F(1, 7)), "D7 positive control")

    print("theorem=THM-2421")
    print("status=PROVED+VERIFIED-EXACT+HOSTILE-AUDITED; independent-audit=PENDING")
    print(f"random_fixtures={fixture_count}")
    print(f"random_chamber_checks={chamber_checks}")
    print(f"random_nonzero_jump_checks={jump_checks}")
    print(f"full_flat_certificate_checks={flat_certificate_checks}")
    print(f"wrap_chamber_checks={wrap_checks}")
    print("coincident_event=1/2:(-1,1,0,0,0,0,0)")
    print(f"coincident_gamma={coincident[0]}; coincident_nonflat_mass={coincident[1]}")
    print(
        "rank_one_clocks="
        + ",".join(map(str, rank_one_clocks))
        + f"; max_clock={rank_one_clocks[-1]}"
    )
    print("rank_one_gamma=0; rank_one_nonflat_mass=0")
    print(f"single_interval_classifier_checks={single_interval_checks}")
    print("seven_sheet_events=1/4:+ones,3/4:-ones; gamma=0")
    print(f"one_cylinder_clock={large_clock}; gamma={one[0]}; nonflat_mass={one[1]}")
    print(
        "divisible_total_profile="
        + ",".join(map(str, divisible[4][0]))
        + f"; total=7; gamma={divisible[0]}"
    )
    print(
        "balanced_integrated_profile="
        + ",".join(map(str, balanced_integral))
        + f"; pointwise_gamma={balanced[0]}"
    )
    print(f"D7_mass={interval_measure(terminal_d7)}; D7_positive_gamma={d7_positive[0]}")
    print("scalar_cover_Ej_Qjsigma_application=OPEN_NO_EXPLICIT_GLOBAL_PACKET")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
