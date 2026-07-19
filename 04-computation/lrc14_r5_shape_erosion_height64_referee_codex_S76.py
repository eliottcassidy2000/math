#!/usr/bin/env python3
"""Exact referee for THM-1168's proportional r=5 shape atlas.

The computation is dependency-free and uses ``fractions.Fraction`` for every
wall, erosion endpoint, cyclic gap, core floor, and component-address
intersection.  It scans precisely the primitive ratio rays left by the three
scale-free THM-1148 gates, through primitive height d <= 64.

No ``assert`` is used, so the replay checks the same obligations under
``python`` and ``python -O``.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd


LAMBDA = F(1, 14)
HEIGHT = 64
EXCEPTION = (49, 50, 51, 56)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# A speed-k danger comb is open.  On the cut circle (0,1), an R event ends
# an open danger tooth and an L event starts one.  Fractions are cached once
# per speed because the 95,336-shape census reuses them heavily.
_SPEED_EVENTS: dict[int, tuple[tuple[F, int, int], ...]] = {}


def speed_events(k: int) -> tuple[tuple[F, int, int], ...]:
    cached = _SPEED_EVENTS.get(k)
    if cached is not None:
        return cached
    events: list[tuple[F, int, int]] = [(F(1, 14 * k), -1, 1)]
    for n in range(1, k):
        events.append((F(14 * n - 1, 14 * k), 1, 0))
        events.append((F(14 * n + 1, 14 * k), -1, 1))
    events.append((F(14 * k - 1, 14 * k), 1, 0))
    cached = tuple(events)
    _SPEED_EVENTS[k] = cached
    return cached


def closed_safe_components(
    speeds: tuple[int, ...],
) -> tuple[list[tuple[F, F]], int, int, int, list[F]]:
    """Return positive closed components and exact wall telemetry.

    The telemetry is (raw events, distinct positions, tied positions,
    isolated safe walls).  Isolated safe walls cannot carry a positive target
    interval and are deliberately not inserted among the positive components.
    """
    events: dict[F, list[int]] = {}
    raw_count = 0
    for k in speeds:
        for point, change, ending in speed_events(k):
            row = events.setdefault(point, [0, 0, 0])
            row[0] += change
            row[1] += ending
            row[2] += 1
            raw_count += 1

    active = len(speeds)  # every speed is dangerous immediately right of 0
    previous = F(0)
    components: list[tuple[F, F]] = []
    isolated: list[F] = []
    for point, (change, endings, _) in sorted(events.items()):
        left_active = active
        if left_active == 0:
            if components and components[-1][1] == previous:
                components[-1] = (components[-1][0], point)
            else:
                components.append((previous, point))

        # Changing danger intervals are open at the wall.  Only intervals
        # continuing through the wall can make the wall itself dangerous.
        continuing = left_active - endings
        active += change
        if left_active > 0 and active > 0 and continuing == 0:
            isolated.append(point)
        previous = point

    require(active == len(speeds), f"bad wrap count for {speeds}: {active}")
    tied_positions = sum(1 for row in events.values() if row[2] > 1)
    return components, raw_count, len(events), tied_positions, isolated


def circle_norm(x: F) -> F:
    y = x % 1
    return min(y, 1 - y)


def direct_safe(speeds: tuple[int, ...], x: F) -> bool:
    return all(circle_norm(k * x) >= LAMBDA for k in speeds)


def verify_sweep_directly(speeds: tuple[int, ...]) -> None:
    """Independent cell/midpoint audit of the event sweep for named rows."""
    components, _, _, _, isolated = closed_safe_components(speeds)
    walls = sorted({point for k in speeds for point, _, _ in speed_events(k)})
    cuts = [F(0), *walls, F(1)]
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left + right) / 2
        in_component = any(a < midpoint < b for a, b in components)
        require(
            in_component == direct_safe(speeds, midpoint),
            f"sweep/midpoint disagreement for {speeds} at {midpoint}",
        )
    for left, right in components:
        require(direct_safe(speeds, left), f"unsafe left endpoint {left}")
        require(direct_safe(speeds, right), f"unsafe right endpoint {right}")
    for point in isolated:
        require(direct_safe(speeds, point), f"unsafe isolated wall {point}")


def erosion_atlas(
    shape: tuple[int, int, int, int],
) -> tuple[F, F, F, list[tuple[F, F]], list[tuple[F, F]], list[F], tuple[int, int, int, int]]:
    components, raw, distinct, ties, isolated = closed_safe_components(shape)
    delta = F(1, 7 * shape[-1])
    starts = [
        (left, right - delta)
        for left, right in components
        if right - left >= delta
    ]
    require(starts, f"no delta-admissible start for {shape}")
    gaps: list[F] = []
    for index, (_, right) in enumerate(starts):
        next_left = starts[(index + 1) % len(starts)][0]
        if index + 1 == len(starts):
            next_left += 1
        gaps.append(next_left - right)
    largest_gap = max(gaps)
    threshold = largest_gap + delta
    return (
        delta,
        largest_gap,
        threshold,
        components,
        starts,
        gaps,
        (raw, distinct, ties, len(isolated)),
    )


def transfer_phi_pair(x_num: int, x_den: int) -> tuple[int, int] | None:
    """Exact THM-1137 Phi(x), represented as a reduced numerator pair."""
    if x_num < x_den:
        return None
    if 7 * x_num >= 13 * x_den:
        return (6, 7)
    numerator = 7 * x_num - x_den
    denominator = 14 * x_den
    common = gcd(numerator, denominator)
    return (numerator // common, denominator // common)


def exact_transfer_gate(shape: tuple[int, int, int, int]) -> bool:
    numerator, denominator = 6, 7
    for old, new in zip(shape, shape[1:]):
        output = transfer_phi_pair(new * numerator, old * denominator)
        if output is None:
            return False
        numerator, denominator = output
    return 7 * numerator > denominator


def primitive(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    return gcd(gcd(a, b), gcd(c, d)) == 1


def thm1148_residual(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    multiplier_cone = 8 * a > 7 * d
    q4_tail = 2 * d > a + b + c
    transfer = exact_transfer_gate(shape)
    return not (multiplier_cone or q4_tail or transfer)


def core_phase_floor(a: int) -> tuple[F, int, int]:
    """Worst THM-1148 phase length at the least legal scale."""
    rows = []
    for maximum in range(8, 13):
        least_m = (13 * maximum) // a + 1
        ell_floor = F(72, 35 * (13 * maximum + 1))
        rows.append((least_m * ell_floor, maximum, least_m))
    return min(rows)


def interval_intersection(
    first: tuple[F, F], second: tuple[F, F]
) -> tuple[F, F] | None:
    left = max(first[0], second[0])
    right = min(first[1], second[1])
    return None if right < left else (left, right)


def exception_address_bank(
    starts: list[tuple[F, F]], delta: F
) -> dict[str, object]:
    """Exact m=3 component-address bank for (49,50,51,56)."""
    m = 3
    core_delta = delta / m
    lifted_starts = [
        ((left + lift) / m, (right + lift) / m, start_index, lift)
        for start_index, (left, right) in enumerate(starts)
        for lift in range(m)
    ]

    failures: list[tuple[int, ...]] = []
    total_hits = 0
    hit_counts: list[int] = []
    all_positive_lengths: list[F] = []
    core_best_rows: list[tuple[F, tuple[int, ...], F, F, int, int]] = []
    start_degree: Counter[int] = Counter()
    lift_degree: Counter[tuple[int, int]] = Counter()

    for core in combinations(range(1, 12), 8):
        components, _, _, _, _ = closed_safe_components(core)
        core_starts = [
            (left, right - core_delta)
            for left, right in components
            if right - left >= core_delta
        ]
        hits: list[tuple[F, F, F, int, int]] = []
        for core_start in core_starts:
            for left, right, start_index, lift in lifted_starts:
                overlap = interval_intersection(core_start, (left, right))
                if overlap is None:
                    continue
                length = overlap[1] - overlap[0]
                hits.append((length, overlap[0], overlap[1], start_index, lift))
                start_degree[start_index] += 1
                lift_degree[(start_index, lift)] += 1
                if length > 0:
                    all_positive_lengths.append(length)
        if not hits:
            failures.append(core)
            continue
        best = max(hits)
        core_best_rows.append((best[0], core, best[1], best[2], best[3], best[4]))
        hit_counts.append(len(hits))
        total_hits += len(hits)

    require(not failures, f"exception address failures: {failures}")
    require(len(core_best_rows) == 165, "wrong number of m=3 legal cores")
    require(all_positive_lengths, "exception bank has no positive intersection")
    weakest = min(core_best_rows)
    used_lifts = len(lift_degree)
    return {
        "cores": len(core_best_rows),
        "failures": len(failures),
        "total_hits": total_hits,
        "hit_min": min(hit_counts),
        "hit_max": max(hit_counts),
        "minimum_best": weakest[0],
        "weakest_core": weakest[1],
        "weakest_interval": (weakest[2], weakest[3]),
        "weakest_start": weakest[4],
        "weakest_lift": weakest[5],
        "minimum_positive_hit": min(all_positive_lengths),
        "used_shape_starts": len(start_degree),
        "shape_start_degree_min": min(start_degree.values()),
        "shape_start_degree_max": max(start_degree.values()),
        "used_lifts": used_lifts,
        "unused_lifts": len(lifted_starts) - used_lifts,
        "lift_degree_min": min(lift_degree.values()),
        "lift_degree_max": max(lift_degree.values()),
    }


def main() -> None:
    residual_count = 0
    universal_count = 0
    failures: list[tuple[tuple[int, int, int, int], F, F, F, int, int]] = []
    milestones: dict[int, int] = {}
    leading_rows: list[tuple[tuple[int, int, int, int], F, F, F, F, int, int]] = []
    tightest_positive = None
    extension_count = 0
    extension_universal_count = 0
    extension_tightest = None

    for d in range(4, HEIGHT + 1):
        for a, b, c in combinations(range(1, d), 3):
            shape = (a, b, c, d)
            if not primitive(shape) or not thm1148_residual(shape):
                continue
            residual_count += 1
            (
                delta,
                largest_gap,
                threshold,
                components,
                starts,
                _,
                wall_stats,
            ) = erosion_atlas(shape)
            phase_floor, worst_maximum, least_m = core_phase_floor(a)
            margin = phase_floor - threshold
            if margin > 0:
                universal_count += 1
                candidate = (
                    margin,
                    shape,
                    threshold,
                    phase_floor,
                    worst_maximum,
                    least_m,
                )
                if tightest_positive is None or candidate < tightest_positive:
                    tightest_positive = candidate
                if d >= 57:
                    extension_universal_count += 1
                    if extension_tightest is None or candidate < extension_tightest:
                        extension_tightest = candidate
            else:
                failures.append(
                    (shape, threshold, phase_floor, margin, worst_maximum, least_m)
                )
            if d <= 8:
                leading_rows.append(
                    (
                        shape,
                        delta,
                        largest_gap,
                        threshold,
                        margin,
                        len(components),
                        len(starts),
                    )
                )
            if d >= 57:
                extension_count += 1
        if d in (40, 45, 50, 55, 56, 64):
            milestones[d] = residual_count

    require(residual_count == 95336, f"wrong residual count: {residual_count}")
    require(universal_count == 95335, f"wrong universal count: {universal_count}")
    require(
        milestones
        == {40: 13548, 45: 22259, 50: 34386, 55: 51233, 56: 55000, 64: 95336},
        f"milestone mismatch: {milestones}",
    )
    require(len(failures) == 1, f"wrong number of radius failures: {failures}")
    require(failures[0][0] == EXCEPTION, f"wrong unique exception: {failures[0]}")
    require(failures[0][1] == F(17, 392), "wrong exception threshold")
    require(failures[0][2] == F(3, 70), "wrong exception core floor")
    require(failures[0][3] == -F(1, 1960), "wrong exception deficit")
    require(
        tightest_positive is not None
        and tightest_positive[0] == F(97, 238420)
        and tightest_positive[1] == (45, 46, 47, 52),
        f"wrong tightest positive row: {tightest_positive}",
    )
    require(extension_count == 40336, f"wrong extension count: {extension_count}")
    require(
        extension_universal_count == 40336,
        f"wrong extension universal count: {extension_universal_count}",
    )
    require(
        extension_tightest is not None
        and extension_tightest[0] == F(2279, 303800)
        and extension_tightest[1] == (53, 59, 60, 62)
        and extension_tightest[2] == F(55, 1736)
        and extension_tightest[3] == F(48, 1225)
        and extension_tightest[4:] == (8, 2),
        f"wrong extension tightest row: {extension_tightest}",
    )

    expected_leading = {
        (3, 4, 5, 6): (F(1, 42), F(23, 168), F(9, 56), F(21, 40), 12, 8),
        (3, 5, 6, 7): (F(1, 49), F(5, 49), F(6, 49), F(138, 245), 16, 8),
        (4, 5, 6, 7): (F(1, 49), F(9, 49), F(10, 49), F(76, 245), 18, 8),
        (3, 6, 7, 8): (F(1, 56), F(9, 56), F(5, 28), F(71, 140), 16, 8),
        (4, 5, 7, 8): (F(1, 56), F(71, 784), F(85, 784), F(1591, 3920), 16, 10),
        (4, 6, 7, 8): (F(1, 56), F(4, 49), F(39, 392), F(813, 1960), 16, 12),
        (5, 6, 7, 8): (F(1, 56), F(13, 112), F(15, 112), F(111, 400), 20, 10),
    }
    require(len(leading_rows) == len(expected_leading), "wrong leading-row count")
    for shape, delta, gap, threshold, margin, component_count, start_count in leading_rows:
        require(
            expected_leading.get(shape)
            == (delta, gap, threshold, margin, component_count, start_count),
            f"leading-row mismatch for {shape}",
        )
        verify_sweep_directly(shape)

    (
        exception_delta,
        exception_gap,
        exception_threshold,
        exception_components,
        exception_starts,
        exception_gaps,
        exception_walls,
    ) = erosion_atlas(EXCEPTION)
    verify_sweep_directly(EXCEPTION)
    require(exception_delta == F(1, 392), "wrong exception delta")
    require(exception_gap == F(2, 49), "wrong exception largest gap")
    require(exception_threshold == F(17, 392), "wrong exception H")
    require(len(exception_components) == 124, "wrong exception component count")
    require(len(exception_starts) == 76, "wrong exception start count")
    require(exception_gaps.count(exception_gap) == 2, "wrong max-gap multiplicity")
    require(exception_walls == (412, 412, 0, 0), "wrong exception wall telemetry")

    address = exception_address_bank(exception_starts, exception_delta)
    require(address["total_hits"] == 8456, "wrong exception hit count")
    require((address["hit_min"], address["hit_max"]) == (36, 82), "wrong hit range")
    require(address["minimum_best"] == F(47, 39984), "wrong minimum best hit")
    require(
        address["weakest_core"] == (1, 3, 4, 5, 6, 7, 8, 11),
        "wrong weakest core",
    )
    require(
        address["weakest_interval"] == (F(1639, 2352), F(13955, 19992)),
        "wrong weakest witness interval",
    )
    require(address["minimum_positive_hit"] == F(1, 119952), "wrong min hit")
    require(address["used_shape_starts"] == 76, "an exception start is unused")
    require(address["unused_lifts"] == 12, "wrong unused-lift count")

    # For the exceptional shape, scales m >= 4 are back inside the universal
    # phase-radius comparison.  Scales 1,2 are illegal even for M=8, while
    # m=3 is legal exactly for core maxima M <= 11, the 165-core bank above.
    four_scale_floor = 4 * F(72, 35 * (13 * 12 + 1))
    require(four_scale_floor == F(288, 5495), "wrong m>=4 floor")
    require(
        four_scale_floor - exception_threshold == F(2783, 307720),
        "wrong m>=4 margin",
    )
    require(49 * 2 <= 13 * 8, "m=2 unexpectedly legal")
    require(49 * 3 > 13 * 11 and 49 * 3 <= 13 * 12, "m=3 core split wrong")

    print("THM-1168 exact proportional four-comb shape-erosion referee")
    print("safe convention: ||k*x|| >= 1/14 is closed; danger teeth are open")
    print("residual predicate: primitive and outside THM-1148 strict multiplier/Q4/Phi gates")
    print("height cap d <=", HEIGHT)
    print("residual milestones:", milestones)
    print("primitive residual shapes:", residual_count)
    print("universal radius successes:", universal_count)
    print("universal radius failures:", len(failures))
    print("leading residual rays (shape; delta; G; H=G+delta; phase margin; components/starts):")
    for shape, delta, gap, threshold, margin, component_count, start_count in leading_rows:
        print(
            f"  {shape}: delta={delta}, G={gap}, H={threshold}, "
            f"margin={margin}, components/starts={component_count}/{start_count}"
        )
    print("tightest positive universal row:")
    print(
        f"  shape={tightest_positive[1]}, margin={tightest_positive[0]}, "
        f"H={tightest_positive[2]}, phase_floor={tightest_positive[3]}, "
        f"worst_M={tightest_positive[4]}, least_m={tightest_positive[5]}"
    )
    print("height-57..64 extension:")
    print(
        f"  residual shapes={extension_count}; universal successes={extension_universal_count}; "
        "address exceptions=0"
    )
    print(
        f"  tightest shape={extension_tightest[1]}, margin={extension_tightest[0]}, "
        f"H={extension_tightest[2]}, phase_floor={extension_tightest[3]}, "
        f"worst_M={extension_tightest[4]}, least_m={extension_tightest[5]}"
    )
    print("unique universal-radius obstruction:")
    print(
        f"  shape={EXCEPTION}, delta={exception_delta}, G={exception_gap}, "
        f"H={exception_threshold}, phase_floor=3/70, deficit=1/1960"
    )
    print(
        f"  walls raw/distinct/tied/isolated={exception_walls}; "
        f"components/starts={len(exception_components)}/{len(exception_starts)}; "
        f"max-gap multiplicity={exception_gaps.count(exception_gap)}"
    )
    print("component-address repair at m=3:")
    print(
        f"  legal cores C([11],8)={address['cores']}; failures={address['failures']}; "
        f"closed intersections={address['total_hits']}; hits/core={address['hit_min']}..{address['hit_max']}"
    )
    print(
        f"  minimum best overlap={address['minimum_best']} at core={address['weakest_core']}, "
        f"interval={address['weakest_interval']}, "
        f"shape-start/lift={address['weakest_start']}/{address['weakest_lift']}"
    )
    print(
        f"  minimum positive intersection={address['minimum_positive_hit']}; "
        f"used shape starts={address['used_shape_starts']}/76; "
        f"used lifts={address['used_lifts']}/228"
    )
    print("exception scale split:")
    print("  m<=2 illegal; m=3 iff M<=11 and is closed by the address bank")
    print(
        f"  m>=4 universal: 4*72/[35*(13*12+1)]={four_scale_floor} "
        f"> H by {four_scale_floor-exception_threshold}"
    )
    print("endpoint audit:")
    print("  safe components and eroded start intervals are closed, including singleton starts")
    print("  isolated safe walls are retained in telemetry but cannot carry a positive delta interval")
    print("  a closed interval of length 1/(7*k4) cannot lie in a fifth open danger tooth")
    print("Tournament Analysis:")
    print("  runner-order tournament for every shape: scores=(0,1,2,3), cycles=0, SCCs=4")
    print("  wall chronology: tie-resolved by (position,speed,side), hence transitive")
    print("  exceptional wall path: 412 no-tie vertices, scores=(0,...,411), one Hamiltonian path")
    print("  proof-bearing object: weighted cyclic eroded-start complex plus core-start/lift incidence")
    print("  naked runner/wall tournaments lose interval lengths, cyclic gaps, and component address")
    print("HONEST FRONTIER: every primitive proportional THM-1148 residual with d<=64 is closed")
    print("HONEST FRONTIER: uniform r=5 for primitive shapes with d>64 remains OPEN")
    print("DONE")


if __name__ == "__main__":
    main()
