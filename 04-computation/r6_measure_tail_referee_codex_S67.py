#!/usr/bin/env python3
"""Exact referee for the unbounded r=6 measure-horn tail.

THM-1102 scanned five removed killers only in the bottom window
``[13*max(P)+1, 13*max(P)+17)`` and found ``max T = 308.4`` there.  An
interior maximizer proves the maximum of that *finite window*, but it does
not control translated, scaled, or otherwise separated quintuples.

This dependency-free verifier does four things, all with ``Fraction``:

1. reproduces THM-1102's displayed window maximizer;
2. gives an exact translated quintuple with ``T > 333`` (and an explicit
   lonely time, so this refutes only the proof extrapolation, not LRC);
3. verifies representative translated/scaled/spaced/highly-composite
   probes and records a general lower bound showing that no constant
   upper bound on ``T`` can hold on the unbounded tail;
4. exactly integrates the rational two-dimensional limit atlas of an
   infinite Covering progression for which ``T/k5 -> 28/27 > 1`` even
   though every member is lonely at ``t=5/16``.

The threshold is the one used in THM-1102:

    T(E) = min(N(E)/(6*measure(E)), 1/(3*largest_component(E))).

Tournament Analysis is diagnostic only.  On the five killer vertices we
orient by the one-killer threshold and break ties by the integer label.
That gauge is transitive (score histogram 0..4, no cycles, singleton SCCs,
one Hamiltonian path), but it destroys the residual-component incidence
and all higher intersections.  The proof-facing object is the exact
interval complex after all five removals, not this tournament quotient.
"""

from fractions import Fraction as F
from itertools import combinations


CORE = (1, 2, 4, 7, 9, 11, 12)
WINDOW_WITNESS = (158, 160, 162, 164, 166)
TAIL_WITNESS = (290, 292, 294, 296, 298)
COVERING_SIXTH = 338
NONCOVERING_SIXTH = 333
OFFSET_PATTERN = (0, 10, 12, 14, 18)


def require(condition: bool, message: str) -> None:
    """Proof-bearing check that remains active under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def norm_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def safe_set(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Exact positive-length components of {t in [0,1]: ||v t|| >= 1/14}."""
    breakpoints = {F(0), F(1)}
    for speed in speeds:
        for j in range(speed + 1):
            radius = F(1, 14 * speed)
            for endpoint in (F(j, speed) - radius, F(j, speed) + radius):
                if 0 <= endpoint <= 1:
                    breakpoints.add(endpoint)

    ordered = sorted(breakpoints)
    result: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(norm_distance(speed * midpoint) >= F(1, 14) for speed in speeds):
            if result and result[-1][1] == left:
                result[-1] = (result[-1][0], right)
            else:
                result.append((left, right))
    return result


def remove_bad(intervals: list[tuple[F, F]], speed: int) -> list[tuple[F, F]]:
    """Subtract {t: ||speed*t|| < 1/14}, retaining exact rational endpoints."""
    result: list[tuple[F, F]] = []
    for left, right in intervals:
        first = (left * speed).__floor__() - 1
        last = (right * speed).__floor__() + 1
        cursor = left
        for j in range(first, last + 1):
            bad_left = F(14 * j - 1, 14 * speed)
            bad_right = F(14 * j + 1, 14 * speed)
            if bad_right <= left or bad_left >= right:
                continue
            bad_left = max(left, bad_left)
            bad_right = min(right, bad_right)
            if bad_left > cursor:
                result.append((cursor, bad_left))
            cursor = max(cursor, bad_right)
        if cursor < right:
            result.append((cursor, right))
    return result


def residual(core: tuple[int, ...], killers: tuple[int, ...]) -> list[tuple[F, F]]:
    intervals = safe_set(core)
    for killer in killers:
        intervals = remove_bad(intervals, killer)
    return intervals


def metrics(core: tuple[int, ...], killers: tuple[int, ...]) -> dict[str, F | int | str]:
    intervals = residual(core, killers)
    require(bool(intervals), f"empty residual for core={core}, killers={killers}")
    measure = sum((right - left for left, right in intervals), F(0))
    components = len(intervals)
    largest = max(right - left for left, right in intervals)
    count_threshold = F(components, 6) / measure
    component_threshold = F(1, 3) / largest
    threshold = min(count_threshold, component_threshold)
    return {
        "components": components,
        "measure": measure,
        "largest": largest,
        "count_threshold": count_threshold,
        "component_threshold": component_threshold,
        "threshold": threshold,
        "branch": "count" if count_threshold <= component_threshold else "component",
        "ratio": threshold / killers[-1],
    }


def least_abs_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def rational_loneliness(speeds: tuple[int, ...], q: int, a: int) -> F:
    return min(F(least_abs_residue(a * speed, q), q) for speed in speeds)


def individual_thresholds(core: tuple[int, ...], killers: tuple[int, ...]) -> dict[int, F]:
    return {killer: metrics(core, (killer,))["threshold"] for killer in killers}


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(speed % q == 0 for speed in speeds) for q in range(2, 15))


def fiber_atlas_metrics(
    core: tuple[int, ...], offsets: tuple[int, ...]
) -> tuple[F, F, F]:
    """Exact area, component intensity, and max horizontal gap of the limit atlas.

    Put x={m t}.  The translated killers m+a are safe precisely when
    ||x+a*t||>=1/14.  Their two-dimensional limit region over the core-safe
    t-set is cut into finitely many rational polygons by the ten bad-arc
    endpoints x=-a*t +/- 1/14.  This routine integrates its vertical-fiber
    measure and component count and maximizes its fiber-component length.
    """
    lam = F(1, 14)
    expressions = tuple((offset, sign) for offset in offsets for sign in (-1, 1))
    breakpoints = {F(0), F(1)}

    # Endpoint wraps: -a*t+s/14 is integral.
    for offset, sign in expressions:
        if offset == 0:
            continue
        for integer in range(-offset - 2, 3):
            t = F(sign, 14 * offset) - F(integer, offset)
            if 0 <= t <= 1:
                breakpoints.add(t)

    # Pair-order changes modulo one.
    for (offset_a, sign_a), (offset_b, sign_b) in combinations(expressions, 2):
        difference = offset_b - offset_a
        if difference == 0:
            continue
        constant = F(sign_a - sign_b, 14)
        for integer in range(-abs(difference) - 2, abs(difference) + 3):
            t = (F(integer) - constant) / difference
            if 0 <= t <= 1:
                breakpoints.add(t)

    core_safe = safe_set(core)
    for left, right in core_safe:
        breakpoints.add(left)
        breakpoints.add(right)

    def core_is_safe(t: F) -> bool:
        return any(left <= t <= right for left, right in core_safe)

    def fiber_is_safe(t: F, x: F) -> bool:
        return all(norm_distance(x + offset * t) >= lam for offset in offsets)

    area = F(0)
    component_intensity = F(0)
    max_gap = F(0)
    ordered_breakpoints = sorted(breakpoints)

    for left, right in zip(ordered_breakpoints, ordered_breakpoints[1:]):
        if left == right:
            continue
        midpoint = (left + right) / 2
        if not core_is_safe(midpoint):
            continue

        # Freeze the integer sheet of each affine endpoint on this t-cell.
        endpoints = []
        for offset, sign in expressions:
            raw_midpoint = -offset * midpoint + F(sign, 14)
            sheet = raw_midpoint.__floor__()
            residue = raw_midpoint - sheet
            endpoints.append((residue, -offset, F(sign, 14) - sheet))
        endpoints.sort()
        require(
            all(endpoints[i][0] != endpoints[i + 1][0] for i in range(9)),
            f"unrecorded endpoint collision on t-cell [{left},{right}]",
        )

        safe_length_left = F(0)
        safe_length_right = F(0)
        safe_components = 0
        for index, endpoint in enumerate(endpoints):
            successor = endpoints[(index + 1) % len(endpoints)]
            wrap = 1 if index + 1 == len(endpoints) else 0
            midpoint_gap = successor[0] + wrap - endpoint[0]
            test_x = (endpoint[0] + midpoint_gap / 2) % 1
            if not fiber_is_safe(midpoint, test_x):
                continue

            slope = successor[1] - endpoint[1]
            intercept = successor[2] + wrap - endpoint[2]
            gap_left = slope * left + intercept
            gap_right = slope * right + intercept
            require(
                gap_left >= 0 and gap_right >= 0,
                f"negative affine fiber gap on t-cell [{left},{right}]",
            )
            safe_length_left += gap_left
            safe_length_right += gap_right
            safe_components += 1
            max_gap = max(max_gap, gap_left, gap_right)

        width = right - left
        area += width * (safe_length_left + safe_length_right) / 2
        component_intensity += width * safe_components

    return area, component_intensity, max_gap


def print_metric(label: str, killers: tuple[int, ...]) -> None:
    data = metrics(CORE, killers)
    print(label)
    print("  killers: " + " ".join(map(str, killers)))
    print(f"  components N: {data['components']}")
    print(f"  measure mu: {data['measure']}")
    print(f"  largest component L: {data['largest']}")
    print(f"  N/(6mu): {data['count_threshold']}")
    print(f"  1/(3L): {data['component_threshold']}")
    print(f"  T ({data['branch']} branch): {data['threshold']}")
    print(f"  R=T/k5: {data['ratio']}")


def main() -> None:
    print("### exact r=6 measure-tail referee ###")
    print("core: " + " ".join(map(str, CORE)))
    print()

    print_metric("THM-1102 finite-window maximizer", WINDOW_WITNESS)
    window = metrics(CORE, WINDOW_WITNESS)
    require(window["components"] == 206, "window component count changed")
    require(window["largest"] == F(67, 61992), "window largest component changed")
    require(window["threshold"] == F(20664, 67), "window threshold changed")
    print()

    print_metric("translated counterexample to max-T <= 333", TAIL_WITNESS)
    tail = metrics(CORE, TAIL_WITNESS)
    require(tail["components"] == 370, "tail component count changed")
    require(
        tail["measure"] == F(214450238441, 1981564300485),
        "tail residual measure changed",
    )
    require(tail["largest"] == F(1, 1043), "tail largest component changed")
    require(tail["threshold"] == F(1043, 3), "tail threshold changed")
    require(tail["threshold"] > 333, "tail no longer exceeds 333")
    require(tail["threshold"] > window["threshold"], "tail no longer beats window")
    require(tail["threshold"] > COVERING_SIXTH, "covering sixth no longer defeats horn")
    print(f"  exact comparison: T-333 = {tail['threshold'] - 333}")
    print(f"  exact comparison: T-338 = {tail['threshold'] - 338}")
    print(f"  exact comparison: T-window_max = {tail['threshold'] - window['threshold']}")
    print()

    noncovering_family = CORE + TAIL_WITNESS + (NONCOVERING_SIXTH,)
    noncovering_lonely = rational_loneliness(noncovering_family, 13, 1)
    require(not is_covering(noncovering_family), "333 variant unexpectedly covering")
    require(
        noncovering_lonely == F(1, 13) > F(1, 14),
        "333 variant lost its q=13 witness",
    )

    covering_family = CORE + TAIL_WITNESS + (COVERING_SIXTH,)
    covering_lonely = rational_loneliness(covering_family, 16, 5)
    require(is_covering(covering_family), "338 countercertificate is not Covering")
    require(
        covering_lonely == F(1, 8) > F(1, 14),
        "338 countercertificate lost its q=16 witness",
    )
    print("scope guardrail: the countercertificate survives Covering")
    print("  sixth killer: 338 (outside THM-1121's certified interval 92..332)")
    print(f"  the THM-1102 measure sufficient condition fails: 338 < T = {tail['threshold']}")
    print("  divisibility carriers q=2..14: " + " ".join(
        str(next(speed for speed in covering_family if speed % q == 0))
        for q in range(2, 15)
    ))
    print(f"  nevertheless t=5/16 has exact minimum distance {covering_lonely}")
    print("  the k6=333 variant is noncovering and is dispatched at t=1/13")
    print("  conclusion: covering proof-route counterexample, NOT an LRC counterexample")
    print()

    print("### deterministic tail probes (all exact) ###")
    probes = (
        ("translated-step-2", (1005, 1007, 1009, 1011, 1013)),
        ("scaled", (855, 865, 875, 885, 895)),
        ("spaced", (400, 800, 1200, 1600, 2000)),
        ("highly-composite", (840, 900, 960, 1020, 1080)),
        ("near-consecutive", (1000, 1001, 1003, 1007, 1009)),
        ("geometric", (200, 400, 800, 1600, 3200)),
    )
    for label, killers in probes:
        data = metrics(CORE, killers)
        print(
            f"{label:18s} k5={killers[-1]:4d} "
            f"T={data['threshold']} R={data['ratio']} "
            f"N={data['components']} L={data['largest']}"
        )
    print()

    # Infinite obstruction to a *constant-T* tail statement.  For m == 4 (mod 13),
    # every speed m,m+2,...,m+8 is nonzero mod 13, so t=1/13 keeps the residual
    # nonempty.  Since the residual is contained in Safe(m), L <= 6/(7m).
    # Also mu <= N*L, hence both terms defining T are >= 1/(6L).
    # Therefore T >= 7m/36 -> infinity.
    print("### infinite translated obstruction to any constant upper bound on T ###")
    for m in (290, 1005, 1733):
        require(m % 13 == 4, f"translated probe m={m} has wrong residue")
        killers = tuple(m + 2 * j for j in range(5))
        data = metrics(CORE, killers)
        lower_bound = F(7 * m, 36)
        require(data["threshold"] >= lower_bound, f"lower bound failed at m={m}")
        require(
            rational_loneliness(CORE + killers, 13, 1) == F(1, 13),
            f"q=13 witness failed at m={m}",
        )
        print(
            f"m={m:4d} T={data['threshold']} >= 7m/36={lower_bound}; "
            f"R={data['ratio']}"
        )
    require(F(7 * 1733, 36) > 333, "chosen lower-bound probe does not exceed 333")
    print("symbolic law: m == 4 (mod 13) => E_m nonempty and T(E_m) >= 7m/36")
    print("therefore sup_m T(E_m) = infinity; a valid tail lemma must compare T to k6")
    print("this step-2 family alone does not decide uniform R; the next family refutes it")
    print()

    print("### infinite covering obstruction to the measure horn itself ###")
    area, component_intensity, max_gap = fiber_atlas_metrics(CORE, OFFSET_PATTERN)
    count_limit = component_intensity / (6 * area)
    component_limit = F(1, 3) / max_gap
    threshold_limit = min(count_limit, component_limit)
    require(area == F(76004881, 627525360), "limit-atlas area changed")
    require(
        component_intensity == F(14813, 13860),
        "limit-atlas component intensity changed",
    )
    require(max_gap == F(9, 28), "limit-atlas maximum fiber gap changed")
    require(count_limit == F(111778898, 76004881), "count-branch limit changed")
    require(
        component_limit == threshold_limit == F(28, 27) > 1,
        "component-branch limit is no longer 28/27",
    )
    print("offsets: " + " ".join(map(str, OFFSET_PATTERN)))
    print(f"limit-atlas area: {area}")
    print(f"limit-atlas component intensity: {component_intensity}")
    print(f"limit-atlas maximum fiber gap: {max_gap}")
    print(f"lim count branch T/m: {count_limit}")
    print(f"lim component branch T/m: {component_limit}")
    print(f"lim T/(m+18): {threshold_limit} > 1")

    for m in (1360, 5000, 8640):
        killers = tuple(m + offset for offset in OFFSET_PATTERN)
        data = metrics(CORE, killers)
        sixth = m + 19
        full = CORE + killers + (sixth,)
        require(is_covering(full), f"finite obstruction m={m} is not Covering")
        require(data["threshold"] > sixth, f"measure horn succeeds at probe m={m}")
        print(
            f"m={m:4d} k6={sixth} T={data['threshold']} "
            f"T-k6={data['threshold'] - sixth}"
        )

    covering_base = 1360
    covering_period = 3640
    covering_base_family = CORE + tuple(
        covering_base + offset for offset in OFFSET_PATTERN
    ) + (covering_base + 19,)
    require(is_covering(covering_base_family), "covering progression base is not Covering")
    require(
        all(covering_period % q == 0 for q in (5, 8, 10, 13, 14)),
        "covering period does not preserve all required carriers",
    )
    print("covering progression: m=1360+3640n")
    print("  every family is covering (carriers are periodic modulo 5,8,10,13,14)")

    progression_base = 5000
    progression_period = 7280
    base_family = CORE + tuple(progression_base + offset for offset in OFFSET_PATTERN) + (
        progression_base + 19,
    )
    require(is_covering(base_family), "lonely subprogression base is not Covering")
    require(
        rational_loneliness(base_family, 16, 5) == F(1, 8),
        "lonely subprogression lost its q=16 witness",
    )
    require(progression_period % 16 == 0, "lonely period does not preserve q=16 residues")
    require(
        all(progression_period % q == 0 for q in (5, 8, 10, 13, 14)),
        "lonely period does not preserve Covering carriers",
    )
    print("lonely subprogression: m=5000+7280n")
    print("  every family is explicitly lonely at t=5/16 with minimum 1/8")
    print("  equidistribution/crossing limit 28/27 implies T>k6 for all sufficiently large n")
    print("  conclusion: infinitely many covering failures of this measure certificate")
    print()

    print("### Tournament Analysis (diagnostic quotient) ###")
    loads = individual_thresholds(CORE, TAIL_WITNESS)
    order = sorted(TAIL_WITNESS, key=lambda killer: (loads[killer], killer))
    print("vertices: the five removed killers")
    print("observable: one-killer threshold difference; tie gauge: integer label")
    print("Hamiltonian path: " + " -> ".join(map(str, order)))
    print("score histogram: 0 1 2 3 4")
    print("directed cycles: 0; SCCs: 5 singleton components; Hamiltonian paths: 1")
    print("destroyed by quotient: residual components and all >=2-way bad-set intersections")
    print("proof-facing object: exact residual interval complex")
    print("DONE")


if __name__ == "__main__":
    main()
