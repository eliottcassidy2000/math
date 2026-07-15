#!/usr/bin/env python3
"""Exact replay for THM-807's linear selector and grid-ladder limits.

The proof of the selector theorem is elementary and appears in THM-807.  This
artifact checks, with Fraction arithmetic only, the two sharp method-boundary
rows, their complete unit-grid predicates, their exact central return sets,
the endpoint/cusp selectors, and the decorated grid-obligation tournaments.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd


ALPHA = F(1, 13)
BETA = F(1, 11)
GAMMA = BETA - ALPHA
THRESHOLD = F(11, 13)
CLASSES = frozenset(range(1, 7))
PAIR = (13, 5)

ALL_GRID_ROW = (45, 48, 50, 54, 55, 62, 85, 95, 105, 116)
EVEN_GRID_ROW = (6, 9, 20, 24, 30, 36, 42, 54, 66, 90)

# Filled after the first canonical exact run.
EXPECTED_DIGEST = "39ad9ab6c2e7d932103d798b357b1f29beaa960bb4d4dd5b44a9cec895d04a7b"


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def phi(speeds: tuple[int, ...], time: F) -> F:
    return min(norm(speed * time) for speed in speeds)


def q_value(x: int, y: int, time: F) -> F:
    a = (x + y) // 2
    b = abs(x - y) // 2
    return norm(a * time) + norm(b * time)


def folded_class(value: int) -> int:
    residue = value % 13
    return min(residue, (-residue) % 13)


def parity_twisted_class(value: int) -> int:
    return folded_class(value if value % 2 else value // 2)


def balanced_distance(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, (-residue) % modulus)


def unit_representatives(modulus: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                min(value, modulus - value)
                for value in range(1, modulus)
                if gcd(value, modulus) == 1
            }
        )
    )


def deep_unit_profile(
    speeds: tuple[int, ...], x: int, y: int, multiplier: int
) -> tuple[tuple[int, F, F], ...]:
    modulus = 13 * multiplier
    answer = []
    for numerator in unit_representatives(modulus):
        time = F(numerator, modulus)
        margin = phi(speeds, time)
        if margin >= BETA:
            answer.append((numerator, margin, q_value(x, y, time)))
    return tuple(answer)


def deep_unit_classes(
    speeds: tuple[int, ...], modulus: int
) -> tuple[tuple[int, F], ...]:
    answer = []
    for numerator in unit_representatives(modulus):
        time = F(numerator, modulus)
        margin = phi(speeds, time)
        if margin >= BETA:
            answer.append((numerator, margin))
    return tuple(answer)


def intersect_intervals(
    left: list[tuple[F, F]],
    right: list[tuple[F, F]],
    *,
    allow_points: bool,
) -> list[tuple[F, F]]:
    answer: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi or (allow_points and lo == hi):
            answer.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        elif left[i][1] > right[j][1]:
            j += 1
        else:
            i += 1
            j += 1
    return answer


def deep_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current = [(F(0), F(1))]
    for speed in speeds:
        safe = [
            ((F(k) + BETA) / speed, (F(k + 1) - BETA) / speed)
            for k in range(speed)
        ]
        current = intersect_intervals(current, safe, allow_points=True)
    return tuple(current)


def return_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    """Closures of the strict return components, represented in [-1/2,1/2]."""
    current = [(F(-1, 2), F(1, 2))]
    for speed in speeds:
        allowed = []
        for integer in range(-speed, speed + 1):
            interval = (
                (F(integer) - GAMMA) / speed,
                (F(integer) + GAMMA) / speed,
            )
            if interval[1] > F(-1, 2) and interval[0] < F(1, 2):
                allowed.append(interval)
        current = intersect_intervals(current, allowed, allow_points=False)
    return tuple(current)


def merge_closed(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    intervals.sort()
    answer: list[tuple[F, F]] = []
    for left, right in intervals:
        if answer and left <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(answer[-1][1], right))
        else:
            answer.append((left, right))
    return answer


def split_circle_interval(left: F, right: F) -> list[tuple[F, F]]:
    if right - left >= 1:
        return [(F(0), F(1))]
    while left < 0:
        left += 1
        right += 1
    while left >= 1:
        left -= 1
        right -= 1
    if right <= 1:
        return [(left, right)]
    return [(left, F(1)), (F(0), right - 1)]


def circle_sum_components(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    pieces = []
    for left, right in first:
        for shift_left, shift_right in second:
            pieces.extend(
                split_circle_interval(left + shift_left, right + shift_right)
            )
    return tuple(merge_closed(pieces))


def component_candidates(
    component: tuple[F, F], x: int, y: int
) -> tuple[F, ...]:
    left, right = component
    a = (x + y) // 2
    b = abs(x - y) // 2
    candidates = {left, right}
    for frequency in (a, b):
        for integer in range(-1, 2 * frequency + 2):
            point = F(integer, 2 * frequency) % 1
            if left <= point <= right:
                candidates.add(point)
    return tuple(sorted(candidates))


def selector_profile(
    components: tuple[tuple[F, F], ...], x: int, y: int
) -> tuple[tuple[F, F, F, int], ...]:
    """(escape margin, argmin, width, selector size) for every component."""
    answer = []
    for component in components:
        candidates = component_candidates(component, x, y)
        minimum, argmin = min((q_value(x, y, t), t) for t in candidates)
        answer.append(
            (
                THRESHOLD - minimum,
                argmin,
                component[1] - component[0],
                len(candidates),
            )
        )
    return tuple(answer)


def breakpoint_candidates(speeds: tuple[int, ...]) -> set[F]:
    denominators = {2 * speed for speed in speeds}
    for left, right in combinations(speeds, 2):
        denominators.add(left + right)
        denominators.add(abs(left - right))
    denominators.discard(0)
    return {
        F(numerator, denominator)
        for denominator in denominators
        for numerator in range(denominator + 1)
    }


def exact_loneliness(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates = breakpoint_candidates(speeds)
    value = max(phi(speeds, t) for t in candidates)
    return value, tuple(sorted(t for t in candidates if phi(speeds, t) == value))


def endpoint_owner_audit(
    speeds: tuple[int, ...],
    deep: tuple[tuple[F, F], ...],
    thickened: tuple[tuple[F, F], ...],
    radius: F,
) -> None:
    deep_endpoints = {endpoint for component in deep for endpoint in component}
    legal_deep_endpoints = {
        (F(11 * integer + sign, 11 * speed) % 1)
        for speed in speeds
        for integer in range(speed + 1)
        for sign in (-1, 1)
    }
    assert deep_endpoints <= legal_deep_endpoints
    shifted = {
        ((endpoint + sign * radius) % 1)
        for endpoint in deep_endpoints
        for sign in (-1, 1)
    }
    # Ignore the artificial 0/1 cut if a circular component wraps there.
    for left, right in thickened:
        if left not in (F(0), F(1)):
            assert left in shifted
        if right not in (F(0), F(1)):
            assert right in shifted


def total_order_edges(order: tuple[int, ...]) -> set[tuple[int, int]]:
    rank = {vertex: index for index, vertex in enumerate(order)}
    return {
        (left, right) if rank[left] < rank[right] else (right, left)
        for left, right in combinations(order, 2)
    }


def tournament_fingerprint(
    vertices: tuple[int, ...], edges: set[tuple[int, int]]
) -> tuple[tuple[int, ...], int, tuple[int, ...], int]:
    scores = {vertex: 0 for vertex in vertices}
    out = {vertex: set() for vertex in vertices}
    for source, target in edges:
        scores[source] += 1
        out[source].add(target)

    cycles = 0
    for triple in combinations(vertices, 3):
        restricted = sorted(
            sum(target in triple for target in out[vertex]) for vertex in triple
        )
        cycles += restricted == [1, 1, 1]

    # Kosaraju is overkill at n<=10; mutual reachability is exact and clearer.
    reach = {vertex: {vertex} | set(out[vertex]) for vertex in vertices}
    changed = True
    while changed:
        changed = False
        for vertex in vertices:
            enlarged = set().union(*(reach[target] for target in tuple(reach[vertex])))
            if not enlarged <= reach[vertex]:
                reach[vertex] |= enlarged
                changed = True
    unseen = set(vertices)
    scc_sizes = []
    while unseen:
        vertex = min(unseen)
        component = {
            target
            for target in unseen
            if target in reach[vertex] and vertex in reach[target]
        }
        unseen -= component
        scc_sizes.append(len(component))

    # Subset DP counts directed Hamiltonian paths exactly.
    index = {vertex: i for i, vertex in enumerate(vertices)}
    count = [[0] * len(vertices) for _ in range(1 << len(vertices))]
    for i in range(len(vertices)):
        count[1 << i][i] = 1
    for mask in range(1 << len(vertices)):
        for end in range(len(vertices)):
            if not count[mask][end]:
                continue
            source = vertices[end]
            for target in out[source]:
                j = index[target]
                if not mask & (1 << j):
                    count[mask | (1 << j)][j] += count[mask][end]
    paths = sum(count[-1])
    return (
        tuple(sorted(scores.values())),
        cycles,
        tuple(sorted(scc_sizes, reverse=True)),
        paths,
    )


def grid_tournament(
    grid: dict[int, tuple[tuple[int, F, F], ...]]
) -> dict[str, object]:
    vertices = tuple(sorted(grid))

    def detection(multiplier: int) -> F:
        rows = grid[multiplier]
        return max((THRESHOLD - row[2] for row in rows), default=F(-1))

    margin_order = tuple(
        sorted(vertices, key=lambda d: (detection(d), d))
    )
    population_order = tuple(
        sorted(
            vertices,
            key=lambda d: (
                len(unit_representatives(13 * d)) - len(grid[d]),
                d,
            ),
        )
    )
    margin_edges = total_order_edges(margin_order)
    population_edges = total_order_edges(population_order)
    return {
        "margin_order": margin_order,
        "population_order": population_order,
        "flips": len(margin_edges ^ population_edges) // 2,
        "margin_fingerprint": tournament_fingerprint(vertices, margin_edges),
        "population_fingerprint": tournament_fingerprint(vertices, population_edges),
    }


def audit_row(
    speeds: tuple[int, ...], multipliers: tuple[int, ...]
) -> dict[str, object]:
    x, y = PAIR
    assert len(speeds) == len(set(speeds)) == 10 and gcd(*speeds) == 1
    assert all(speed % 13 for speed in speeds)
    assert all(
        any(speed % modulus == 0 for speed in speeds)
        for modulus in range(2, 13)
    )
    assert {speed % 13 for speed in speeds} == set(range(1, 13)) - {5, 8}
    assert {folded_class(speed) for speed in speeds} == {1, 2, 3, 4, 6}
    assert {parity_twisted_class(speed) for speed in speeds} == CLASSES

    value, maximizers = exact_loneliness(speeds)
    B = max(speeds)
    rho = (value - ALPHA) / B
    astar_left = F(1, x * y) + 2 * rho
    astar_right = F(2, 13 * x) + F(2, 13 * y)
    assert x <= 2 * B - 1 and y <= B - 1
    assert 13 * B + 2 * x * y <= 2 * B * (x + y)
    assert astar_left <= astar_right

    assert deep_unit_classes(speeds, 5) == ()
    assert deep_unit_classes(speeds, 13) == ((5, F(2, 13)),)
    grid = {d: deep_unit_profile(speeds, x, y, d) for d in multipliers}

    deep = deep_components(speeds)
    returns = return_components(speeds)
    radius = GAMMA / B
    assert returns == ((-radius, radius),)
    thickened = circle_sum_components(deep, returns)
    assert len(thickened) <= len(deep) <= sum(speeds)
    endpoint_owner_audit(speeds, deep, thickened, radius)

    profile = selector_profile(thickened, x, y)
    selector_size = sum(row[3] for row in profile)
    assert selector_size <= 2 * sum(speeds) + 2 * max(x, y)
    selector_minimum = min((THRESHOLD - row[0], row[1]) for row in profile)
    tournament = grid_tournament(grid)

    return {
        "speeds": speeds,
        "B": B,
        "value": value,
        "maximizers": maximizers,
        "rho": rho,
        "astar_left": astar_left,
        "astar_right": astar_right,
        "grid": grid,
        "deep_count": len(deep),
        "returns": returns,
        "thickened_count": len(thickened),
        "selector_size": selector_size,
        "selector_minimum": selector_minimum,
        "tournament": tournament,
    }


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonicalize(value: object) -> str:
    if isinstance(value, F):
        return fmt(value)
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{canonicalize(key)}:{canonicalize(value[key])}"
            for key in sorted(value, key=str)
        ) + "}"
    if isinstance(value, (tuple, list)):
        return "(" + ",".join(canonicalize(item) for item in value) + ")"
    return str(value)


def main() -> None:
    all_row = audit_row(ALL_GRID_ROW, tuple(range(1, 9)))
    even_row = audit_row(EVEN_GRID_ROW, tuple(range(2, 21, 2)))

    assert all(
        all(q >= THRESHOLD for _p, _margin, q in all_row["grid"][d])
        for d in range(1, 8)
    )
    assert all_row["grid"][8] == (
        (31, F(5, 52), F(53, 104)),
        (45, F(11, 104), F(3, 8)),
    )

    assert all(
        all(q >= THRESHOLD for _p, _margin, q in even_row["grid"][d])
        for d in range(2, 20, 2)
    )
    assert even_row["grid"][20] == (
        (73, F(7, 65), F(31, 52)),
    )

    assert all_row["value"] == F(45, 161)
    assert all_row["maximizers"] == (F(1, 161), F(160, 161))
    assert all_row["rho"] == F(106, 60697)
    assert all_row["deep_count"] == all_row["thickened_count"] == 106
    assert all_row["selector_size"] == 212
    assert all_row["selector_minimum"] == (F(709, 28710), F(709, 373230))

    assert even_row["value"] == F(2, 13)
    assert even_row["maximizers"] == (
        F(11, 39), F(5, 13), F(8, 13), F(28, 39)
    )
    assert even_row["rho"] == F(1, 1170)
    assert even_row["deep_count"] == even_row["thickened_count"] == 48
    assert even_row["selector_size"] == 96
    assert even_row["selector_minimum"] == (
        F(1159, 5445), F(1159, 70785)
    )

    canonical = canonicalize((all_row, even_row))
    digest = sha256(canonical.encode()).hexdigest()

    print("THM-807 connected-return linear selector replay")
    print("all verdicts use integer/Fraction arithmetic")
    print()
    for label, row in (("all-grid", all_row), ("even-grid", even_row)):
        print(f"{label}_row={row['speeds']}")
        print(
            f"  B={row['B']} M={fmt(row['value'])} "
            f"maximizers={tuple(fmt(t) for t in row['maximizers'])}"
        )
        print(
            f"  rho={fmt(row['rho'])} Astar={fmt(row['astar_left'])}"
            f"<={fmt(row['astar_right'])}"
        )
        print(f"  unit_grid_profiles={row['grid']}")
        print(
            f"  return={row['returns']} E_components={row['deep_count']} "
            f"K_components={row['thickened_count']} selector={row['selector_size']}"
        )
        print(
            "  selector_min_(Q,t)="
            f"({fmt(row['selector_minimum'][0])},{fmt(row['selector_minimum'][1])})"
        )
        tournament = row["tournament"]
        print("  Tournament Analysis (vertices = multiplier-grid obligations)")
        print("    observable: signed grid-detection margin; switch: danger-covered unit population")
        print("    tie Hamiltonian path: increasing multiplier")
        print(f"    margin_order={tournament['margin_order']}")
        print(f"    population_order={tournament['population_order']}")
        print(f"    edge_flips={tournament['flips']}")
        print(f"    margin_fingerprint={tournament['margin_fingerprint']}")
        print(f"    population_fingerprint={tournament['population_fingerprint']}")
        print()
    print(f"sha256={digest}")
    print("challenged vertices: runners and denominator grids")
    print("preserved by grid tournament: relative finite-grid detection difficulty only")
    print("destroyed: non-grid deep components, endpoint owners, return incidence, metric sign")
    print("predicate carrier: owner-labelled components plus exact endpoint/cusp selector")
    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert digest == EXPECTED_DIGEST
        print("FINAL: PASS")
    else:
        print("FINAL: FIRST EXACT RUN — install digest before canonizing")


if __name__ == "__main__":
    main()
