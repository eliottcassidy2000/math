#!/usr/bin/env python3
"""Exact referee for THM-2091's centered triple-energy obstruction.

The analytic inputs are the THM-1234 seven-comb pair floor and THM-2080's
mixed-overlap formula.  Everything new here is audited with exact rational
arithmetic: Cayley-tree incidence, the odd reduced-charge ledger, the
centered third moment, and the complete rank-seven terminal shell through
speed 24.  Runtime checks remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction as F
from heapq import heapify, heappop, heappush
from itertools import combinations, product
from math import gcd


N_VERTICES = 7
MAX_CORE_SPEED = 24
DANGER_RADIUS = F(1, 14)
GUARD_RADIUS = F(1, 7)
PAIR_FLOOR = F(22, 65)
CHARGE_FLOOR = -F(41, 495)
BASE_MARGIN = F(2059, 315315)
ENERGY_THRESHOLD = F(2059, 90090)
MODULAR_THRESHOLDS = (
    F(2059, 90090),
    F(396587, 8828820),
    F(608227, 8828820),
    F(273289, 2942940),
    F(368527, 2942940),
)

FULL_DIVISOR_MASK = (1 << 13) - 1
DIVISOR_MASK = {
    q: sum(1 << (d - 2) for d in range(2, 15) if q % d == 0)
    for q in range(1, MAX_CORE_SPEED + 1)
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def comb_boundaries(speed: int, radius: F) -> set[F]:
    points: set[F] = set()
    for k in range(speed + 1):
        center = F(k, speed)
        for sign in (-1, 1):
            point = center + sign * radius / speed
            if 0 <= point <= 1:
                points.add(point)
    return points


def in_comb(speed: int, radius: F, t: F) -> bool:
    return circle_distance(speed * t) < radius


def atomized_two_comb(
    first_speed: int,
    first_radius: F,
    second_speed: int,
    second_radius: F,
) -> F:
    boundaries = {F(0), F(1)}
    boundaries.update(comb_boundaries(first_speed, first_radius))
    boundaries.update(comb_boundaries(second_speed, second_radius))
    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if in_comb(first_speed, first_radius, midpoint) and in_comb(
            second_speed, second_radius, midpoint
        ):
            measure += right - left
    return measure


def atomized_centered_triple(h: int, p: int, q: int) -> tuple[F, F]:
    """Return outside-guard overlap and centered triple correlation."""
    boundaries = {F(0), F(1)}
    for speed, radius in (
        (h, GUARD_RADIUS),
        (p, DANGER_RADIUS),
        (q, DANGER_RADIUS),
    ):
        boundaries.update(comb_boundaries(speed, radius))
    points = sorted(boundaries)
    outside_overlap = F(0)
    centered = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        guard = in_comb(h, GUARD_RADIUS, midpoint)
        first = in_comb(p, DANGER_RADIUS, midpoint)
        second = in_comb(q, DANGER_RADIUS, midpoint)
        length = right - left
        if first and second and not guard:
            outside_overlap += length
        centered += (
            length
            * (F(int(first)) - F(1, 7))
            * (F(int(second)) - F(1, 7))
            * (F(int(guard)) - F(2, 7))
        )
    return outside_overlap, centered


def atomized_centered_moment(Q: tuple[int, ...], h: int) -> F:
    boundaries = {F(0), F(1)}
    boundaries.update(comb_boundaries(h, GUARD_RADIUS))
    for q in Q:
        boundaries.update(comb_boundaries(q, DANGER_RADIUS))
    points = sorted(boundaries)
    total = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        guard = int(in_comb(h, GUARD_RADIUS, midpoint))
        multiplicity = sum(in_comb(q, DANGER_RADIUS, midpoint) for q in Q)
        polynomial = (
            multiplicity * multiplicity
            - F(19, 7) * multiplicity
            + F(6, 7)
        )
        total += (
            (right - left)
            * (F(guard) - F(2, 7))
            * polynomial
            / 2
        )
    return total


def mixed_overlap_formula(h: int, q: int) -> F:
    g = gcd(h, q)
    a, b = h // g, q // g
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, 49) + F(2, a * b) * correction


def global_pair_formula(p: int, q: int) -> F:
    g = gcd(p, q)
    a, b = p // g, q // g

    def fold(r: int) -> int:
        residue = r % 14
        return residue * (14 - residue)

    return F(1, 49) + F(fold(a + b) - fold(b - a), 196 * a * b)


def reduced_charge(a: int, b: int) -> F:
    require(a > 0 and b > 0 and gcd(a, b) == 1, "charge pair not reduced")
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, a * b) * correction


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def prufer_edges(word: tuple[int, ...], n: int = N_VERTICES) -> tuple[tuple[int, int], ...]:
    degrees = [1] * n
    for vertex in word:
        degrees[vertex] += 1
    leaves = [vertex for vertex, degree in enumerate(degrees) if degree == 1]
    heapify(leaves)
    edges: list[tuple[int, int]] = []
    for vertex in word:
        leaf = heappop(leaves)
        edges.append(tuple(sorted((leaf, vertex))))
        degrees[leaf] -= 1
        degrees[vertex] -= 1
        if degrees[vertex] == 1:
            heappush(leaves, vertex)
    last = tuple(sorted((heappop(leaves), heappop(leaves))))
    edges.append(last)
    return tuple(edges)


def cayley_edge_census() -> tuple[int, dict[tuple[int, int], int]]:
    counts = {edge: 0 for edge in combinations(range(N_VERTICES), 2)}
    tree_count = 0
    for word in product(range(N_VERTICES), repeat=N_VERTICES - 2):
        edges = prufer_edges(word)
        require(len(edges) == N_VERTICES - 1, "Pruefer decoder edge count")
        require(len(set(edges)) == N_VERTICES - 1, "Pruefer decoder repeated edge")
        for edge in edges:
            counts[edge] += 1
        tree_count += 1
    return tree_count, counts


def maximum_spanning_tree(
    weights: dict[tuple[int, int], F], n: int = N_VERTICES
) -> tuple[F, tuple[tuple[F, int, int], ...]]:
    parent = list(range(n))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F(0)
    chosen: list[tuple[F, int, int]] = []
    for (i, j), weight in sorted(weights.items(), key=lambda row: (row[1], row[0]), reverse=True):
        first, second = root(i), root(j)
        if first == second:
            continue
        parent[first] = second
        total += weight
        chosen.append((weight, i, j))
        if len(chosen) == n - 1:
            break
    require(len(chosen) == n - 1, "Kruskal search did not span")
    return total, tuple(chosen)


def get_mixed(
    h: int, q: int, cache: dict[tuple[int, int], F]
) -> F:
    key = (h, q)
    if key not in cache:
        cache[key] = atomized_two_comb(h, GUARD_RADIUS, q, DANGER_RADIUS)
    value = mixed_overlap_formula(h, q)
    require(cache[key] == value, f"mixed formula mismatch at {(h, q)}")
    return value


def get_global(
    p: int, q: int, cache: dict[tuple[int, int], F]
) -> F:
    key = tuple(sorted((p, q)))
    if key not in cache:
        cache[key] = atomized_two_comb(p, DANGER_RADIUS, q, DANGER_RADIUS)
    value = global_pair_formula(p, q)
    require(cache[key] == value, f"global pair formula mismatch at {(p, q)}")
    return value


def get_triple(
    h: int,
    p: int,
    q: int,
    cache: dict[tuple[int, int, int], tuple[F, F]],
) -> tuple[F, F]:
    key = (h, min(p, q), max(p, q))
    if key not in cache:
        cache[key] = atomized_centered_triple(*key)
    return cache[key]


def exact_packet(
    Q: tuple[int, ...],
    h: int,
    edge_counts: dict[tuple[int, int], int],
    tree_count: int,
    mixed_cache: dict[tuple[int, int], F],
    global_cache: dict[tuple[int, int], F],
    triple_cache: dict[tuple[int, int, int], tuple[F, F]],
) -> dict[str, object]:
    intersections = tuple(get_mixed(h, q, mixed_cache) for q in Q)
    charges = tuple(value - F(2, 49) for value in intersections)
    deficit = -sum(charges, F(0))
    global_sum = F(0)
    energy = F(0)
    outside_sum = F(0)
    outside_weights: dict[tuple[int, int], F] = {}
    edge_identity_checks = 0
    for i, j in combinations(range(N_VERTICES), 2):
        rho = get_global(Q[i], Q[j], global_cache)
        outside, delta = get_triple(h, Q[i], Q[j], triple_cache)
        reconstructed = (
            F(5, 7) * rho - (charges[i] + charges[j]) / 7 - delta
        )
        require(
            outside == reconstructed,
            f"centered edge identity failed for Q={Q}, h={h}, edge={(i, j)}",
        )
        global_sum += rho
        energy += delta
        outside_sum += outside
        outside_weights[(i, j)] = outside
        edge_identity_checks += 1

    moment = atomized_centered_moment(Q, h)
    require(moment == energy, f"centered moment identity failed for Q={Q}, h={h}")

    cayley_weight = sum(
        (F(edge_counts[edge], tree_count) * outside_weights[edge] for edge in outside_weights),
        F(0),
    )
    require(cayley_weight == F(2, 7) * outside_sum, "Cayley average drift")
    average_margin = cayley_weight - deficit
    formula_margin = (
        F(10, 49) * global_sum
        + F(37, 49) * sum(charges, F(0))
        - F(2, 7) * energy
    )
    require(average_margin == formula_margin, "Cayley centered-energy identity failed")
    maximum_tree, tree_edges = maximum_spanning_tree(outside_weights)
    return {
        "charges": charges,
        "charge_sum": sum(charges, F(0)),
        "deficit": deficit,
        "global_sum": global_sum,
        "energy": energy,
        "outside_sum": outside_sum,
        "cayley_weight": cayley_weight,
        "average_margin": average_margin,
        "maximum_tree": maximum_tree,
        "maximum_margin": maximum_tree - deficit,
        "tree_edges": tree_edges,
        "edge_identity_checks": edge_identity_checks,
    }


def audit_charge_census() -> tuple[int, tuple[tuple[int, int, F], ...]]:
    # THM-2080's pointwise fold bound, checked on the complete residue grid.
    corrections = []
    for a_residue in range(14):
        for b_residue in range(7):
            x, y = F(a_residue, 14), F(b_residue, 7)
            correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
            corrections.append(correction)
    require(min(corrections) == -F(6, 49), "discrete fold correction floor changed")
    require(min(corrections) >= -F(1, 8), "analytic fold correction floor failed")

    rows: list[tuple[int, int, F]] = []
    reduced_count = 0
    for a in range(1, 30, 2):
        for b in range(1, 29 // a + 1):
            if a * b > 29 or gcd(a, b) != 1:
                continue
            charge = reduced_charge(a, b)
            if charge < 0:
                rows.append((a, b, -charge))
            reduced_count += 1
    rows.sort(key=lambda row: (-row[2], row[0], row[1]))
    expected_top = (
        (1, 6, F(5, 294)),
        (11, 1, F(8, 539)),
        (1, 5, F(3, 245)),
        (3, 5, F(3, 245)),
        (9, 1, F(4, 441)),
        (9, 2, F(4, 441)),
        (11, 2, F(9, 1078)),
    )
    require(tuple(rows[:7]) == expected_top, "top seven odd-a charges changed")
    require(sum((row[2] for row in rows[:7]), F(0)) == F(41, 495), "charge sum")
    require(F(1, 120) < rows[6][2], "ab>=30 tail can enter top seven")
    return reduced_count, tuple(rows)


def audit_pointwise_moment() -> int:
    checks = 0
    for multiplicity in range(N_VERTICES + 1):
        flags = [1] * multiplicity + [0] * (N_VERTICES - multiplicity)
        pair_sum = sum(
            (
                (F(flags[i]) - F(1, 7))
                * (F(flags[j]) - F(1, 7))
                for i, j in combinations(range(N_VERTICES), 2)
            ),
            F(0),
        )
        polynomial = F(1, 2) * (
            multiplicity * multiplicity
            - F(19, 7) * multiplicity
            + F(6, 7)
        )
        require(pair_sum == polynomial, f"pointwise moment failed at N={multiplicity}")
        checks += 1
    return checks


def main() -> None:
    require(
        F(10, 49) * PAIR_FLOOR + F(37, 49) * CHARGE_FLOOR == BASE_MARGIN,
        "universal base constant changed",
    )
    require(F(7, 2) * BASE_MARGIN == ENERGY_THRESHOLD, "energy threshold changed")

    tree_count, edge_counts = cayley_edge_census()
    expected_occurrences = 2 * N_VERTICES ** (N_VERTICES - 3)
    require(tree_count == N_VERTICES ** (N_VERTICES - 2), "Cayley tree count")
    require(set(edge_counts.values()) == {expected_occurrences}, "Cayley edge incidence")
    require(F(expected_occurrences, tree_count) == F(2, 7), "Cayley edge probability")

    reduced_count, negative_rows = audit_charge_census()
    moment_point_checks = audit_pointwise_moment()
    computed_modular_thresholds = tuple(
        F(7, 2)
        * (
            F(10, 49) * PAIR_FLOOR
            - F(37, 49)
            * sum((row[2] for row in negative_rows[: 7 - k]), F(0))
        )
        for k in range(5)
    )
    require(
        computed_modular_thresholds == MODULAR_THRESHOLDS,
        "7-divisible modular thresholds changed",
    )

    mixed_cache: dict[tuple[int, int], F] = {}
    global_cache: dict[tuple[int, int], F] = {}
    triple_cache: dict[tuple[int, int, int], tuple[F, F]] = {}
    core_count = 0
    pair_count = 0
    scalar_closures = 0
    scalar_survivors = 0
    average_closures = 0
    threshold_closures = 0
    edge_identity_checks = 0
    moment_checks = 0
    modular_zero_checks = 0
    modular_charge_checks = 0
    minimum_average = None

    for Q in combinations(range(1, MAX_CORE_SPEED + 1), N_VERTICES):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        core_count += 1
        for h in range(1, (8 * max(Q) - 1) // 3 + 1, 2):
            if 6 * h >= 16 * max(Q):
                continue
            pair_count += 1
            intersections = tuple(get_mixed(h, q, mixed_cache) for q in Q)
            charges = tuple(value - F(2, 49) for value in intersections)
            deficit = -sum(charges, F(0))
            require(sum(charges, F(0)) >= CHARGE_FLOOR, "seven-charge floor failed")
            if h % 7 != 0:
                seven_multiplicity = sum(q % 7 == 0 for q in Q)
                for q, charge in zip(Q, charges):
                    if q % 7 == 0:
                        require(charge == 0, f"7-divisible speed has charge at {(Q, h, q)}")
                        modular_zero_checks += 1
                if seven_multiplicity <= 4:
                    modular_cap = sum(
                        (row[2] for row in negative_rows[: 7 - seven_multiplicity]),
                        F(0),
                    )
                    require(
                        sum(charges, F(0)) >= -modular_cap,
                        f"modular charge floor failed at {(Q, h)}",
                    )
                    modular_charge_checks += 1
            if deficit < 0:
                scalar_closures += 1
                continue
            scalar_survivors += 1
            packet = exact_packet(
                Q,
                h,
                edge_counts,
                tree_count,
                mixed_cache,
                global_cache,
                triple_cache,
            )
            require(packet["global_sum"] >= PAIR_FLOOR, "seven-pair floor failed")
            require(packet["average_margin"] > 0, f"average tree failed at {(Q, h)}")
            average_closures += 1
            edge_identity_checks += int(packet["edge_identity_checks"])
            moment_checks += 1
            if packet["energy"] < ENERGY_THRESHOLD:
                threshold_closures += 1
            row = (
                packet["average_margin"],
                Q,
                h,
                packet["energy"],
                packet["global_sum"],
                packet["charge_sum"],
            )
            if minimum_average is None or row < minimum_average:
                minimum_average = row

    require(core_count == 131, f"terminal core count changed: {core_count}")
    require(pair_count == 4120, f"terminal pair count changed: {pair_count}")
    require(scalar_closures == 2798, f"scalar closure count changed: {scalar_closures}")
    require(scalar_survivors == 1322, "scalar survivor count changed")
    require(average_closures == scalar_survivors, "average certificate left a survivor")
    require(threshold_closures == 29, "threshold-only closure count changed")
    require(moment_checks == scalar_survivors, "moment audit count changed")
    require(edge_identity_checks == 27762, "edge identity count changed")
    require(minimum_average is not None, "empty terminal bank")
    expected_minimum = (
        F(99413, 2461368),
        (1, 9, 10, 11, 13, 14, 24),
        23,
        F(5448845, 40612572),
        F(25327, 63063),
        -F(49781, 11603592),
    )
    require(minimum_average == expected_minimum, "minimum average row changed")

    hostile_Q = (4, 5, 6, 7, 11, 13, 27)
    hostile_h = 1
    hostile = exact_packet(
        hostile_Q,
        hostile_h,
        edge_counts,
        tree_count,
        mixed_cache,
        global_cache,
        triple_cache,
    )
    require(hostile["charge_sum"] == -F(181007, 3783780), "hostile charges")
    require(hostile["global_sum"] == F(1716859, 3783780), "hostile pair sum")
    require(hostile["energy"] == F(29135, 147147), "hostile energy")
    require(hostile["deficit"] == F(181007, 3783780), "hostile deficit")
    require(hostile["cayley_weight"] == F(57481, 1203930), "hostile average weight")
    require(hostile["average_margin"] == -F(2467, 26486460), "hostile average margin")
    require(hostile["maximum_tree"] == F(118147, 1261260), "hostile maximum tree")
    require(hostile["maximum_margin"] == F(86717, 1891890), "hostile maximum margin")
    require(hostile["energy"] > ENERGY_THRESHOLD, "hostile energy should clear threshold")
    require(hostile["average_margin"] < 0 < hostile["maximum_margin"], "average/max guardrail")

    print("THM-2091 CENTERED TRIPLE ENERGY EXACT REFEREE")
    print(f"Cayley labelled K7 trees: {tree_count}")
    print(f"occurrences of each edge: {expected_occurrences}; probability 2/7")
    print(f"odd-a reduced pairs with ab<=29: {reduced_count}")
    print(f"negative-charge rows: {len(negative_rows)}")
    print("top seven negative charges:")
    for a, b, charge in negative_rows[:7]:
        print(f"  (a,b)=({a},{b}) -> {charge}")
    print(f"top-seven charge sum: {sum((row[2] for row in negative_rows[:7]), F(0))}")
    print(f"universal average-margin base: {BASE_MARGIN}")
    print(f"necessary containment energy: D >= {ENERGY_THRESHOLD}")
    print("7-divisible-speed thresholds k=0..4:")
    print("  " + ", ".join(str(value) for value in MODULAR_THRESHOLDS))
    print(f"pointwise multiplicity-moment checks: {moment_point_checks}")
    print(f"hereditary divisor-complete cores through 24: {core_count}")
    print(f"allowed odd terminal core/guard pairs: {pair_count}")
    print(f"negative-deficit scalar closures: {scalar_closures}")
    print(f"scalar survivors sent to centered energy: {scalar_survivors}")
    print(f"positive average-tree closures: {average_closures}")
    print(f"universal-threshold-only closures: {threshold_closures}")
    print(f"centered edge identity checks: {edge_identity_checks}")
    print(f"atomized multiplicity-moment checks: {moment_checks}")
    print(f"7-divisible zero-charge checks: {modular_zero_checks}")
    print(f"modular charge-floor packet checks: {modular_charge_checks}")
    print(
        "exact atom caches (mixed/global/triple): "
        f"{len(mixed_cache)}/{len(global_cache)}/{len(triple_cache)}"
    )
    print(
        "minimum terminal average margin: "
        f"{minimum_average[0]} at Q={minimum_average[1]}, h={minimum_average[2]}"
    )
    print(
        "minimum-row (D,R,sum e): "
        f"{minimum_average[3]}, {minimum_average[4]}, {minimum_average[5]}"
    )
    print(f"average-vs-max hostile control: h={hostile_h}, Q={hostile_Q}")
    print(
        "hostile (D,Delta,average margin,maximum margin): "
        f"{hostile['energy']}, {hostile['deficit']}, "
        f"{hostile['average_margin']}, {hostile['maximum_margin']}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
