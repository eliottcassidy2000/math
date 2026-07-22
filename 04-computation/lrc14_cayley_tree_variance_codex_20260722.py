#!/usr/bin/env python3
"""Exact referee for the Cayley-tree variance refinement of relative Hunter.

The script is dependency-free.  It verifies the labelled-tree inclusion
probabilities directly from all Prufer words, derives the exact variance from
those counts, replays the complete THM-2091 terminal bank through height 24,
and checks the hostile row on which the mean tree fails but the second-moment
ratio succeeds.  All arithmetic after interval atomization is rational.
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


def atomized_outside_overlap(h: int, p: int, q: int) -> F:
    boundaries = {F(0), F(1)}
    for speed, radius in (
        (h, GUARD_RADIUS),
        (p, DANGER_RADIUS),
        (q, DANGER_RADIUS),
    ):
        boundaries.update(comb_boundaries(speed, radius))
    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if (
            not in_comb(h, GUARD_RADIUS, midpoint)
            and in_comb(p, DANGER_RADIUS, midpoint)
            and in_comb(q, DANGER_RADIUS, midpoint)
        ):
            measure += right - left
    return measure


def mixed_overlap_formula(h: int, q: int) -> F:
    common = gcd(h, q)
    a, b = h // common, q // common
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, 49) + F(2, a * b) * correction


def global_pair_formula(p: int, q: int) -> F:
    common = gcd(p, q)
    a, b = p // common, q // common

    def fold(residue: int) -> int:
        value = residue % 14
        return value * (14 - value)

    return F(1, 49) + F(fold(a + b) - fold(b - a), 196 * a * b)


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def prufer_edges(
    word: tuple[int, ...], n: int = N_VERTICES
) -> tuple[tuple[int, int], ...]:
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
    edges.append(tuple(sorted((heappop(leaves), heappop(leaves)))))
    return tuple(edges)


def cayley_census() -> tuple[
    tuple[tuple[tuple[int, int], ...], ...],
    dict[tuple[int, int], int],
    dict[tuple[tuple[int, int], tuple[int, int]], int],
]:
    edges = tuple(combinations(range(N_VERTICES), 2))
    edge_counts = {edge: 0 for edge in edges}
    pair_counts = {pair: 0 for pair in combinations(edges, 2)}
    trees: list[tuple[tuple[int, int], ...]] = []
    for word in product(range(N_VERTICES), repeat=N_VERTICES - 2):
        tree = prufer_edges(word)
        require(len(tree) == N_VERTICES - 1, "Prufer edge count")
        require(len(set(tree)) == N_VERTICES - 1, "Prufer repeated edge")
        trees.append(tree)
        for edge in tree:
            edge_counts[edge] += 1
        for pair in combinations(sorted(tree), 2):
            pair_counts[pair] += 1
    return tuple(trees), edge_counts, pair_counts


def maximum_spanning_tree(
    weights: dict[tuple[int, int], F], n: int = N_VERTICES
) -> F:
    parent = list(range(n))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F(0)
    edge_count = 0
    for (i, j), weight in sorted(
        weights.items(), key=lambda row: (row[1], row[0]), reverse=True
    ):
        first, second = root(i), root(j)
        if first == second:
            continue
        parent[first] = second
        total += weight
        edge_count += 1
        if edge_count == n - 1:
            break
    require(edge_count == n - 1, "Kruskal search did not span")
    return total


def get_mixed(
    h: int, q: int, cache: dict[tuple[int, int], F]
) -> F:
    key = (h, q)
    if key not in cache:
        exact = atomized_two_comb(h, GUARD_RADIUS, q, DANGER_RADIUS)
        formula = mixed_overlap_formula(h, q)
        require(exact == formula, f"mixed formula mismatch at {key}")
        cache[key] = exact
    return cache[key]


def get_global(
    p: int, q: int, cache: dict[tuple[int, int], F]
) -> F:
    key = tuple(sorted((p, q)))
    if key not in cache:
        exact = atomized_two_comb(p, DANGER_RADIUS, q, DANGER_RADIUS)
        formula = global_pair_formula(p, q)
        require(exact == formula, f"global formula mismatch at {key}")
        cache[key] = exact
    return cache[key]


def get_outside(
    h: int,
    p: int,
    q: int,
    cache: dict[tuple[int, int, int], F],
) -> F:
    key = (h, min(p, q), max(p, q))
    if key not in cache:
        cache[key] = atomized_outside_overlap(*key)
    return cache[key]


def cayley_mean_variance(
    weights: dict[tuple[int, int], F]
) -> tuple[F, F, F]:
    total = sum(weights.values(), F(0))
    mean = F(2, 7) * total
    square_sum = sum((weight * weight for weight in weights.values()), F(0))
    vertex_sums = [
        sum(
            (
                weight
                for (i, j), weight in weights.items()
                if i == vertex or j == vertex
            ),
            F(0),
        )
        for vertex in range(N_VERTICES)
    ]
    variance = (
        12 * square_sum - sum((value * value for value in vertex_sums), F(0))
    ) / 49
    require(variance >= 0, "negative Cayley variance")
    second_ratio = mean + variance / mean if mean > 0 else F(0)
    return mean, variance, second_ratio


def exact_packet(
    Q: tuple[int, ...],
    h: int,
    mixed_cache: dict[tuple[int, int], F],
    global_cache: dict[tuple[int, int], F],
    outside_cache: dict[tuple[int, int, int], F],
) -> dict[str, object]:
    mixed = tuple(get_mixed(h, q, mixed_cache) for q in Q)
    charges = tuple(value - F(2, 49) for value in mixed)
    charge_sum = sum(charges, F(0))
    deficit = -charge_sum
    weights = {
        (i, j): get_outside(h, Q[i], Q[j], outside_cache)
        for i, j in combinations(range(N_VERTICES), 2)
    }
    mean, variance, second_ratio = cayley_mean_variance(weights)
    global_sum = sum(
        (get_global(Q[i], Q[j], global_cache) for i, j in combinations(range(7), 2)),
        F(0),
    )
    outside_sum = sum(weights.values(), F(0))
    energy = F(5, 7) * global_sum - F(6, 7) * charge_sum - outside_sum
    average_margin = mean - deficit
    require(mean > 0, "zero Cayley mean on a scalar survivor")
    formula_margin = (
        F(10, 49) * global_sum
        + F(37, 49) * charge_sum
        - F(2, 7) * energy
    )
    require(average_margin == formula_margin, "centered mean identity")
    exact_energy_floor = (
        F(5, 7) * global_sum
        + F(37, 14) * charge_sum
        + F(7, 2) * variance / mean
    )
    variance_margin = second_ratio - deficit
    require(
        variance_margin == F(2, 7) * (exact_energy_floor - energy),
        "variance/energy threshold algebra",
    )
    maximum = maximum_spanning_tree(weights)
    require(second_ratio <= maximum, "second-moment ratio exceeded maximum tree")
    return {
        "weights": weights,
        "charge_sum": charge_sum,
        "deficit": deficit,
        "global_sum": global_sum,
        "energy": energy,
        "mean": mean,
        "variance": variance,
        "second_ratio": second_ratio,
        "average_margin": average_margin,
        "variance_margin": variance_margin,
        "maximum": maximum,
        "maximum_margin": maximum - deficit,
    }


def main() -> None:
    require(
        F(10, 49) * PAIR_FLOOR + F(37, 49) * CHARGE_FLOOR
        == BASE_MARGIN,
        "universal base margin",
    )
    require(F(7, 2) * BASE_MARGIN == ENERGY_THRESHOLD, "energy threshold")

    trees, edge_counts, pair_counts = cayley_census()
    tree_count = len(trees)
    require(tree_count == 7**5, "Cayley tree count")
    require(set(edge_counts.values()) == {2 * 7**4}, "single-edge incidence")
    adjacent_counts = {
        count for (first, second), count in pair_counts.items() if set(first) & set(second)
    }
    disjoint_counts = {
        count for (first, second), count in pair_counts.items() if not set(first) & set(second)
    }
    require(adjacent_counts == {3 * 7**3}, "adjacent-pair incidence")
    require(disjoint_counts == {4 * 7**3}, "disjoint-pair incidence")
    require(
        F(next(iter(adjacent_counts)), tree_count) == F(3, 49),
        "adjacent-pair probability",
    )
    require(
        F(next(iter(disjoint_counts)), tree_count) == F(4, 49),
        "disjoint-pair probability",
    )

    mixed_cache: dict[tuple[int, int], F] = {}
    global_cache: dict[tuple[int, int], F] = {}
    outside_cache: dict[tuple[int, int, int], F] = {}
    core_count = 0
    pair_count = 0
    scalar_closures = 0
    scalar_survivors = 0
    average_closures = 0
    variance_closures = 0
    threshold_closures = 0
    refined_threshold_closures = 0
    positive_variances = 0
    minimum_average = None
    minimum_variance = None
    minimum_maximum_slack = None

    for Q in combinations(range(1, MAX_CORE_SPEED + 1), N_VERTICES):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        core_count += 1
        for h in range(1, (8 * max(Q) - 1) // 3 + 1, 2):
            if 6 * h >= 16 * max(Q):
                continue
            pair_count += 1
            charge_sum = sum(
                (get_mixed(h, q, mixed_cache) - F(2, 49) for q in Q), F(0)
            )
            require(charge_sum >= CHARGE_FLOOR, "charge floor failed")
            if -charge_sum < 0:
                scalar_closures += 1
                continue
            scalar_survivors += 1
            packet = exact_packet(
                Q, h, mixed_cache, global_cache, outside_cache
            )
            require(packet["global_sum"] >= PAIR_FLOOR, "pair floor failed")
            if packet["average_margin"] > 0:
                average_closures += 1
            if packet["variance_margin"] > 0:
                variance_closures += 1
            if packet["variance"] > 0:
                positive_variances += 1
            if packet["energy"] < ENERGY_THRESHOLD:
                threshold_closures += 1
            refined_threshold = (
                ENERGY_THRESHOLD
                + F(7, 2) * packet["variance"] / packet["mean"]
            )
            if packet["energy"] < refined_threshold:
                refined_threshold_closures += 1
            average_row = (packet["average_margin"], Q, h)
            variance_row = (packet["variance_margin"], Q, h)
            maximum_slack_row = (
                packet["maximum"] - packet["second_ratio"], Q, h
            )
            if minimum_average is None or average_row < minimum_average:
                minimum_average = average_row
            if minimum_variance is None or variance_row < minimum_variance:
                minimum_variance = variance_row
            if (
                minimum_maximum_slack is None
                or maximum_slack_row < minimum_maximum_slack
            ):
                minimum_maximum_slack = maximum_slack_row

    require(core_count == 131, f"terminal core count changed: {core_count}")
    require(pair_count == 4120, f"terminal pair count changed: {pair_count}")
    require(scalar_closures == 2798, "scalar closure count changed")
    require(scalar_survivors == 1322, "scalar survivor count changed")
    require(average_closures == 1322, "average closure count changed")
    require(variance_closures == 1322, "variance closure count changed")
    require(positive_variances == 1322, "zero variance entered bounded bank")
    require(threshold_closures == 29, "threshold-only closure count changed")
    require(
        refined_threshold_closures == 69,
        "variance-refined threshold closure count changed",
    )
    require(minimum_average is not None, "empty average ledger")
    require(minimum_variance is not None, "empty variance ledger")
    require(minimum_maximum_slack is not None, "empty maximum ledger")
    require(
        minimum_average
        == (F(99413, 2461368), (1, 9, 10, 11, 13, 14, 24), 23),
        "minimum average row changed",
    )
    require(
        minimum_variance
        == (
            F(22391180465191, 526381866410400),
            (1, 9, 10, 11, 13, 14, 24),
            23,
        ),
        "minimum variance row changed",
    )
    require(
        minimum_maximum_slack
        == (
            F(3701706622037, 313124176749384),
            (9, 10, 11, 13, 14, 17, 24),
            49,
        ),
        "minimum maximum-slack row changed",
    )

    hostile_Q = (4, 5, 6, 7, 11, 13, 27)
    hostile_h = 1
    hostile = exact_packet(
        hostile_Q, hostile_h, mixed_cache, global_cache, outside_cache
    )
    require(hostile["mean"] == F(57481, 1203930), "hostile mean")
    require(hostile["deficit"] == F(181007, 3783780), "hostile deficit")
    require(
        hostile["average_margin"] == -F(2467, 26486460),
        "hostile average margin",
    )
    require(
        hostile["variance"] == F(919435861, 5314640631300),
        "hostile variance",
    )
    require(
        hostile["variance_margin"] == F(488619049, 138406200660),
        "hostile variance margin",
    )
    require(
        hostile["maximum_margin"] == F(86717, 1891890),
        "hostile maximum margin",
    )
    require(hostile["energy"] == F(29135, 147147), "hostile energy")
    require(
        hostile["average_margin"] < 0 < hostile["variance_margin"],
        "hostile mean/variance split",
    )

    hostile_tree_weights = tuple(
        sum((hostile["weights"][edge] for edge in tree), F(0))
        for tree in trees
    )
    empirical_mean = sum(hostile_tree_weights, F(0)) / tree_count
    empirical_second = (
        sum((value * value for value in hostile_tree_weights), F(0)) / tree_count
    )
    require(empirical_mean == hostile["mean"], "empirical hostile mean")
    require(
        empirical_second - empirical_mean * empirical_mean == hostile["variance"],
        "empirical hostile variance",
    )
    require(max(hostile_tree_weights) == hostile["maximum"], "hostile tree maximum")

    print("CAYLEY-TREE VARIANCE EXACT REFEREE")
    print(f"labelled K7 trees: {tree_count}")
    print(f"single-edge occurrences: {next(iter(set(edge_counts.values())))} = probability 2/7")
    print(
        "two-edge occurrences: "
        f"adjacent {next(iter(adjacent_counts))} = 3/49; "
        f"disjoint {next(iter(disjoint_counts))} = 4/49"
    )
    print("variance formula: (12 sum_e w_e^2 - sum_i s_i^2)/49")
    print("second-moment tree bound: tau >= mu + Var/mu")
    print(
        "refined universal energy obstruction: "
        "D >= 2059/90090 + (7/2) Var/mu"
    )
    print(f"hereditary divisor-complete cores through 24: {core_count}")
    print(f"allowed odd terminal core/guard pairs: {pair_count}")
    print(f"negative-deficit scalar closures: {scalar_closures}")
    print(f"scalar survivors: {scalar_survivors}")
    print(f"positive average-tree closures: {average_closures}")
    print(f"positive variance-tree closures: {variance_closures}")
    print(f"strictly positive tree variances: {positive_variances}")
    print(f"universal energy-threshold closures: {threshold_closures}")
    print(
        "variance-refined universal-threshold closures: "
        f"{refined_threshold_closures}"
    )
    print(
        "minimum average margin: "
        f"{minimum_average[0]} at Q={minimum_average[1]}, h={minimum_average[2]}"
    )
    print(
        "minimum variance margin: "
        f"{minimum_variance[0]} at Q={minimum_variance[1]}, h={minimum_variance[2]}"
    )
    print(
        "minimum (maximum tree - second ratio): "
        f"{minimum_maximum_slack[0]} at "
        f"Q={minimum_maximum_slack[1]}, h={minimum_maximum_slack[2]}"
    )
    print(
        "exact atom caches (mixed/global/outside): "
        f"{len(mixed_cache)}/{len(global_cache)}/{len(outside_cache)}"
    )
    print(f"hostile row: h={hostile_h}, Q={hostile_Q}")
    print(
        "hostile (mu, Delta, Var): "
        f"{hostile['mean']}, {hostile['deficit']}, {hostile['variance']}"
    )
    print(
        "hostile (average margin, variance margin, maximum margin): "
        f"{hostile['average_margin']}, {hostile['variance_margin']}, "
        f"{hostile['maximum_margin']}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
