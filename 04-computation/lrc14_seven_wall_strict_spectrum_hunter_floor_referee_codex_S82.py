#!/usr/bin/env python3
"""Exact referee for THM-1221, the sharpened seven-wall Hunter floor.

The script is dependency-free and uses ``Fraction`` throughout.  It replays
the THM-1166 reduced pair spectrum and sharp three-speed sum provider, proves
the finite low-graph facts used by the new 1+6-cut argument, and checks every
branch of the spanning-tree consumer.  No bounded speed box is used as a
substitute for an analytic tail: the two product cutoffs are asserted and
checked explicitly before their finite banks are inspected.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path


ONE = F(1)
PAIR_FLOOR = F(1, 91)
HIGH_BAR = F(1, 63)
TRIPLE_SUM = F(51, 1183)
NEXT_HIGH = F(5, 308)
HEAVY_CIRCUIT_CUTOFF = F(60, 637)
ALIGNMENT_FLOOR = F(265, 2772)
TREE_FLOOR = F(15, 154)

CROSS_A = F(1, 70)
CROSS_B = F(1, 77)
CROSS_C = F(1, 84)
CROSS_D = F(1, 91)
SPECIAL_A_EDGE = F(32, 1575)


def require(condition: bool, message: str) -> None:
    """Optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def fmt(value: F) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fold(residue: int, modulus: int) -> int:
    residue %= modulus
    return residue * (modulus - residue)


def tent_sum(z: F) -> F:
    q = z.numerator // z.denominator
    return (2 * q + 1) * z - q * (q + 1)


def rho_fold(a: int, b: int) -> F:
    """Haar measure of D_a intersect D_b at radius 1/14."""
    require(a > 0 and b > 0, "pair speeds must be positive")
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    modulus = 14 * g
    numerator = 4 * a * b + fold(a + b, modulus) - fold(b - a, modulus)
    return F(numerator, 196 * a * b)


def rho_tent(a: int, b: int) -> F:
    """Independent reduced tent evaluation of the same pair mass."""
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    a //= g
    b //= g
    return (tent_sum(F(a + b, 14)) - tent_sum(F(b - a, 14))) / (a * b)


def rho_ratio(x: F, y: F) -> F:
    ratio = x / y
    return rho_fold(abs(ratio.numerator), abs(ratio.denominator))


def reduced_bank(product_cap: int, predicate) -> list[tuple[int, int, F]]:
    return [
        (a, b, rho_fold(a, b))
        for a in range(1, product_cap + 1)
        for b in range(a + 1, product_cap + 1)
        if a * b <= product_cap and gcd(a, b) == 1 and predicate(rho_fold(a, b))
    ]


def oriented_reduced_ratios(product_cap: int) -> list[F]:
    return sorted(
        {
            F(a, b)
            for a in range(1, product_cap + 1)
            for b in range(1, product_cap + 1)
            if a != b and a * b <= product_cap and gcd(a, b) == 1
        }
    )


def maximum_spanning_tree(values: tuple[F, ...]) -> F:
    n = len(values)
    edges = sorted(
        ((rho_ratio(values[i], values[j]), i, j) for i in range(n) for j in range(i + 1, n)),
        reverse=True,
    )
    parent = list(range(n))

    def root(v: int) -> int:
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    total = F(0)
    used = 0
    for weight, i, j in edges:
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
            total += weight
            used += 1
            if used == n - 1:
                return total
    raise RuntimeError("complete pair graph did not connect")


def graph_connected(values: tuple[F, ...], predicate) -> bool:
    """Connectivity of the pair graph selected by ``predicate(weight)``."""
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(len(values)):
            if j not in seen and predicate(rho_ratio(values[i], values[j])):
                seen.add(j)
                stack.append(j)
    return len(seen) == len(values)


def tournament_fingerprint(values: tuple[F, ...]) -> tuple:
    """Low edges point down the speed order; all other edges point up."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    rank = {vertex: i for i, vertex in enumerate(order)}
    edge: set[tuple[int, int]] = set()
    low_flips = 0
    for i in range(n):
        for j in range(i + 1, n):
            small, large = (i, j) if rank[i] < rank[j] else (j, i)
            if rho_ratio(values[i], values[j]) < HIGH_BAR:
                edge.add((large, small))
                low_flips += 1
            else:
                edge.add((small, large))

    scores = tuple(sorted(sum((i, j) in edge for j in range(n) if j != i) for i in range(n)))
    cycles = sum(
        (i, j) in edge and (j, k) in edge and (k, i) in edge
        for i, j, k in permutations(range(n), 3)
    ) // 3

    path: list[int] = []
    for vertex in range(n):
        pos = 0
        while pos < len(path) and (path[pos], vertex) in edge:
            pos += 1
        path.insert(pos, vertex)
    require(all((path[i], path[i + 1]) in edge for i in range(n - 1)), "bad Redei path")
    ham_count = sum(
        all((p[i], p[i + 1]) in edge for i in range(n - 1))
        for p in permutations(range(n))
    )

    reach = [[i == j or (i, j) in edge for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc_sizes: list[int] = []
    while unseen:
        v = min(unseen)
        component = {w for w in unseen if reach[v][w] and reach[w][v]}
        unseen -= component
        scc_sizes.append(len(component))
    return scores, cycles, tuple(sorted(scc_sizes)), low_flips, tuple(path), ham_count


def main() -> None:
    print("THM-1221 sharpened seven-wall strict-spectrum Hunter-floor referee")

    # Formula audit, including nonprimitive pairs.
    formula_checks = 0
    for a in range(1, 121):
        for b in range(a + 1, 121):
            require(rho_fold(a, b) == rho_tent(a, b), f"pair formula mismatch at {(a, b)}")
            formula_checks += 1
    print(f"folded-vs-tent pair checks: {formula_checks}")

    # THM-1166 strict and closed low-channel providers.  The folded defect
    # bound rho>=1/49-1/(4ab) puts every product ab>=56 strictly above 1/63.
    require(F(1, 49) - F(1, 4 * 56) > HIGH_BAR, "low-channel tail cutoff failed")
    strict_low = reduced_bank(55, lambda value: value < HIGH_BAR)
    expected_strict = [
        (1, 10, CROSS_A),
        (1, 11, CROSS_B),
        (1, 12, CROSS_C),
        (1, 13, CROSS_D),
        (2, 11, CROSS_B),
        (3, 10, CROSS_A),
        (3, 11, CROSS_B),
    ]
    require(strict_low == expected_strict, "strict-low bank mismatch")

    closed_low = reduced_bank(55, lambda value: value <= HIGH_BAR)
    expected_closed = [
        (1, 9, HIGH_BAR),
        (1, 10, CROSS_A),
        (1, 11, CROSS_B),
        (1, 12, CROSS_C),
        (1, 13, CROSS_D),
        (1, 27, HIGH_BAR),
        (2, 9, HIGH_BAR),
        (2, 11, CROSS_B),
        (3, 10, CROSS_A),
        (3, 11, CROSS_B),
        (4, 9, HIGH_BAR),
        (5, 9, HIGH_BAR),
    ]
    require(closed_low == expected_closed, "closed-low bank mismatch")
    print("analytic low-channel tail: ab>=56 implies rho>1/63")
    print("strict-low reduced channels (7): " + ", ".join(f"({a},{b}):{fmt(v)}" for a, b, v in strict_low))
    print("low-or-equal reduced channels: 12 (five equality channels)")

    # Sharp three-speed aggregate floor, replayed rather than merely trusted.
    first_threshold = F(3, 4) / (3 * F(1, 49) - TRIPLE_SUM)
    second_threshold = F(1, 2) / (PAIR_FLOOR + 2 * F(1, 49) - TRIPLE_SUM)
    require(first_threshold == F(8281, 200) and 41 < first_threshold < 42, "first triple cutoff")
    require(second_threshold == F(8281, 144) and 57 < second_threshold < 58, "second triple cutoff")
    first_ratios = oriented_reduced_ratios(41)
    second_ratios = oriented_reduced_ratios(57)
    candidates = [
        (
            rho_ratio(first, ONE) + rho_ratio(second, ONE) + rho_ratio(first, second),
            first,
            second,
        )
        for first in first_ratios
        for second in second_ratios
        if first != second
    ]
    require(len(first_ratios) == 124 and len(second_ratios) == 184, "ratio bank sizes")
    require(len(candidates) == 22692, "sharp triple candidate count")
    triple_min = min(value for value, _, _ in candidates)
    triple_minimizers = {(x, y) for value, x, y in candidates if value == triple_min}
    require(triple_min == TRIPLE_SUM, "sharp triple sum mismatch")
    require(triple_minimizers == {(F(1, 13), F(13)), (F(13), F(1, 13))}, "triple minimizers")
    print("sharp triple provider: 124*184-124=22692 exact configurations")
    print("sharp three-speed sum: 51/1183, unique ratio triple (1,13,169)")

    # The graph rho<=1/63 has clique number THREE, not two.  There are two
    # reciprocal-dual triangle shapes, (1,12,27) and (4,9,108); neither has a
    # fourth common low-or-equal neighbor.  This is all the connected-branch
    # proof needs for a seven-packet.
    closed_positive = {F(b, a) for a, b, _ in closed_low}
    closed_products = sorted(
        (r, s, r * s) for r in closed_positive for s in closed_positive if r * s in closed_positive
    )
    require(
        closed_products == [(F(9, 4), F(12), F(27)), (F(12), F(9, 4), F(27))],
        "low-or-equal triangle classification mismatch",
    )
    closed_oriented = closed_positive | {1 / r for r in closed_positive}
    closed_triangles = ((F(1), F(12), F(27)), (F(4), F(9), F(108)))
    for triangle in closed_triangles:
        common_to_triangle = set.intersection(
            *({x * r for r in closed_oriented} for x in triangle)
        )
        require(not common_to_triangle, f"closed triangle {triangle} has a fourth neighbor")
    print(
        "rho<=1/63 graph: clique number 3; reciprocal-dual triangles "
        "(1,12,27), (4,9,108); no K4"
    )

    # First spectral value strictly above the bar.  Products >=60 are handled
    # analytically; products <=59 are the complete finite core.
    require(F(1, 49) - F(1, 4 * 60) > NEXT_HIGH, "next-high analytic tail cutoff")
    above_core = reduced_bank(59, lambda value: value > HIGH_BAR)
    next_value = min(value for _, _, value in above_core)
    next_rows = [(a, b) for a, b, value in above_core if value == next_value]
    require(next_value == NEXT_HIGH and next_rows == [(4, 11)], "next-high spectrum mismatch")
    require(NEXT_HIGH + 5 * HIGH_BAR == ALIGNMENT_FLOOR, "aligned tree-floor identity")
    print("first pair mass above 1/63: 5/308, unique reduced ratio (4,11)")
    print("first connected-high ledger: 5/308+5/63=265/2772")

    # Strict-low triangle-freeness and the exact common-neighbor obstruction.
    strict_positive = {F(b, a) for a, b, _ in strict_low}
    require(
        not any(r * s in strict_positive for r in strict_positive for s in strict_positive),
        "strict-low triangle found",
    )
    strict_oriented = strict_positive | {1 / r for r in strict_positive}
    intersection_rows: dict[F, tuple[F, ...]] = {}
    for ratio in sorted({r / s for r in strict_oriented for s in strict_oriented} - {ONE}):
        common = tuple(sorted(z for z in strict_oriented if z / ratio in strict_oriented))
        if common:
            intersection_rows[ratio] = common
    histogram = Counter(map(len, intersection_rows.values()))
    require(histogram == Counter({2: 72, 1: 14, 4: 6}), "common-neighbor histogram")
    max_rows = {ratio: common for ratio, common in intersection_rows.items() if len(common) == 4}
    expected_max_ratios = {F(3, 110), F(1, 3), F(10, 11), F(11, 10), F(3), F(110, 3)}
    require(set(max_rows) == expected_max_ratios, "four-common-neighbor ratio list")
    for ratio, common in max_rows.items():
        centers = set.intersection(*({neighbor / r for r in strict_oriented} for neighbor in common))
        require(centers == {ONE, ratio}, f"third common center at ratio {ratio}")
    print("strict-low graph: triangle-free")
    print("distinct-center common-low-neighbor histogram: 1:14, 2:72, 4:6; maximum 4")
    print("every four-neighbor row has exactly its original two centers")
    print("therefore a disconnected high graph on seven vertices has only a 1+6 cut")

    # A singleton's 14 possible strict-low ratios and exact weight capacities.
    crossing_weights = Counter(rho_ratio(ONE, ratio) for ratio in strict_oriented)
    require(
        crossing_weights == Counter({CROSS_A: 4, CROSS_B: 6, CROSS_C: 2, CROSS_D: 2}),
        "singleton crossing capacities",
    )
    print("singleton crossing capacities: 1/70:4, 1/77:6, 1/84:2, 1/91:2")

    # If m<=3 crossings have weight 1/70, the other endpoints have crossing
    # mass <=1/77.  Their mutual edge has weight at least u by the triple
    # floor; all remaining internal K6 edges have weight at least 1/63.
    u = TRIPLE_SUM - 2 * CROSS_B
    require(u == F(223, 13013) and u > HIGH_BAR, "non-A internal credit")
    aligned_singleton_floors = {0: CROSS_B + 5 * u}
    for m in range(1, 4):
        aligned_singleton_floors[m] = CROSS_A + (5 - m) * u + m * HIGH_BAR
    expected_aligned_floors = {
        0: F(1284, 13013),
        1: F(115601, 1171170),
        2: F(16303, 167310),
        3: F(37547, 390390),
    }
    require(aligned_singleton_floors == expected_aligned_floors, "m<=3 aligned floor table")
    require(
        all(value > ALIGNMENT_FLOOR for value in aligned_singleton_floors.values()),
        "m<=3 misses aligned target",
    )

    # Four 1/70 crossings force all four oriented 1/70 ratios.  Every edge
    # among those four endpoints has the stronger exact floor 32/1575.
    a_ratios = sorted(ratio for ratio in strict_oriented if rho_ratio(ONE, ratio) == CROSS_A)
    require(a_ratios == [F(1, 10), F(3, 10), F(10, 3), F(10)], "1/70 ratio set")
    a_internal = [rho_ratio(x, y) for x, y in combinations(a_ratios, 2)]
    require(min(a_internal) == SPECIAL_A_EDGE, "four-A internal floor")
    aligned_m4_floor = CROSS_A + SPECIAL_A_EDGE + 4 * HIGH_BAR
    require(
        aligned_m4_floor == F(103, 1050) and aligned_m4_floor > ALIGNMENT_FLOOR,
        "m=4 aligned floor",
    )

    # Strong analytic singleton tree.  The m A-endpoints form a tree with
    # (m-1) edges of weight at least s=32/1575; the 6-m non-A endpoints form
    # a tree with 5-m edges of weight at least u.  One high connector joins
    # the endpoint trees, and an A crossing joins the singleton.  At m=0 use
    # the forced crossing B and the five-edge non-A tree.
    strong_singleton_floors = {0: CROSS_B + 5 * u}
    for m in range(1, 5):
        strong_singleton_floors[m] = (
            CROSS_A + (m - 1) * SPECIAL_A_EDGE + (5 - m) * u + HIGH_BAR
        )
    expected_strong_floors = {
        0: F(1284, 13013),
        1: F(115601, 1171170),
        2: F(28411, 278850),
        3: F(615257, 5855850),
        4: F(633883, 5855850),
    }
    require(strong_singleton_floors == expected_strong_floors, "strong singleton floor table")
    require(all(value > TREE_FLOOR for value in strong_singleton_floors.values()), "strong singleton floor")
    require(min(strong_singleton_floors.values()) == F(1284, 13013), "strong singleton minimum")
    require(F(1284, 13013) - TREE_FLOOR == F(3, 2366), "strong singleton margin")
    print(f"non-A mutual-edge credit: u=51/1183-2/77={fmt(u)} >1/63")
    for m in range(4):
        print(
            f"singleton aligned cut m={m}: tree>={fmt(aligned_singleton_floors[m])}; "
            f"margin over 265/2772={fmt(aligned_singleton_floors[m]-ALIGNMENT_FLOOR)}"
        )
    print(f"singleton cut m=4: four 1/70 ratios force internal rho>={fmt(SPECIAL_A_EDGE)}")
    print(
        f"singleton aligned cut m=4: tree>={fmt(aligned_m4_floor)}; "
        f"margin over 265/2772={fmt(aligned_m4_floor-ALIGNMENT_FLOOR)}"
    )
    for m in range(5):
        print(
            f"singleton strong cut m={m}: tree>={fmt(strong_singleton_floors[m])}; "
            f"margin over 15/154={fmt(strong_singleton_floors[m]-TREE_FLOOR)}"
        )

    # Exhaust all feasible six-ratio singleton packets as an independent
    # diagnostic.  The strong analytic ledger above already proves this branch.
    ledger_checks = 0
    exact_min: tuple[F, tuple[F, ...]] | None = None
    for chosen in combinations(sorted(strict_oriented), 6):
        m = sum(rho_ratio(ONE, ratio) == CROSS_A for ratio in chosen)
        ledger = strong_singleton_floors[m]
        require(ledger > TREE_FLOOR, f"ledger failure at {chosen}")
        exact = maximum_spanning_tree((ONE,) + chosen)
        require(exact >= ledger, f"actual MST below ledger at {chosen}")
        if exact_min is None or exact < exact_min[0]:
            exact_min = (exact, chosen)
        ledger_checks += 1
    require(ledger_checks == 3003 and exact_min is not None, "six-ratio census size")
    require(exact_min[0] > TREE_FLOOR, "disconnected high-graph census misses strong floor")
    print(f"all C(14,6) singleton ratio packets: {ledger_checks}; ledger violations: 0")
    print(f"diagnostic exact disconnected MST minimum: {fmt(exact_min[0])} at {exact_min[1]}")

    # Sharpen the connected rho>=1/63 branch.  Let H be the strict graph
    # rho>1/63.  Its components have all cross-weights <=1/63, so choosing
    # one point from each component embeds a clique in the closed-low graph.
    # The no-K4 calculation above gives at most three components.  The two
    # possible representative triangles are the reciprocal-dual shapes
    # above.  A second point in each component must be a common closed
    # neighbor of the other two; all six exact intersections have size two,
    # so seven vertices cannot be split among three strict-high components.
    opposite_common: list[list[tuple[F, ...]]] = []
    for triangle in closed_triangles:
        rows: list[tuple[F, ...]] = []
        for omitted in range(3):
            other = [triangle[i] for i in range(3) if i != omitted]
            common = tuple(
                sorted(
                    {other[0] * r for r in closed_oriented}
                    & {other[1] * r for r in closed_oriented}
                )
            )
            rows.append(common)
        opposite_common.append(rows)
    require(
        opposite_common
        == [
            [(F(1), F(324)), (F(9, 4), F(12)), (F(4, 9), F(27))],
            [(F(4), F(243)), (F(9), F(48)), (F(1, 3), F(108))],
        ],
        "three-component common-neighbor bank",
    )
    require(
        all(sum(len(row) for row in triangle_rows) == 6 for triangle_rows in opposite_common),
        "three components could hold seven vertices",
    )
    print(
        "strict-high component count <=3; both three-component shapes "
        "have at most 2+2+2 vertices"
    )

    # Two strict-high components give a complete bipartite closed-low cut.
    # Normalize one vertex to 1.  The three possible cut sizes are exhausted
    # below; scale invariance makes these finite classifications global.
    closed_ratios = tuple(sorted(closed_oriented))
    closed_internal_strict = [
        (rho_ratio(x, y), x, y)
        for x, y in combinations(closed_ratios, 2)
        if rho_ratio(x, y) > HIGH_BAR
    ]
    require(min(value for value, _, _ in closed_internal_strict) == F(1, 56), "closed-Q strict edge floor")

    center_ratios = sorted({r / s for r in closed_ratios for s in closed_ratios} - {ONE})
    closed_center_rows: dict[F, tuple[F, ...]] = {}
    for center_ratio in center_ratios:
        common = tuple(
            sorted(
                vertex
                for vertex in closed_ratios
                if vertex / center_ratio in closed_oriented
            )
        )
        if common:
            closed_center_rows[center_ratio] = common
    require(
        Counter(map(len, closed_center_rows.values()))
        == Counter({2: 198, 4: 26, 1: 22, 6: 4, 3: 2}),
        "closed common-neighbor histogram",
    )
    large_center_rows = {
        ratio: common for ratio, common in closed_center_rows.items() if len(common) >= 5
    }
    require(set(large_center_rows) == {F(1, 3), F(1, 2), F(2), F(3)}, "large center rows")
    require(
        min(rho_ratio(ONE, ratio) for ratio in large_center_rows) == F(1, 21),
        "large-row center edge floor",
    )

    triple_common_histogram: Counter[int] = Counter()
    for second, third in combinations(center_ratios, 2):
        centers = (ONE, second, third)
        if not graph_connected(centers, lambda weight: weight > HIGH_BAR):
            continue
        common = tuple(
            vertex
            for vertex in closed_ratios
            if vertex / second in closed_oriented and vertex / third in closed_oriented
        )
        if common:
            triple_common_histogram[len(common)] += 1
    require(
        triple_common_histogram == Counter({1: 5532, 2: 243, 3: 12}),
        "closed connected-triple common-neighbor histogram",
    )
    require(max(triple_common_histogram) == 3, "closed connected triple has four common neighbors")

    one_six_analytic = 5 * F(1, 56) + HIGH_BAR
    two_five_analytic = F(1, 21) + 4 * NEXT_HIGH + HIGH_BAR
    require(one_six_analytic == F(53, 504) and one_six_analytic > TREE_FLOOR, "analytic 1+6 floor")
    require(two_five_analytic == F(89, 693) and two_five_analytic > TREE_FLOOR, "analytic 2+5 floor")
    print("strict-high two-component analytic cut proof:")
    print("  1+6: strict edge floor inside closed-Q bank=1/56; tree>=5/56+1/63=53/504")
    print("  2+5: five common neighbors force center ratio in {1/3,1/2,2,3}")
    print("       center edge>=1/21; tree>=1/21+4*(5/308)+1/63=89/693")
    print("  3+4: connected normalized center triples have at most three common closed-Q neighbors")

    # The larger full cut census is an independent replay, not a proof
    # dependency after the preceding small-bank analytic construction.
    two_component_counts: dict[str, int] = {}
    two_component_minima: dict[str, tuple[F, tuple[F, ...]]] = {}

    valid_packets: set[tuple[F, ...]] = set()
    for right in combinations(closed_ratios, 6):
        if not graph_connected(right, lambda weight: weight > HIGH_BAR):
            continue
        if not any(rho_ratio(ONE, vertex) == HIGH_BAR for vertex in right):
            continue
        packet = (ONE,) + right
        normal = tuple(value / min(packet) for value in sorted(packet))
        valid_packets.add(normal)
    one_six = [(maximum_spanning_tree(packet), packet) for packet in valid_packets]
    require(one_six, "empty closed 1+6 cut census")
    two_component_counts["1+6"] = len(one_six)
    two_component_minima["1+6"] = min(one_six)

    valid_packets = set()
    for second_center in center_ratios:
        left = (ONE, second_center)
        if rho_ratio(*left) <= HIGH_BAR:
            continue
        common = tuple(vertex for vertex in closed_ratios if vertex / second_center in closed_oriented)
        for right in combinations(common, 5):
            if not graph_connected(right, lambda weight: weight > HIGH_BAR):
                continue
            if not any(rho_ratio(a, b) == HIGH_BAR for a in left for b in right):
                continue
            packet = tuple(sorted(left + right))
            normal = tuple(value / packet[0] for value in packet)
            valid_packets.add(normal)
    two_five = [(maximum_spanning_tree(packet), packet) for packet in valid_packets]
    require(two_five, "empty closed 2+5 cut census")
    two_component_counts["2+5"] = len(two_five)
    two_component_minima["2+5"] = min(two_five)

    valid_packets = set()
    for right in combinations(closed_ratios, 4):
        centers: set[F] | None = None
        for vertex in right:
            possible = {vertex / ratio for ratio in closed_oriented}
            centers = possible if centers is None else centers & possible
        require(centers is not None, "empty center initialization")
        if ONE not in centers:
            continue
        for other_centers in combinations(sorted(centers - {ONE}), 2):
            left = (ONE,) + other_centers
            if not graph_connected(left, lambda weight: weight > HIGH_BAR):
                continue
            if not graph_connected(right, lambda weight: weight > HIGH_BAR):
                continue
            if not any(rho_ratio(a, b) == HIGH_BAR for a in left for b in right):
                continue
            packet = tuple(sorted(left + right))
            normal = tuple(value / packet[0] for value in packet)
            valid_packets.add(normal)
    require(not valid_packets, "unexpected closed 3+4 strict-component cut")
    two_component_counts["3+4"] = 0

    require(two_component_counts == {"1+6": 131593, "2+5": 12, "3+4": 0}, "closed cut counts")
    require(two_component_minima["1+6"][0] == F(2228, 18711), "closed 1+6 minimum")
    require(two_component_minima["2+5"][0] == F(15473, 103950), "closed 2+5 minimum")
    require(all(value[0] > TREE_FLOOR for value in two_component_minima.values()), "two components miss strong floor")
    print("strict-high two-component global cut census:")
    print(
        f"  1+6 valid normalized packets={two_component_counts['1+6']}; "
        f"MST min={fmt(two_component_minima['1+6'][0])}"
    )
    print(
        f"  2+5 valid normalized packets={two_component_counts['2+5']}; "
        f"MST min={fmt(two_component_minima['2+5'][0])}"
    )
    print("  3+4 valid normalized packets=0")

    # The remaining connected branch has a connected strict-high graph, so
    # it has a six-edge tree and every edge is at least the first strict-high
    # spectral value 5/308.
    require(6 * NEXT_HIGH == TREE_FLOOR, "strong tree-floor identity")
    require(TREE_FLOOR - ALIGNMENT_FLOOR == F(5, 2772), "strong-vs-aligned margin")
    print("connected strict-high graph: six edges >=5/308, hence tree>=15/154")
    print("all other threshold-component branches have larger exact minima")

    require(TREE_FLOOR > HEAVY_CIRCUIT_CUTOFF, "does not cross THM-1203 cutoff")
    require(TREE_FLOOR - HEAVY_CIRCUIT_CUTOFF == F(45, 14014), "heavy-cutoff margin")
    require(TREE_FLOOR > F(110, 1183), "does not improve THM-1166 tree floor")
    require(TREE_FLOOR - F(110, 1183) == F(115, 26026), "old-floor margin")
    require(1 - TREE_FLOOR == F(139, 154), "covered-period complement")
    require(7 * (1 - TREE_FLOOR) == F(139, 22), "protected-needle constant")
    print(f"universal maximum-spanning-tree / Hunter safe floor: {fmt(TREE_FLOOR)}")
    print(f"margin over 60/637: {fmt(TREE_FLOOR-HEAVY_CIRCUIT_CUTOFF)}")
    print(f"margin over 110/1183: {fmt(TREE_FLOOR-F(110,1183))}")
    print("common-dilate covered-period bound: g|I|<=139/154")
    print("protected-needle consequence: g/max(core)<=139/22")

    # Threshold-aligned incidence budgets requested by the THM-1203/1218
    # composition.  These deliberately retain the weaker 265/2772 subfloor
    # so their denominators expose the exact old top-stratum comparison.
    require(ALIGNMENT_FLOOR - F(2, 21) == F(1, 2772), "AP incidence budget")
    require(
        ALIGNMENT_FLOOR - HEAVY_CIRCUIT_CUTOFF == F(355, 252252),
        "non-AP incidence budget",
    )
    require(TREE_FLOOR - F(2, 21) == F(1, 462), "strong AP incidence budget")
    require(
        TREE_FLOOR - HEAVY_CIRCUIT_CUTOFF == F(45, 14014),
        "strong non-AP incidence budget",
    )
    print("THM-1203 incidence budgets:")
    print("  aligned floor minus arbitrary BAD quartet: 265/2772-2/21=1/2772")
    print("  aligned floor minus non-AP BAD quartet: 265/2772-60/637=355/252252")
    print("  strong floor minus arbitrary BAD quartet: 15/154-2/21=1/462")
    print("  strong floor minus non-AP BAD quartet: 15/154-60/637=45/14014")

    # Tournament telemetry.  Scale the exact diagnostic minimizer to integer
    # speeds so the fingerprint remains an honest runner packet.
    denominator = lcm(*(ratio.denominator for ratio in (ONE,) + exact_min[1]))
    integer_packet = tuple(int(denominator * ratio) for ratio in (ONE,) + exact_min[1])
    fingerprint = tournament_fingerprint(tuple(F(value) for value in integer_packet))
    print("tournament audit (vertices=runners; low edge reversed from speed order)")
    print(f"sample integer packet: {integer_packet}")
    print(f"score histogram: {fingerprint[0]}; directed triangles: {fingerprint[1]}")
    print(f"SCC sizes: {fingerprint[2]}; strict-low edge flips: {fingerprint[3]}")
    print(f"tie Hamiltonian path: {fingerprint[4]}; Hamiltonian-path count: {fingerprint[5]}")
    print("preserved: strict-spectrum cut, high connectivity, and global Hunter credit")
    print("destroyed: interval phase, gcd positioning error, tooth ownership, and BAD quartet mass")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={source_hash}")


if __name__ == "__main__":
    main()
