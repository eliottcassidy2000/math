#!/usr/bin/env python3
"""Exact referee for the universal seven-comb projective MST bound.

For D_x={t: ||xt|| <= 1/13}, put p=2/13 and

    h(x,y) = mu(D_x intersect D_y) - p^2.

If x=ga, y=gb and gcd(a,b)=1, THM-856 gives

    h(x,y) = (Q(a+b)-Q(b-a))/(169ab),
    Q(r) = (r mod 13)(13-(r mod 13)).

This script checks the complete finite ratio tables used in the proof that
MST(h) >= -237/7436 on every seven distinct positive integer speeds.  It also
audits two exact packets.  The toothpick packet is only a search candidate;
no global extremality claim is made for it.

Tournament Analysis uses the 21 pair-overlap obligations as vertices.  Weight
comparison is the pairwise observable, lexicographic endpoint order is the tie
gauge, and the resulting descending edge order is the tie Hamiltonian path.
The graphic-matroid rank derivative must remain attached: the tournament by
itself forgets incidence and numerical weight gaps and cannot evaluate MST.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd
import heapq


P = Fraction(2, 13)
THETA = Fraction(-3, 676)
MIN_EDGE = Fraction(-11, 1014)
SECOND_EDGE = Fraction(-18, 1859)
UNIVERSAL_MST_BOUND = Fraction(-237, 7436)
PROJECTIVE_MARGIN = Fraction(19, 572)
HARMONIC_7 = Fraction(363, 140)
COMMON_DILATE_COST = Fraction(1203, 140)
COMMON_DILATE_ENDPOINT_COEFFICIENT = Fraction(1203, 70)


def Q(value: int) -> int:
    residue = value % 13
    return residue * (13 - residue)


def numerator_defect(a: int, b: int) -> int:
    return Q(a + b) - Q(b - a)


def reduced_pair(x: int, y: int) -> tuple[int, int]:
    common = gcd(x, y)
    a, b = x // common, y // common
    return (a, b) if a <= b else (b, a)


def h(x: int, y: int) -> Fraction:
    a, b = reduced_pair(x, y)
    return Fraction(numerator_defect(a, b), 169 * a * b)


def reduced_pairs(product_cutoff: int):
    for a in range(1, product_cutoff + 1):
        for b in range(a, product_cutoff // a + 1):
            if gcd(a, b) == 1:
                yield a, b


def complete_bad_ratio_table():
    # h < -3/676 is D/(ab) < -3/4.  Since D >= -42, ab >= 56
    # cannot contribute, including the equality boundary ab=56.
    rows = []
    for a, b in reduced_pairs(55):
        quotient = Fraction(numerator_defect(a, b), a * b)
        if quotient < Fraction(-3, 4):
            rows.append((a, b, quotient))
    expected = [
        (1, 9, Fraction(-10, 9)),
        (1, 10, Fraction(-7, 5)),
        (1, 11, Fraction(-18, 11)),
        (1, 12, Fraction(-11, 6)),
        (1, 25, Fraction(-22, 25)),
        (2, 9, Fraction(-10, 9)),
        (2, 11, Fraction(-18, 11)),
        (3, 10, Fraction(-7, 5)),
        (3, 11, Fraction(-28, 33)),
        (4, 9, Fraction(-10, 9)),
    ]
    assert rows == expected
    return rows


def triangle_free_ratio_audit(rows):
    ratios = set()
    for a, b, _ in rows:
        ratios.add(Fraction(a, b))
        ratios.add(Fraction(b, a))
    triangles = []
    for r, s in combinations(sorted(ratios), 2):
        if r / s in ratios or s / r in ratios:
            triangles.append((r, s))
    assert len(ratios) == 20
    assert triangles == []
    return ratios


def minimum_edge_audit():
    # To beat -11/6 before the factor 1/169 one needs ab <= 22,
    # because -42/23 > -11/6.  To beat the proposed second value
    # -18/11 after removing the minimum one needs ab <= 25, because
    # -42/26 > -18/11.
    min_rows = []
    for a, b in reduced_pairs(22):
        value = Fraction(numerator_defect(a, b), 169 * a * b)
        min_rows.append((value, a, b))
    minimum = min(value for value, _, _ in min_rows)
    minimum_at = [(a, b) for value, a, b in min_rows if value == minimum]
    assert minimum == MIN_EDGE
    assert minimum_at == [(1, 12)]

    second_rows = []
    for a, b in reduced_pairs(25):
        value = Fraction(numerator_defect(a, b), 169 * a * b)
        if value != minimum:
            second_rows.append((value, a, b))
    second = min(value for value, _, _ in second_rows)
    second_at = [(a, b) for value, a, b in second_rows if value == second]
    assert second == SECOND_EDGE
    assert second_at == [(1, 11), (2, 11)]
    return minimum_at, second_at


class DSU:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.components = size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> bool:
        left, right = self.find(left), self.find(right)
        if left == right:
            return False
        self.parent[left] = right
        self.components -= 1
        return True


def complete_edges(speeds):
    return [
        (h(speeds[i], speeds[j]), i, j)
        for i, j in combinations(range(len(speeds)), 2)
    ]


def kruskal_mst(speeds):
    edges = complete_edges(speeds)
    ordered = sorted(edges, key=lambda edge: (-edge[0], edge[1], edge[2]))
    dsu = DSU(len(speeds))
    chosen = []
    total = Fraction(0)
    for weight, left, right in ordered:
        if dsu.union(left, right):
            chosen.append((weight, left, right))
            total += weight
    assert len(chosen) == len(speeds) - 1
    return total, chosen


def threshold_rank_word(speeds):
    by_weight = {}
    for weight, left, right in complete_edges(speeds):
        by_weight.setdefault(weight, []).append((left, right))
    dsu = DSU(len(speeds))
    word = []
    for weight in sorted(by_weight, reverse=True):
        before = dsu.components
        for left, right in by_weight[weight]:
            dsu.union(left, right)
        increment = before - dsu.components
        if increment:
            word.append((weight, increment))
    assert sum(increment for _, increment in word) == len(speeds) - 1
    return word


def prufer_tree_edges(sequence, size: int):
    degree = [1] * size
    for vertex in sequence:
        degree[vertex] += 1
    leaves = [vertex for vertex in range(size) if degree[vertex] == 1]
    heapq.heapify(leaves)
    edges = []
    for vertex in sequence:
        leaf = heapq.heappop(leaves)
        edges.append((min(leaf, vertex), max(leaf, vertex)))
        degree[leaf] -= 1
        degree[vertex] -= 1
        if degree[vertex] == 1:
            heapq.heappush(leaves, vertex)
    left = heapq.heappop(leaves)
    right = heapq.heappop(leaves)
    edges.append((min(left, right), max(left, right)))
    return edges


def exhaustive_tree_mst(speeds):
    size = len(speeds)
    weights = {
        (left, right): h(speeds[left], speeds[right])
        for left, right in combinations(range(size), 2)
    }
    best = None
    maximizers = 0
    for sequence in product(range(size), repeat=size - 2):
        total = sum(
            (weights[edge] for edge in prufer_tree_edges(sequence, size)),
            Fraction(0),
        )
        if best is None or total > best:
            best = total
            maximizers = 1
        elif total == best:
            maximizers += 1
    assert best is not None
    return best, maximizers


def tournament_fingerprint(speeds):
    # Pair-overlap obligations, rather than runners, are tournament vertices.
    edges = complete_edges(speeds)
    ordered = sorted(edges, key=lambda edge: (-edge[0], edge[1], edge[2]))
    size = len(ordered)
    scores = [size - 1 - rank for rank in range(size)]
    score_histogram = Counter(scores)
    assert score_histogram == Counter(range(size))
    # Both lexicographic and reverse-lexicographic tie gauges are transitive;
    # their orientations differ on exactly the within-level pairs.
    tie_blocks = Counter(weight for weight, _, _ in edges)
    tie_pairs = sum(count * (count - 1) // 2 for count in tie_blocks.values())
    return {
        "vertices": size,
        "score_histogram": score_histogram,
        "directed_3cycles": 0,
        "sccs": size,
        "hamiltonian_paths": 1,
        "tie_blocks": sorted(Counter(tie_blocks.values()).items()),
        "tie_gauge_flips": tie_pairs,
    }


def packet_audit(name, speeds, expected_mst, expected_word=None, all_negative=False):
    kruskal_total, chosen = kruskal_mst(speeds)
    prufer_total, maximizers = exhaustive_tree_mst(speeds)
    word = threshold_rank_word(speeds)
    assert kruskal_total == prufer_total == expected_mst
    if expected_word is not None:
        assert tuple(increment for _, increment in word) == expected_word
    negative_edges = sum(weight < 0 for weight, _, _ in complete_edges(speeds))
    if all_negative:
        assert negative_edges == 21
    print(f"packet={name}")
    print(f"  speeds={speeds}")
    print(f"  negative_edges={negative_edges}/21")
    print(f"  mst={kruskal_total} maximizers={maximizers}")
    print(
        "  threshold_rank_word="
        + str([(str(weight), increment) for weight, increment in word])
    )
    print(
        "  chosen="
        + str(
            [
                (speeds[left], speeds[right], str(weight))
                for weight, left, right in chosen
            ]
        )
    )
    print(f"  tournament={tournament_fingerprint(speeds)}")


def main():
    rows = complete_bad_ratio_table()
    ratios = triangle_free_ratio_audit(rows)
    minimum_at, second_at = minimum_edge_audit()

    print("SEVEN-COMB UNIVERSAL PROJECTIVE MST REFEREE")
    print("method=integer/Fraction exact; no floating point")
    print("covariance_identity=h_xy=<1_Dx-2/13,1_Dy-2/13>_L2")
    print("covariance_scope=all positive integer speeds; no common-scale hypothesis")
    print("bad_threshold=-3/676; numerator threshold=-3/4")
    print("cutoff_proof=D>=-42 and ab>=56 imply D/ab>=-3/4")
    print("bad_ratio_table:")
    for a, b, quotient in rows:
        print(f"  ({a},{b}) D/ab={quotient} h={quotient/169}")
    print(f"bad_directed_ratios={len(ratios)}")
    print("bad_ratio_triangles=0")
    print("bad_threshold_graph=triangle-free")
    print("minimum_cutoff=22 because -42/23>-11/6")
    print(f"minimum_h={MIN_EDGE} ratios={minimum_at}")
    print("second_cutoff=25 because -42/26>-18/11")
    print(f"second_distinct_h={SECOND_EDGE} ratios={second_at}")
    print("minimum_edge_graph=disjoint subgraph of x--12x; maximum_degree=2")

    connected_bound = 6 * THETA
    disconnected_bound = SECOND_EDGE + 5 * THETA
    assert connected_bound == Fraction(-9, 338)
    assert disconnected_bound == UNIVERSAL_MST_BOUND
    # For every nontrivial two-block partition of seven vertices the cross cut
    # has at least six edges and cannot lie in a graph of maximum degree two.
    for left_size in range(1, 7):
        right_size = 7 - left_size
        assert left_size * right_size >= 6
        assert max(left_size, right_size) > 2
    assert Fraction(11, 169) + UNIVERSAL_MST_BOUND == PROJECTIVE_MARGIN
    assert sum((Fraction(1, index) for index in range(1, 8)), Fraction(0)) == HARMONIC_7
    assert HARMONIC_7 + 6 == COMMON_DILATE_COST
    assert 2 * COMMON_DILATE_COST == COMMON_DILATE_ENDPOINT_COEFFICIENT
    print("component_argument=bad triangle-free => good graph has at most 2 components")
    print(f"connected_good_graph_bound={connected_bound}")
    print("two_component_cross_cut=at least 6 pairs, so not all minimum edges")
    print(f"two_component_bound={disconnected_bound}")
    print(f"universal_MST_bound={UNIVERSAL_MST_BOUND}")
    print(f"seven_comb_projective_margin={PROJECTIVE_MARGIN}")
    print(f"distinct_reduced_speed_reciprocal_sum<=H_7={HARMONIC_7}")
    print(f"common_dilate_reciprocal_cost<={COMMON_DILATE_COST}/g")
    print(
        "common_dilate_Hunter_bound="
        f"19e/572-({COMMON_DILATE_ENDPOINT_COEFFICIENT})c_E/g"
    )

    all_negative = [4, 9, 21, 32, 70, 170, 189]
    packet_audit(
        "all_negative_counterexample",
        all_negative,
        Fraction(-505, 1447992),
        all_negative=True,
    )

    toothpick = [5, 8, 22, 36, 64, 176, 288]
    base = {8, 22, 36}
    assert set(toothpick) == {5} | base | {8 * value for value in base}
    packet_audit(
        "toothpick_candidate",
        toothpick,
        Fraction(-941, 334620),
        expected_word=(1, 3, 1, 1),
    )
    assert Fraction(1, 20280) - 3 * Fraction(1, 14872) - Fraction(
        2, 5577
    ) - Fraction(7, 3042) == Fraction(-941, 334620)
    print("toothpick_shape={5} union A union 8A, A={8,22,36}")
    print("toothpick_status=SEARCH_ONLY; exact value audited, global minimum NOT claimed")

    print("TOURNAMENT ANALYSIS")
    print("vertices=21 pair-overlap obligations (alternate carrier, not runners)")
    print("pairwise_observable=sign(h_e-h_f)")
    print("tie_gauge=lexicographic endpoint order")
    print("tie_Hamiltonian_path=descending (h_e, gauge) edge order")
    print("rank_derivative=graphic rank increase at each h-level")
    print("preserves=Kruskal order only when weights, labels, and graphic incidence remain attached")
    print("destroys=pure tournament forgets numerical gaps, endpoints, cuts, and MST value")
    print("assumption_challenge=tournament vertices need not be runners; here they are proof edges")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
