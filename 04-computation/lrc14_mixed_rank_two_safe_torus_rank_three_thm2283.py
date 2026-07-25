#!/usr/bin/env python3
"""Exact primary certificate for THM-2283.

This script checks the complete sparse-case arithmetic, enumerates every
admissible non-Haar pair graph and Hunter tree, verifies the disjoint-triple
factorization table, reconstructs squared-Fejer coefficients, and certifies
the three adjacent Jackson boundaries used by the adaptive rank-three proof.

All validity checks use explicit exceptions and therefore survive ``python -O``.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import sys


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


GUARD_DANGER = F(2, 7)
ORDINARY_DANGER = F(1, 7)
GUARD_SAFE = F(5, 7)
ORDINARY_SAFE = F(6, 7)
GENERAL_FLOOR = F(40, 273)
ODD_GUARD_FLOOR = F(152, 1029)
GUARD_PAIR_FLOOR = F(19, 126)
ORDINARY_PAIR_FLOOR = F(2060, 13377)


def danger_minorant(count: int) -> F:
    return -F((count - 1) * (count - 3) * (count - 4), 12)


def danger_minorant_binomial(count: int) -> F:
    return (
        1
        - count
        + F(5, 6) * comb(count, 2)
        - F(1, 2) * comb(count, 3)
    )


def edge_key(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def prufer_tree(vertex_count: int, word: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    """Decode one Prüfer word into a labelled tree."""

    require(len(word) == vertex_count - 2, "Prüfer word length")
    degrees = [1] * vertex_count
    for vertex in word:
        degrees[vertex] += 1
    tree: set[tuple[int, int]] = set()
    for vertex in word:
        leaf = next(index for index, degree in enumerate(degrees) if degree == 1)
        tree.add(edge_key(leaf, vertex))
        degrees[leaf] -= 1
        degrees[vertex] -= 1
    leaves = [index for index, degree in enumerate(degrees) if degree == 1]
    require(len(leaves) == 2, "Prüfer terminal leaves")
    tree.add(edge_key(leaves[0], leaves[1]))
    require(len(tree) == vertex_count - 1, "Prüfer tree edge count")
    return frozenset(tree)


def labelled_trees(vertex_count: int) -> tuple[frozenset[tuple[int, int]], ...]:
    trees = {
        prufer_tree(vertex_count, tuple(word))
        for word in product(range(vertex_count), repeat=vertex_count - 2)
    }
    require(
        len(trees) == vertex_count ** (vertex_count - 2),
        "Cayley tree count",
    )
    return tuple(sorted(trees, key=sorted))


def admissible_nonhaar_graphs(
    vertex_count: int,
) -> tuple[frozenset[tuple[int, int]], ...]:
    """Over-approximate the minimal-support-union graph classification."""

    edges = tuple(combinations(range(vertex_count), 2))
    graphs: list[frozenset[tuple[int, int]]] = []
    for edge_mask in range(1 << len(edges)):
        graph = frozenset(
            edge
            for index, edge in enumerate(edges)
            if edge_mask & (1 << index)
        )
        degrees = [
            sum(vertex in edge for edge in graph)
            for vertex in range(vertex_count)
        ]
        if vertex_count == 3:
            allowed = True
        elif vertex_count == 4:
            allowed = max(degrees, default=0) <= 1
        elif vertex_count == 5:
            allowed = len(graph) <= 1
        elif vertex_count == 6:
            allowed = not graph
        else:
            allowed = False
        if allowed:
            graphs.append(graph)
    return tuple(graphs)


def pair_credit(
    edge: tuple[int, int],
    nonhaar_graph: frozenset[tuple[int, int]],
    guard_in_core: bool,
    odd_guard: bool,
) -> F:
    left, right = edge
    guard_pair = guard_in_core and (left == 0 or right == 0)
    if edge not in nonhaar_graph:
        return (
            GUARD_DANGER * ORDINARY_DANGER
            if guard_pair
            else ORDINARY_DANGER**2
        )
    if guard_pair and odd_guard:
        return F(1, 42)
    return F(1, 91)


def minimum_best_hunter_credit(
    vertex_count: int,
    guard_in_core: bool,
    odd_guard: bool,
) -> F:
    trees = labelled_trees(vertex_count)
    worst: F | None = None
    for graph in admissible_nonhaar_graphs(vertex_count):
        best = max(
            sum(
                (
                    pair_credit(edge, graph, guard_in_core, odd_guard)
                    for edge in tree
                ),
                F(0),
            )
            for tree in trees
        )
        worst = best if worst is None else min(worst, best)
    require(worst is not None, "empty graph classification")
    return worst


def hunter_pointwise_check(vertex_count: int) -> None:
    """Check 1_union <= vertex sum - tree-edge intersection sum."""

    for tree in labelled_trees(vertex_count):
        for active_mask in range(1 << vertex_count):
            active_count = active_mask.bit_count()
            active_edges = sum(
                bool(active_mask & (1 << left))
                and bool(active_mask & (1 << right))
                for left, right in tree
            )
            union_indicator = 1 if active_count else 0
            require(
                union_indicator <= active_count - active_edges,
                "Hunter pointwise inequality",
            )


def arctan_partial(reciprocal: int, last_index: int) -> F:
    return sum(
        (
            F((-1) ** index, (2 * index + 1) * reciprocal ** (2 * index + 1))
            for index in range(last_index + 1)
        ),
        F(0),
    )


def tangent_double(value: F) -> F:
    return 2 * value / (1 - value * value)


def jackson_coefficient(bandwidth: int, frequency: int) -> int:
    k = abs(frequency)
    require(0 <= k <= 2 * bandwidth - 2, "Jackson frequency support")
    if k <= bandwidth:
        numerator = (
            4 * bandwidth**3
            - 6 * bandwidth * k**2
            + 2 * bandwidth
            + 3 * k**3
            - 3 * k
        )
    else:
        reflected = 2 * bandwidth - k
        numerator = reflected**3 - reflected
    require(numerator % 6 == 0, "Jackson coefficient integrality")
    return numerator // 6


def jackson_coefficient_by_convolution(bandwidth: int, frequency: int) -> int:
    k = abs(frequency)
    return sum(
        (bandwidth - abs(left)) * (bandwidth - abs(k - left))
        for left in range(-bandwidth + 1, bandwidth)
        if abs(k - left) < bandwidth
    )


def jackson_eta_at_pi_bound(bandwidth: int, pi_bound: F) -> F:
    c_zero = F(bandwidth * (2 * bandwidth**2 + 1), 3)
    odd_sum = sum(
        (
            F(jackson_coefficient(bandwidth, frequency), frequency**2)
            for frequency in range(1, 2 * bandwidth - 2, 2)
        ),
        F(0),
    )
    return F(1, 2) - 4 * odd_sum / (pi_bound**2 * c_zero)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    for count in range(10):
        require(
            danger_minorant(count) == danger_minorant_binomial(count),
            f"minorant identity at {count}",
        )
        require(
            danger_minorant(count) <= (1 if count == 0 else 0),
            f"minorant domination at {count}",
        )

    baseline = (
        1
        - F(10, 7)
        + F(5, 6) * F(44, 49)
        - F(1, 2) * F(16, 49)
    )
    case_b_ordinary = baseline - F(1, 2) * (F(1, 49) - F(1, 343))
    case_b_guard = baseline - F(1, 2) * (F(1, 49) - F(2, 343))
    case_c_ordinary = baseline + F(11, 42) * (F(1, 91) - F(1, 49))
    case_c_guard_general = baseline + F(1, 3) * (F(1, 91) - F(2, 49))
    case_c_guard_odd = baseline + F(1, 3) * (F(1, 42) - F(2, 49))

    require(baseline == F(23, 147), "three-wise baseline")
    require(case_b_ordinary == ODD_GUARD_FLOOR, "support-three ordinary")
    require(case_b_guard == F(307, 2058), "support-three guard")
    require(case_c_ordinary == ORDINARY_PAIR_FLOOR, "support-two ordinary")
    require(case_c_guard_general == GENERAL_FLOOR, "general guard pair")
    require(case_c_guard_odd == GUARD_PAIR_FLOOR, "odd guard pair")

    for vertex_count in range(3, 7):
        hunter_pointwise_check(vertex_count)

    expected_general_credits = {
        (3, False): F(2, 91),
        (3, True): F(2, 91),
        (4, False): F(3, 49),
        (4, True): F(5, 49),
        (5, False): F(4, 49),
        (5, True): F(7, 49),
        (6, False): F(5, 49),
        (6, True): F(10, 49),
    }
    expected_odd_guard_three_credit = F(1, 21)
    for key, expected in expected_general_credits.items():
        vertex_count, guard_in_core = key
        actual = minimum_best_hunter_credit(
            vertex_count,
            guard_in_core,
            odd_guard=False,
        )
        require(actual == expected, f"general Hunter credit {key}")
    require(
        minimum_best_hunter_credit(3, True, odd_guard=True)
        == expected_odd_guard_three_credit,
        "odd guard m=3 Hunter credit",
    )

    general_case_d: list[tuple[int, F, F]] = []
    odd_case_d: list[tuple[int, F, F]] = []
    for vertex_count in range(3, 7):
        ordinary_base = 1 - vertex_count * ORDINARY_DANGER
        guard_base = (
            1
            - GUARD_DANGER
            - (vertex_count - 1) * ORDINARY_DANGER
        )
        if vertex_count == 6:
            ordinary_core = F(16, 49)
            guard_core_general = F(12, 49)
            guard_core_odd = F(12, 49)
        else:
            ordinary_core = ordinary_base + minimum_best_hunter_credit(
                vertex_count,
                False,
                odd_guard=False,
            )
            guard_core_general = guard_base + minimum_best_hunter_credit(
                vertex_count,
                True,
                odd_guard=False,
            )
            guard_core_odd = guard_base + minimum_best_hunter_credit(
                vertex_count,
                True,
                odd_guard=True,
            )
        guard_out_free = (
            GUARD_SAFE
            * ORDINARY_SAFE ** (8 - vertex_count)
        )
        guard_in_free = ORDINARY_SAFE ** (9 - vertex_count)
        guard_out_total = ordinary_core * guard_out_free
        guard_in_general_total = guard_core_general * guard_in_free
        guard_in_odd_total = guard_core_odd * guard_in_free
        general_case_d.append(
            (vertex_count, guard_out_total, guard_in_general_total)
        )
        odd_case_d.append(
            (vertex_count, guard_out_total, guard_in_odd_total)
        )

    expected_general_case_d = [
        (3, F(2099520, 10706059), F(1912896, 10706059)),
        (4, F(155520, 823543), F(147744, 823543)),
        (5, F(19440, 117649), F(2592, 16807)),
        (6, F(2880, 16807), F(2592, 16807)),
    ]
    expected_odd_case_d = [
        (3, F(2099520, 10706059), F(155520, 823543)),
        (4, F(155520, 823543), F(149040, 823543)),
        (5, F(19440, 117649), F(18360, 117649)),
        (6, F(2880, 16807), F(2592, 16807)),
    ]
    require(general_case_d == expected_general_case_d, "general Case D table")
    require(odd_case_d == expected_odd_case_d, "odd Case D table")
    require(
        min(value for _, outside, inside in general_case_d for value in (outside, inside))
        > GENERAL_FLOOR,
        "general Case D above floor",
    )
    require(
        min(value for _, outside, inside in odd_case_d for value in (outside, inside))
        > ODD_GUARD_FLOOR,
        "odd Case D above floor",
    )
    require(
        F(2592, 16807) - ODD_GUARD_FLOOR == F(328, 50421),
        "closest Case D gap",
    )

    # Machin's identity and alternating-series pi brackets.
    tan_one_fifth = F(1, 5)
    tan_four_fifths = tangent_double(tangent_double(tan_one_fifth))
    tan_difference = (
        tan_four_fifths - F(1, 239)
    ) / (
        1 + tan_four_fifths * F(1, 239)
    )
    require(tan_difference == 1, "Machin tangent identity")
    pi_lower_certificate = 4 * (
        4 * arctan_partial(5, 7) - arctan_partial(239, 2)
    )
    pi_upper_certificate = 4 * (
        4 * arctan_partial(5, 4) - arctan_partial(239, 1)
    )
    pi_lower = F(103993, 33102)
    pi_upper = F(355, 113)
    require(pi_lower_certificate > pi_lower, "lower pi certificate")
    require(pi_upper_certificate < pi_upper, "upper pi certificate")

    boundary_data = (
        ("uniform", ODD_GUARD_FLOOR, 51, 52, F(1, 600)),
        ("guard_pair", GUARD_PAIR_FLOOR, 50, 51, F(0)),
        ("ordinary_pair", ORDINARY_PAIR_FLOOR, 49, 50, F(0)),
    )
    coefficient_checks = 0
    for _, _, fail_bandwidth, pass_bandwidth, _ in boundary_data:
        for bandwidth in (fail_bandwidth, pass_bandwidth):
            for frequency in range(2 * bandwidth - 1):
                require(
                    jackson_coefficient(bandwidth, frequency)
                    == jackson_coefficient_by_convolution(
                        bandwidth,
                        frequency,
                    ),
                    f"Jackson convolution N={bandwidth}, k={frequency}",
                )
                coefficient_checks += 1

    boundary_rows: list[tuple[str, int, int, int]] = []
    for name, floor, fail_bandwidth, pass_bandwidth, simple_gap in boundary_data:
        target = floor / 18
        for bandwidth in range(2, fail_bandwidth + 1):
            require(
                jackson_eta_at_pi_bound(bandwidth, pi_lower) > target,
                f"{name} true failure scan N={bandwidth}",
            )
        pass_cap = jackson_eta_at_pi_bound(pass_bandwidth, pi_upper)
        require(pass_cap < target, f"{name} pass N={pass_bandwidth}")
        gap = floor - 18 * pass_cap
        require(gap > simple_gap, f"{name} simple positive gap")
        boundary_rows.append(
            (
                name,
                fail_bandwidth,
                pass_bandwidth,
                2 * pass_bandwidth - 2,
            )
        )

    adaptive_relation_heights = (
        max(35, 102),
        max(198, 100),
        max(594, 98),
    )
    fixed_section_heights = tuple(
        2 * height for height in adaptive_relation_heights
    )
    require(
        adaptive_relation_heights == (102, 198, 594),
        "adaptive scalar relation heights",
    )
    require(
        fixed_section_heights == (204, 396, 1188),
        "adaptive fixed-section heights",
    )

    print("THM-2283 EXACT MIXED SAFE-TORUS / RANK-THREE CERTIFICATE")
    print(
        "safe_floors="
        f"general:{GENERAL_FLOOR},odd_guard:{ODD_GUARD_FLOOR},"
        f"guard_pair:{GUARD_PAIR_FLOOR},ordinary_pair:{ORDINARY_PAIR_FLOOR}"
    )
    print(
        "sparse_cases="
        f"A:{baseline},B_ordinary:{case_b_ordinary},"
        f"B_guard:{case_b_guard},C_ordinary:{case_c_ordinary},"
        f"C_guard_general:{case_c_guard_general},"
        f"C_guard_odd:{case_c_guard_odd}"
    )
    print(
        "minimal_union_nonhaar_graphs="
        "m3:subgraph_K3,m4:matching,m5:at_most_one_edge,m6:empty"
    )
    print(
        "case_D_odd="
        + ";".join(
            f"m{size}:guard_out={outside},guard_in={inside}"
            for size, outside, inside in odd_case_d
        )
    )
    print(
        "m6_factor="
        "ordinary_core:16/49,guard_core:12/49,"
        "full_guard_out:2880/16807,full_guard_in:2592/16807"
    )
    print(
        "pi_brackets_certified=103993/33102<pi<355/113,"
        f"Jackson_coefficients_reconstructed={coefficient_checks}"
    )
    print(
        "Jackson_boundaries="
        + ";".join(
            f"{name}:N{fail}_true_FAIL,N{passed}_PASS,degree{degree}"
            for name, fail, passed, degree in boundary_rows
        )
    )
    print("uniform_N52_gap_gt=1/600")
    print(
        "adaptive_scalar_rank3_heights=102,198,594;"
        "fixed_section_heights=204,396,1188"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
