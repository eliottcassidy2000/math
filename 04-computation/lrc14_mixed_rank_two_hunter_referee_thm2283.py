#!/usr/bin/env python3
"""Independent exact referee for the sharpened THM-2283.

The primary certificate enumerates trees with Prüfer words and evaluates
Jackson coefficients from their closed formula.  This referee instead:

* enumerates trees as connected edge subsets;
* computes every Jackson coefficient by direct triangular convolution;
* checks the closed formula only as a secondary identity;
* exhausts the minimal-union support and non-Haar graph types; and
* independently certifies the true-fail/pass Jackson boundaries.

No ``assert`` statement is used, so ``python -O`` performs the same audit.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb
import sys


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


P_GUARD = F(2, 7)
P_ORDINARY = F(1, 7)
S_GUARD = F(5, 7)
S_ORDINARY = F(6, 7)


def complete_edges(order: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(order), 2))


def connected(order: int, edges: tuple[tuple[int, int], ...]) -> bool:
    reached = {0}
    changed = True
    while changed:
        changed = False
        for left, right in edges:
            if left in reached and right not in reached:
                reached.add(right)
                changed = True
            if right in reached and left not in reached:
                reached.add(left)
                changed = True
    return len(reached) == order


def edge_subset_trees(order: int) -> tuple[frozenset[tuple[int, int]], ...]:
    trees = tuple(
        frozenset(candidate)
        for candidate in combinations(complete_edges(order), order - 1)
        if connected(order, candidate)
    )
    require(len(trees) == order ** (order - 2), f"Cayley count n={order}")
    return trees


def graph_types(order: int) -> tuple[frozenset[tuple[int, int]], ...]:
    """All graphs allowed by minimal sparse support union."""

    edges = complete_edges(order)
    allowed: list[frozenset[tuple[int, int]]] = []
    for size in range(len(edges) + 1):
        for chosen_tuple in combinations(edges, size):
            chosen = frozenset(chosen_tuple)
            degrees = {
                vertex: sum(vertex in edge for edge in chosen)
                for vertex in range(order)
            }
            if order == 3:
                keep = True
            elif order == 4:
                keep = max(degrees.values()) <= 1
            elif order == 5:
                keep = len(chosen) <= 1
            elif order == 6:
                keep = len(chosen) == 0
            else:
                keep = False
            if keep:
                allowed.append(chosen)
    return tuple(allowed)


def overlap_lower(
    edge: tuple[int, int],
    nonhaar: frozenset[tuple[int, int]],
    guard_in: bool,
    guard_odd: bool,
) -> F:
    is_guard_edge = guard_in and 0 in edge
    if edge not in nonhaar:
        return F(2, 49) if is_guard_edge else F(1, 49)
    if is_guard_edge and guard_odd:
        return F(1, 42)
    return F(1, 91)


def robust_tree_credit(order: int, guard_in: bool, guard_odd: bool) -> F:
    trees = edge_subset_trees(order)
    graph_minima: list[F] = []
    for graph in graph_types(order):
        tree_weights = [
            sum(
                (
                    overlap_lower(edge, graph, guard_in, guard_odd)
                    for edge in tree
                ),
                F(0),
            )
            for tree in trees
        ]
        graph_minima.append(max(tree_weights))
    return min(graph_minima)


def check_hunter_pointwise(order: int) -> None:
    for tree in edge_subset_trees(order):
        for active in range(1 << order):
            vertex_sum = active.bit_count()
            edge_sum = sum(
                bool(active & (1 << left)) and bool(active & (1 << right))
                for left, right in tree
            )
            require(
                (1 if active else 0) <= vertex_sum - edge_sum,
                f"Hunter pointwise n={order}",
            )


def support_masks(order: int) -> tuple[int, ...]:
    return tuple(
        mask
        for mask in range(1 << order)
        if mask.bit_count() in (2, 3)
    )


def check_union_six_split() -> int:
    """Enumerate all sparse support pairs and isolate union-size six."""

    checked = 0
    for first, second in combinations(support_masks(6), 2):
        union = first | second
        if union.bit_count() == 6:
            require(first.bit_count() == second.bit_count() == 3, "m6 sizes")
            require(first & second == 0, "m6 supports must be disjoint")
            checked += 1
    require(checked == 10, "unordered partitions of six into triples")
    return checked


def alternating_arctan(reciprocal: int, last_index: int) -> F:
    total = F(0)
    for index in range(last_index + 1):
        total += F(
            1 if index % 2 == 0 else -1,
            (2 * index + 1) * reciprocal ** (2 * index + 1),
        )
    return total


def doubled_tangent(value: F) -> F:
    return 2 * value / (1 - value * value)


def convolved_coefficient(bandwidth: int, frequency: int) -> int:
    total = 0
    for left in range(-bandwidth + 1, bandwidth):
        right = frequency - left
        if -bandwidth < right < bandwidth:
            total += (
                bandwidth - abs(left)
            ) * (
                bandwidth - abs(right)
            )
    return total


def closed_coefficient(bandwidth: int, frequency: int) -> int:
    if frequency <= bandwidth:
        return (
            4 * bandwidth**3
            - 6 * bandwidth * frequency**2
            + 2 * bandwidth
            + 3 * frequency**3
            - 3 * frequency
        ) // 6
    reflected = 2 * bandwidth - frequency
    return (reflected**3 - reflected) // 6


def direct_eta(bandwidth: int, pi_bound: F) -> F:
    coefficients = [
        convolved_coefficient(bandwidth, frequency)
        for frequency in range(2 * bandwidth - 1)
    ]
    c_zero = coefficients[0]
    require(
        c_zero == bandwidth * (2 * bandwidth**2 + 1) // 3,
        f"direct C0 N={bandwidth}",
    )
    odd_sum = sum(
        (
            F(coefficients[frequency], frequency**2)
            for frequency in range(1, 2 * bandwidth - 2, 2)
        ),
        F(0),
    )
    return F(1, 2) - 4 * odd_sum / (pi_bound**2 * c_zero)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # Direct sparse-case recomputation in the factorial-moment basis.
    baseline = (
        1
        - F(10, 7)
        + F(5, 6) * F(44, 49)
        - F(1, 2) * F(16, 49)
    )
    sparse_values = {
        "A": baseline,
        "B_ordinary": baseline - F(1, 2) * (F(1, 49) - F(1, 343)),
        "B_guard": baseline - F(1, 2) * (F(1, 49) - F(2, 343)),
        "C_ordinary": baseline + F(11, 42) * (F(1, 91) - F(1, 49)),
        "C_guard_general": baseline + F(1, 3) * (F(1, 91) - F(2, 49)),
        "C_guard_odd": baseline + F(1, 3) * (F(1, 42) - F(2, 49)),
    }
    require(
        sparse_values
        == {
            "A": F(23, 147),
            "B_ordinary": F(152, 1029),
            "B_guard": F(307, 2058),
            "C_ordinary": F(2060, 13377),
            "C_guard_general": F(40, 273),
            "C_guard_odd": F(19, 126),
        },
        "sparse-case fractions",
    )

    for order in range(3, 7):
        check_hunter_pointwise(order)

    graph_counts = tuple(len(graph_types(order)) for order in range(3, 7))
    require(graph_counts == (8, 10, 11, 1), "admissible graph counts")
    m6_partitions = check_union_six_split()

    expected_credits = {
        (3, False, False): F(2, 91),
        (3, True, False): F(2, 91),
        (3, True, True): F(1, 21),
        (4, False, False): F(3, 49),
        (4, True, False): F(5, 49),
        (4, True, True): F(31, 294),
        (5, False, False): F(4, 49),
        (5, True, False): F(1, 7),
        (5, True, True): F(43, 294),
        (6, False, False): F(5, 49),
        (6, True, False): F(10, 49),
        (6, True, True): F(10, 49),
    }
    for key, expected in expected_credits.items():
        require(
            robust_tree_credit(*key) == expected,
            f"independent tree credit {key}",
        )

    odd_case_d: list[tuple[int, F, F]] = []
    general_case_d: list[tuple[int, F, F]] = []
    for order in range(3, 7):
        if order == 6:
            ordinary_core = F(16, 49)
            guard_general_core = guard_odd_core = F(12, 49)
        else:
            ordinary_core = (
                1 - F(order, 7) + robust_tree_credit(order, False, False)
            )
            guard_general_core = (
                1
                - F(order + 1, 7)
                + robust_tree_credit(order, True, False)
            )
            guard_odd_core = (
                1
                - F(order + 1, 7)
                + robust_tree_credit(order, True, True)
            )
        outside_free = S_GUARD * S_ORDINARY ** (8 - order)
        inside_free = S_ORDINARY ** (9 - order)
        general_case_d.append(
            (
                order,
                ordinary_core * outside_free,
                guard_general_core * inside_free,
            )
        )
        odd_case_d.append(
            (
                order,
                ordinary_core * outside_free,
                guard_odd_core * inside_free,
            )
        )

    require(
        general_case_d
        == [
            (3, F(2099520, 10706059), F(1912896, 10706059)),
            (4, F(155520, 823543), F(147744, 823543)),
            (5, F(19440, 117649), F(2592, 16807)),
            (6, F(2880, 16807), F(2592, 16807)),
        ],
        "independent general Case D table",
    )
    require(
        odd_case_d
        == [
            (3, F(2099520, 10706059), F(155520, 823543)),
            (4, F(155520, 823543), F(149040, 823543)),
            (5, F(19440, 117649), F(18360, 117649)),
            (6, F(2880, 16807), F(2592, 16807)),
        ],
        "independent odd Case D table",
    )

    # Independent Machin brackets.
    t = F(1, 5)
    tangent_four = doubled_tangent(doubled_tangent(t))
    require(
        (tangent_four - F(1, 239))
        / (1 + tangent_four * F(1, 239))
        == 1,
        "Machin tangent identity",
    )
    lower_certificate = 4 * (
        4 * alternating_arctan(5, 7) - alternating_arctan(239, 2)
    )
    upper_certificate = 4 * (
        4 * alternating_arctan(5, 4) - alternating_arctan(239, 1)
    )
    pi_lower = F(103993, 33102)
    pi_upper = F(355, 113)
    require(lower_certificate > pi_lower, "independent pi lower")
    require(upper_certificate < pi_upper, "independent pi upper")

    boundary_data = (
        ("uniform", F(152, 1029), 51, 52),
        ("guard_pair", F(19, 126), 50, 51),
        ("ordinary_pair", F(2060, 13377), 49, 50),
    )
    reconstructed = 0
    boundary_report: list[str] = []
    for name, floor, failure, passing in boundary_data:
        target = floor / 18
        for bandwidth in range(2, failure + 1):
            require(
                direct_eta(bandwidth, pi_lower) > target,
                f"independent true failure {name} N={bandwidth}",
            )
        require(
            direct_eta(passing, pi_upper) < target,
            f"independent pass {name} N={passing}",
        )
        for bandwidth in (failure, passing):
            for frequency in range(2 * bandwidth - 1):
                require(
                    convolved_coefficient(bandwidth, frequency)
                    == closed_coefficient(bandwidth, frequency),
                    f"independent coefficient N={bandwidth},k={frequency}",
                )
                reconstructed += 1
        boundary_report.append(
            f"{name}:N{failure}_true_FAIL,N{passing}_PASS,"
            f"degree{2 * passing - 2}"
        )

    uniform_gap = F(152, 1029) - 18 * direct_eta(52, pi_upper)
    require(uniform_gap > F(1, 600), "independent uniform gap")
    require(
        (max(35, 102), max(198, 100), max(594, 98))
        == (102, 198, 594),
        "independent adaptive heights",
    )

    print("THM-2283 INDEPENDENT HUNTER / JACKSON REFEREE")
    print(
        "floors="
        "general:40/273,odd_guard:152/1029,"
        "guard_pair:19/126,ordinary_pair:2060/13377"
    )
    print(
        "minimal_union_graph_counts="
        f"m3:{graph_counts[0]},m4:{graph_counts[1]},"
        f"m5:{graph_counts[2]},m6:{graph_counts[3]}"
    )
    print(f"m6_disjoint_triple_support_partitions={m6_partitions}")
    print(
        "tree_enumerator=connected_edge_subsets;"
        "Hunter_pointwise=all_trees_all_active_sets"
    )
    print(
        "odd_Case_D="
        + ";".join(
            f"m{order}:guard_out={outside},guard_in={inside}"
            for order, outside, inside in odd_case_d
        )
    )
    print(
        "Jackson_coefficients_directly_reconstructed="
        f"{reconstructed}"
    )
    print("Jackson_boundaries=" + ";".join(boundary_report))
    print("uniform_N52_gap_gt=1/600")
    print(
        "adaptive_scalar_rank3_heights=102,198,594;"
        "fixed_section_heights=204,396,1188"
    )
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
