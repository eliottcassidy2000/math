#!/usr/bin/env python3
"""Exact referee for THM-2793.

The script verifies:

1. one oriented ramp moment reconstructs every interval endpoint in all
   one-of-each-length families through rank eight;
2. on every graceful labeling of every nonisomorphic tree through seven
   vertices, the reconstructed ordered gap matrix agrees literally;
3. anchored mod-two cut inversion reconstructs every threshold set and the
   complete vertex labeling;
4. every marked gap-tail predicate is preserved; and
5. Gram-only order loss and root-anchor loss are both real.

All truth-bearing checks use explicit exceptions, so optimized execution
performs the same verification.
"""

from __future__ import annotations

from itertools import permutations, product
from math import factorial

import networkx as nx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def interval_mask(start_zero: int, length: int) -> int:
    return ((1 << length) - 1) << start_zero


def interval_moment(mask: int, ambient_size: int) -> int:
    """First moment against the one-based ramp (1,...,ambient_size)."""
    return sum(
        point + 1 for point in range(ambient_size) if (mask >> point) & 1
    )


def decode_interval(length: int, moment: int, ambient_size: int) -> int:
    """Recover a nonempty interval from length and one-based first moment."""
    require(length > 0, "ramp decoder requires a nonempty interval")
    numerator = 2 * moment - length * (length - 1)
    denominator = 2 * length
    require(
        numerator % denominator == 0,
        "interval ramp moment did not decode to an integral left endpoint",
    )
    left_one = numerator // denominator
    right_one = left_one + length - 1
    require(
        1 <= left_one <= right_one <= ambient_size,
        "decoded interval left the ambient line",
    )
    return interval_mask(left_one - 1, length)


def gram(rows: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple((left & right).bit_count() for right in rows)
        for left in rows
    )


def generic_interval_audit(max_rank: int = 8):
    by_rank = {}
    total = 0
    for rank in range(1, max_rank + 1):
        families = 0
        for starts in product(
            *(range(rank - length + 1) for length in range(1, rank + 1))
        ):
            rows = tuple(
                interval_mask(start, length)
                for length, start in enumerate(starts, start=1)
            )
            matrix = gram(rows)
            moments = tuple(
                interval_moment(row, rank) for row in rows
            )
            decoded = tuple(
                decode_interval(matrix[index][index], moments[index], rank)
                for index in range(rank)
            )
            require(decoded == rows, "ramp failed to recover interval order")
            families += 1
        require(
            families == factorial(rank),
            "one-of-each-length family count changed",
        )
        by_rank[rank] = families
        total += families
    return by_rank, total


def is_graceful(edges, labels) -> bool:
    edge_count = len(edges)
    return (
        sorted(labels) == list(range(edge_count + 1))
        and sorted(abs(labels[u] - labels[v]) for u, v in edges)
        == list(range(1, edge_count + 1))
    )


def graceful_gap_rows(edges, labels) -> tuple[int, ...]:
    edge_count = len(edges)
    rows = [None] * edge_count
    for u, v in edges:
        low = min(labels[u], labels[v])
        high = max(labels[u], labels[v])
        difference = high - low
        require(rows[difference - 1] is None, "edge difference repeated")
        # Bit t-1 represents the internal threshold t in {1,...,m}.
        rows[difference - 1] = sum(
            1 << (threshold - 1)
            for threshold in range(low + 1, high + 1)
        )
    require(all(row is not None for row in rows), "edge spectrum incomplete")
    return tuple(rows)


def edge_by_difference(edges, labels):
    result = {}
    for edge in edges:
        u, v = edge
        difference = abs(labels[u] - labels[v])
        require(difference not in result, "edge difference repeated")
        result[difference] = edge
    return result


def invert_cut_column(
    vertex_count: int,
    edge_rows,
    rows,
    threshold: int,
    anchor: int,
):
    """Invert one tree-cut column, choosing the side containing ``anchor``."""
    adjacency = [[] for _ in range(vertex_count)]
    for difference, (u, v) in edge_rows.items():
        crosses = bool((rows[difference - 1] >> (threshold - 1)) & 1)
        adjacency[u].append((v, crosses))
        adjacency[v].append((u, crosses))

    membership = [None] * vertex_count
    membership[anchor] = True
    order = [anchor]
    for u in order:
        for v, crosses in adjacency[u]:
            candidate = membership[u] != crosses
            if membership[v] is None:
                membership[v] = candidate
                order.append(v)
            else:
                require(
                    membership[v] == candidate,
                    "tree-cut inversion became inconsistent",
                )
    require(len(order) == vertex_count, "cut inversion missed a vertex")
    return tuple(bool(value) for value in membership)


def suffix_start(bits: tuple[bool, ...]):
    first = None
    seen = False
    for index, bit in enumerate(bits, start=1):
        if bit:
            if first is None:
                first = index
            seen = True
        elif seen:
            return None
    return first if first is not None else len(bits) + 1


def gap_lift(label: int, gap: int) -> int:
    return label if label < gap else label + 1


def gap_tail_predicate(rows, labels, parent, gap):
    edge_count = len(rows)
    if gap in (0, edge_count + 1):
        bits = (False,) * edge_count
    else:
        bits = tuple(
            bool((row >> (gap - 1)) & 1) for row in rows
        )
    missing = suffix_start(bits)
    new_difference = abs(gap - gap_lift(labels[parent], gap))
    return missing is not None and missing == new_difference


def graceful_reconstruction_audit():
    expected_counts = {2: 2, 3: 4, 4: 16, 5: 68, 6: 376, 7: 2184}
    counts = {}
    total_labelings = 0
    total_gap_tests = 0
    total_cut_tests = 0
    for vertex_count in range(2, 8):
        local = 0
        for graph in nx.generators.nonisomorphic_trees(vertex_count):
            edges = tuple(graph.edges())
            edge_count = vertex_count - 1
            for labels in permutations(range(vertex_count)):
                if not is_graceful(edges, labels):
                    continue
                local += 1
                total_labelings += 1
                rows = graceful_gap_rows(edges, labels)
                matrix = gram(rows)
                moments = tuple(
                    interval_moment(row, edge_count) for row in rows
                )
                reconstructed_rows = tuple(
                    decode_interval(
                        matrix[index][index],
                        moments[index],
                        edge_count,
                    )
                    for index in range(edge_count)
                )
                require(
                    reconstructed_rows == rows,
                    "ramp reconstruction changed a graceful gap row",
                )

                by_difference = edge_by_difference(edges, labels)
                anchor = labels.index(0)
                threshold_sets = []
                for threshold in range(1, edge_count + 1):
                    membership = invert_cut_column(
                        vertex_count,
                        by_difference,
                        reconstructed_rows,
                        threshold,
                        anchor,
                    )
                    literal = tuple(
                        label < threshold for label in labels
                    )
                    require(
                        membership == literal,
                        "anchored cut inversion chose the wrong threshold side",
                    )
                    require(
                        sum(membership) == threshold,
                        "threshold set lost its rank",
                    )
                    if threshold_sets:
                        require(
                            all(
                                not threshold_sets[-1][vertex]
                                or membership[vertex]
                                for vertex in range(vertex_count)
                            ),
                            "reconstructed threshold sets stopped nesting",
                        )
                    threshold_sets.append(membership)
                    total_cut_tests += 1

                reconstructed_labels = tuple(
                    edge_count
                    - sum(
                        membership[vertex]
                        for membership in threshold_sets
                    )
                    for vertex in range(vertex_count)
                )
                require(
                    reconstructed_labels == labels,
                    "rooted threshold sets lost the graceful labeling",
                )

                # The ramp moment is the edge-sum sidecar in closed form.
                for difference, (u, v) in by_difference.items():
                    require(
                        2 * moments[difference - 1]
                        == difference * (labels[u] + labels[v] + 1),
                        "edge-sum/ramp identity changed",
                    )

                for parent in range(vertex_count):
                    for gap in range(edge_count + 2):
                        original = gap_tail_predicate(
                            rows, labels, parent, gap
                        )
                        reconstructed = gap_tail_predicate(
                            reconstructed_rows,
                            reconstructed_labels,
                            parent,
                            gap,
                        )
                        require(
                            reconstructed == original,
                            "marked gap-tail decision changed after reconstruction",
                        )
                        total_gap_tests += 1
        counts[vertex_count] = local
    require(counts == expected_counts, "graceful-labeling census changed")
    return counts, total_labelings, total_cut_tests, total_gap_tests


def sharp_boundaries():
    # Same Gram and same column multiset, different non-reflected order.
    first = (0b001, 0b011, 0b111)
    second = (0b010, 0b011, 0b111)
    require(gram(first) == gram(second), "Gram-only hostile changed")
    first_moments = tuple(interval_moment(row, 3) for row in first)
    second_moments = tuple(interval_moment(row, 3) for row in second)
    require(
        first_moments != second_moments,
        "oriented ramp stopped separating the Gram-only hostile",
    )
    require(
        tuple(decode_interval(row.bit_count(), moment, 3)
              for row, moment in zip(first, first_moments))
        == first
        and tuple(decode_interval(row.bit_count(), moment, 3)
                  for row, moment in zip(second, second_moments))
        == second,
        "ramp hostile controls failed to reconstruct",
    )

    # One-edge tree: the same gap row has two complementary labelings.
    edge = ((0, 1),)
    labels_left = (0, 1)
    labels_right = (1, 0)
    rows_left = graceful_gap_rows(edge, labels_left)
    rows_right = graceful_gap_rows(edge, labels_right)
    require(
        rows_left == rows_right == (0b1,),
        "root-anchor hostile lost its common cut column",
    )
    left_membership = invert_cut_column(
        2, {1: (0, 1)}, rows_left, 1, anchor=0
    )
    right_membership = invert_cut_column(
        2, {1: (0, 1)}, rows_right, 1, anchor=1
    )
    require(
        left_membership == (True, False)
        and right_membership == (False, True),
        "root anchor stopped resolving the complementary cut gauge",
    )
    return first_moments, second_moments


def main() -> None:
    interval_counts, interval_total = generic_interval_audit()
    (
        graceful_counts,
        graceful_total,
        cut_tests,
        gap_tests,
    ) = graceful_reconstruction_audit()
    first_moments, second_moments = sharp_boundaries()

    print("THM2793 ORIENTED RAMP / ROOTED GAP EXACT REFEREE")
    print("interval families:", interval_counts, "total=", interval_total)
    print(
        "graceful labelings:",
        graceful_counts,
        "total=",
        graceful_total,
    )
    print("rooted threshold-cut inversions:", cut_tests)
    print("marked gap-tail decisions:", gap_tests)
    print(
        "Gram-only ordered hostile ramp moments:",
        first_moments,
        second_moments,
    )
    print("one-edge root-anchor hostile: common_gap=1 two_complements=1")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
