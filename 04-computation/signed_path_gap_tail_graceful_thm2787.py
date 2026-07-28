#!/usr/bin/env python3
"""Exact referee for THM-2787.

The script checks four logically separate layers.

1.  On every labelled tree through five vertices and every signed
    permutation of its edge magnitudes, the signed-path, bounded-injective,
    consecutive-potential, and graceful conditions agree exactly.
2.  On every nonisomorphic tree through seven vertices, every graceful
    labeling, every label gap, and every marked parent, the Ferrers/suffix
    column criterion is equivalent to literal graceful leaf insertion.
3.  Every marked vertex of every nonisomorphic tree through nine vertices
    has at least one graceful labeling and one valid gap-tail insertion.
4.  The three sharp boundary examples in the theorem are checked directly.

All truth-bearing checks use explicit exceptions, so ``python -O`` executes
the same verification.
"""

from __future__ import annotations

import hashlib
from itertools import permutations, product

import networkx as nx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def prufer_trees(n: int):
    """Yield every labelled tree on range(n), once, by Prüfer words."""
    if n == 2:
        yield ((0, 1),)
        return
    for word in product(range(n), repeat=n - 2):
        degree = [1] * n
        for value in word:
            degree[value] += 1
        edges = []
        for value in word:
            leaf = next(i for i, d in enumerate(degree) if d == 1)
            edges.append((leaf, value))
            degree[leaf] -= 1
            degree[value] -= 1
        last = [i for i, d in enumerate(degree) if d == 1]
        require(len(last) == 2, "Prüfer decoder lost its last edge")
        edges.append((last[0], last[1]))
        yield tuple(edges)


def root_potentials(n, edges, weights, signs, root=0):
    """Integrate signed edge values away from ``root``."""
    adjacency = [[] for _ in range(n)]
    for edge_index, (u, v) in enumerate(edges):
        adjacency[u].append((v, edge_index))
        adjacency[v].append((u, edge_index))
    parent = [-1] * n
    parent_edge = [-1] * n
    parent[root] = root
    order = [root]
    for u in order:
        for v, edge_index in adjacency[u]:
            if parent[v] < 0:
                parent[v] = u
                parent_edge[v] = edge_index
                order.append(v)
    require(len(order) == n, "tree integration did not reach every vertex")
    values = [0] * n
    for v in order[1:]:
        edge_index = parent_edge[v]
        values[v] = (
            values[parent[v]] + signs[edge_index] * weights[edge_index]
        )
    return tuple(values)


def is_graceful(edges, labels):
    m = len(edges)
    if sorted(labels) != list(range(m + 1)):
        return False
    differences = sorted(abs(labels[u] - labels[v]) for u, v in edges)
    return differences == list(range(1, m + 1))


def signed_path_equivalence_audit():
    assignment_count = 0
    by_n = {}
    for n in range(2, 6):
        local_count = 0
        m = n - 1
        for edges in prufer_trees(n):
            for weights in permutations(range(1, m + 1)):
                for signs in product((-1, 1), repeat=m):
                    values = root_potentials(n, edges, weights, signs)
                    pair_differences = [
                        abs(values[u] - values[v])
                        for u in range(n)
                        for v in range(u + 1, n)
                    ]
                    path_condition = all(
                        1 <= difference <= m
                        for difference in pair_differences
                    )
                    injective_bounded = (
                        len(set(values)) == n
                        and max(values) - min(values) <= m
                    )
                    consecutive = sorted(values) == list(
                        range(min(values), min(values) + n)
                    )
                    translated = tuple(
                        value - min(values) for value in values
                    )
                    graceful = is_graceful(edges, translated)
                    require(
                        path_condition
                        == injective_bounded
                        == consecutive
                        == graceful,
                        "signed-path/graceful equivalence failed",
                    )
                    assignment_count += 1
                    local_count += 1
        by_n[n] = local_count
    require(
        (assignment_count, by_n)
        == (48794, {2: 2, 3: 24, 4: 768, 5: 48000}),
        "signed-path assignment census changed",
    )
    return assignment_count, by_n


def gap_lift(label, gap):
    return label if label < gap else label + 1


def difference_edges(edges, labels):
    by_difference = {}
    for u, v in edges:
        difference = abs(labels[u] - labels[v])
        require(
            difference not in by_difference,
            "graceful edge difference repeated",
        )
        by_difference[difference] = (u, v)
    return by_difference


def gap_column(edges, labels, gap):
    """Return C_f(d,gap), in increasing old edge-difference order."""
    by_difference = difference_edges(edges, labels)
    m = len(edges)
    require(
        set(by_difference) == set(range(1, m + 1)),
        "gap column requested for a nongraceful labeling",
    )
    bits = []
    for difference in range(1, m + 1):
        u, v = by_difference[difference]
        bits.append(
            min(labels[u], labels[v])
            < gap
            <= max(labels[u], labels[v])
        )
    return tuple(bits)


def suffix_start(bits):
    """Return k for 0^(k-1)1^(m-k+1), with all-zero k=m+1."""
    first_one = None
    seen_one = False
    for index, bit in enumerate(bits, start=1):
        if bit:
            if first_one is None:
                first_one = index
            seen_one = True
        elif seen_one:
            return None
    return first_one if first_one is not None else len(bits) + 1


def literal_leaf_extension(edges, labels, parent, gap):
    lifted = tuple(gap_lift(label, gap) for label in labels)
    new_label = gap
    new_labels = lifted + (new_label,)
    new_vertex = len(labels)
    new_edges = tuple(edges) + ((parent, new_vertex),)
    return is_graceful(new_edges, new_labels), new_labels, new_edges


def predicted_leaf_extension(edges, labels, parent, gap):
    bits = gap_column(edges, labels, gap)
    missing = suffix_start(bits)
    new_difference = abs(gap - gap_lift(labels[parent], gap))
    return missing is not None and missing == new_difference


def valid_gap_witnesses(edges, labels, parent):
    m = len(edges)
    witnesses = []
    for gap in range(m + 2):
        if predicted_leaf_extension(edges, labels, parent, gap):
            bits = gap_column(edges, labels, gap)
            witnesses.append(
                (
                    gap,
                    suffix_start(bits),
                    "".join("1" if bit else "0" for bit in bits),
                )
            )
    return tuple(witnesses)


def gap_iff_audit():
    graceful_labeling_count = 0
    insertion_case_count = 0
    cut_case_count = 0
    row_case_count = 0
    by_n = {}
    for n in range(2, 8):
        local_graceful = 0
        for graph in nx.generators.nonisomorphic_trees(n):
            edges = tuple(graph.edges())
            m = n - 1
            for labels in permutations(range(n)):
                if not is_graceful(edges, labels):
                    continue
                graceful_labeling_count += 1
                local_graceful += 1
                by_difference = difference_edges(edges, labels)

                for difference in range(1, m + 1):
                    u, v = by_difference[difference]
                    row = [
                        min(labels[u], labels[v])
                        < gap
                        <= max(labels[u], labels[v])
                        for gap in range(m + 2)
                    ]
                    require(
                        sum(row) == difference,
                        "consecutive-ones row lost its edge-length sum",
                    )
                    occupied = [i for i, bit in enumerate(row) if bit]
                    require(
                        occupied
                        == list(
                            range(
                                min(labels[u], labels[v]) + 1,
                                max(labels[u], labels[v]) + 1,
                            )
                        ),
                        "edge-gap row lost consecutive-ones support",
                    )
                    row_case_count += 1

                for gap in range(m + 2):
                    low = {
                        vertex
                        for vertex, label in enumerate(labels)
                        if label < gap
                    }
                    cut_bits = tuple(
                        (u in low) != (v in low)
                        for _, (u, v) in sorted(by_difference.items())
                    )
                    require(
                        cut_bits == gap_column(edges, labels, gap),
                        "gap column no longer equals the threshold cut",
                    )
                    cut_case_count += 1
                    for parent in range(n):
                        predicted = predicted_leaf_extension(
                            edges, labels, parent, gap
                        )
                        literal, _, _ = literal_leaf_extension(
                            edges, labels, parent, gap
                        )
                        require(
                            predicted == literal,
                            "gap-tail criterion disagrees with literal insertion",
                        )
                        insertion_case_count += 1
        by_n[n] = local_graceful
    require(
        by_n
        == {
            2: 2,
            3: 4,
            4: 16,
            5: 68,
            6: 376,
            7: 2184,
        },
        f"nonisomorphic graceful-labeling census changed: {by_n}",
    )
    return (
        graceful_labeling_count,
        insertion_case_count,
        cut_case_count,
        row_case_count,
        by_n,
    )


def marked_vertex_census():
    expected_shapes = {
        2: 1,
        3: 1,
        4: 2,
        5: 3,
        6: 6,
        7: 11,
        8: 23,
        9: 47,
    }
    summary = {}
    witness_rows = []
    for n in range(2, 10):
        graphs = list(nx.generators.nonisomorphic_trees(n))
        require(
            len(graphs) == expected_shapes[n],
            "nonisomorphic-tree count changed",
        )
        marked_count = 0
        graceful_seen = 0
        for shape_index, graph in enumerate(graphs):
            edges = tuple(graph.edges())
            witnesses = [None] * n
            for labels in permutations(range(n)):
                if not is_graceful(edges, labels):
                    continue
                graceful_seen += 1
                for parent in range(n):
                    if witnesses[parent] is not None:
                        continue
                    candidates = valid_gap_witnesses(edges, labels, parent)
                    if candidates:
                        gap, missing, bits = candidates[0]
                        witnesses[parent] = (
                            tuple(labels),
                            gap,
                            missing,
                            bits,
                        )
                if all(witness is not None for witness in witnesses):
                    break
            require(
                all(witness is not None for witness in witnesses),
                f"marked gap-tail census failed at n={n}, shape={shape_index}",
            )
            for parent, witness in enumerate(witnesses):
                marked_count += 1
                witness_rows.append(
                    (n, shape_index, tuple(sorted(edges)), parent, witness)
                )
        summary[n] = (
            len(graphs),
            marked_count,
            graceful_seen,
        )
    expected_summary = {
        2: (1, 2, 1),
        3: (1, 3, 1),
        4: (2, 8, 3),
        5: (3, 15, 7),
        6: (6, 36, 20),
        7: (11, 77, 98),
        8: (23, 184, 359),
        9: (47, 423, 2403),
    }
    require(summary == expected_summary, "marked-vertex summary changed")
    digest = hashlib.sha256(repr(witness_rows).encode("ascii")).hexdigest()
    require(
        digest
        == "a838da45aa50a01ade774984c619bbd51c556d1ac8eacc8d5191d06aee26718d",
        f"marked witness bank digest changed: {digest}",
    )
    return summary, digest


def boundary_controls():
    # Nonzero path sums without the upper cap: P3, +1,+2.
    path3 = ((0, 1), (1, 2))
    values = root_potentials(3, path3, (1, 2), (1, 1))
    require(values == (0, 1, 3), "P3 upper-cap hostile changed")
    require(len(set(values)) == 3, "P3 hostile unexpectedly collided")
    require(max(values) - min(values) == 3, "P3 hostile span changed")

    # The upper cap without nonvanishing: P4, +1,-3,+2.
    path4 = ((0, 1), (1, 2), (2, 3))
    values = root_potentials(4, path4, (1, 3, 2), (1, -1, 1))
    require(values == (0, 1, -2, 0), "P4 zero-path hostile changed")
    pair_differences = [
        abs(values[u] - values[v])
        for u in range(4)
        for v in range(u + 1, 4)
    ]
    require(max(pair_differences) == 3, "P4 hostile exceeded its cap")
    require(0 in pair_differences, "P4 hostile lost its zero path")

    # Root paths alone miss a large leaf-to-leaf path.
    rooted_star = ((0, 1), (0, 2))
    values = root_potentials(3, rooted_star, (1, 2), (1, -1))
    require(values == (0, 1, -2), "root-only hostile changed")
    require(
        all(0 < abs(values[v]) <= 2 for v in (1, 2)),
        "root-only hostile lost its bounded root paths",
    )
    require(
        abs(values[1] - values[2]) == 3,
        "root-only hostile lost its nonroot violation",
    )

    # A fixed graceful K1,3 labeling need not extend at a fixed marked leaf.
    star4 = ((0, 1), (0, 2), (0, 3))
    star_labels = (0, 1, 2, 3)
    require(is_graceful(star4, star_labels), "K1,3 control is nongraceful")
    require(
        valid_gap_witnesses(star4, star_labels, 2) == (),
        "fixed-label hostile unexpectedly gap-extended",
    )

    # First obstruction to the extreme-parent (zero-rotatable) shortcut.
    tree6 = ((1, 0), (1, 2), (1, 3), (0, 4), (4, 5))
    graceful6 = [
        labels
        for labels in permutations(range(6))
        if is_graceful(tree6, labels)
    ]
    require(len(graceful6) == 12, "six-vertex hostile census changed")
    require(
        all(labels[0] not in (0, 5) for labels in graceful6),
        "marked six-vertex parent became extreme",
    )
    labels6 = (1, 4, 2, 3, 5, 0)
    witnesses6 = valid_gap_witnesses(tree6, labels6, 0)
    require(
        (3, 2, "01111") in witnesses6,
        "six-vertex gap-tail repair disappeared",
    )
    literal, extension_labels, extension_edges = literal_leaf_extension(
        tree6, labels6, 0, 3
    )
    require(literal, "six-vertex repaired insertion is nongraceful")
    require(
        extension_labels == (1, 5, 2, 4, 6, 0, 3),
        "six-vertex repaired labels changed",
    )
    require(
        sorted(
            abs(extension_labels[u] - extension_labels[v])
            for u, v in extension_edges
        )
        == list(range(1, 7)),
        "six-vertex repaired edge spectrum changed",
    )

    # Naively deleting/compressing a leaf from a graceful P5 can fail.
    path5 = ((1, 0), (1, 2), (0, 3), (3, 4))
    labels5 = (0, 2, 3, 4, 1)
    require(is_graceful(path5, labels5), "P5 descent hostile is nongraceful")
    deleted_leaf = 2
    deleted_label = labels5[deleted_leaf]
    remaining_edges = tuple(
        edge for edge in path5 if deleted_leaf not in edge
    )
    compressed = {
        vertex: (
            labels5[vertex]
            if labels5[vertex] < deleted_label
            else labels5[vertex] - 1
        )
        for vertex in range(5)
        if vertex != deleted_leaf
    }
    compressed_differences = sorted(
        abs(compressed[u] - compressed[v]) for u, v in remaining_edges
    )
    require(
        compressed_differences == [2, 2, 3],
        "P5 naive descent hostile changed",
    )
    return {
        "p3": values,
        "p4": (0, 1, -2, 0),
        "root_only": (0, 1, -2),
        "fixed_star": (star_labels, 2),
        "tree6_graceful": len(graceful6),
        "tree6_witness": (labels6, 3, 2, extension_labels),
        "descent": compressed_differences,
    }


def main():
    assignment_count, assignment_by_n = signed_path_equivalence_audit()
    (
        graceful_labeling_count,
        insertion_case_count,
        cut_case_count,
        row_case_count,
        gap_by_n,
    ) = gap_iff_audit()
    summary, witness_digest = marked_vertex_census()
    controls = boundary_controls()

    print("THM2787 SIGNED PATH / GAP-TAIL EXACT REFEREE")
    print("signed assignments:", assignment_count, assignment_by_n)
    print(
        "gap iff:",
        graceful_labeling_count,
        "graceful labelings;",
        insertion_case_count,
        "marked gap cases",
    )
    print(
        "matrix checks:",
        row_case_count,
        "consecutive-one rows;",
        cut_case_count,
        "threshold columns",
    )
    print("gap graceful counts:", gap_by_n)
    for n in sorted(summary):
        shapes, marked, graceful_seen = summary[n]
        print(
            f"marked n={n}: shapes={shapes} marked={marked} "
            f"search_graceful={graceful_seen} all_extend=1"
        )
    print("marked witness sha256:", witness_digest)
    print("P3 nonzero/no-cap:", (0, 1, 3), "span=3>2")
    print("P4 cap/zero-path:", controls["p4"], "span=3 zero=1")
    print("root-only hostile:", controls["root_only"], "leaf-gap=3")
    print("fixed K1,3 marked-leaf gap witnesses: 0")
    print(
        "six-vertex extreme hostile:",
        f"graceful={controls['tree6_graceful']} marked_extreme=0",
    )
    print("six-vertex gap repair:", controls["tree6_witness"])
    print("naive P5 leaf-compression differences:", controls["descent"])
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
