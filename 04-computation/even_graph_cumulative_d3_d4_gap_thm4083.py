#!/usr/bin/env python3
"""Proof-local exact audit for THM-4083.

This path uses no Fourier transform and no relabeling quotient.  It gauges a
signed complete graph at vertex zero, directly counts negative simple cycles
for every nontrivial character through n=7, and audits the balanced-deletion
classification used in the proof.  The n=8 row is then certified by the two
deletion-average bridges in the theorem.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import comb, factorial


N_DIRECT_MAX = 7
EXPECTED_SEMANTIC_SHA256 = "62a6f0c37659e26d36ea7dbeb98d59b355f30c4b9eddc86ea42b3fad198c91c7"
EXPECTED_B_COUNTS = {
    4: {"b0": 1, "b1": 0, "b_ge_2": 6},
    5: {"b0": 38, "b1": 15, "b_ge_2": 10},
    6: {"b0": 948, "b1": 60, "b_ge_2": 15},
    7: {"b0": 32571, "b1": 175, "b_ge_2": 21},
}
EXPECTED_D3 = {4: (4, 7), 5: (9, 10), 6: (16, 15), 7: (25, 21), 8: (36, 28)}
EXPECTED_D4 = {5: (15, 10), 6: (40, 15), 7: (85, 21), 8: (156, 28)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gauge_system(n: int) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    edges = list(itertools.combinations(range(1, n), 2))
    return edges, {edge: index for index, edge in enumerate(edges)}


def gauge_edge_mask(a: int, b: int, edge_index: dict[tuple[int, int], int]) -> int:
    if a == 0 or b == 0:
        return 0
    edge = (a, b) if a < b else (b, a)
    return 1 << edge_index[edge]


def cycle_mask(order: tuple[int, ...], edge_index: dict[tuple[int, int], int]) -> int:
    mask = 0
    for index, a in enumerate(order):
        b = order[(index + 1) % len(order)]
        mask ^= gauge_edge_mask(a, b, edge_index)
    return mask


def simple_cycles(n: int, length: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    cycles: list[int] = []
    for vertices in itertools.combinations(range(n), length):
        first = vertices[0]
        for tail in itertools.permutations(vertices[1:]):
            if tail[0] > tail[-1]:
                continue
            cycles.append(cycle_mask((first,) + tail, edge_index))
    expected = factorial(n) // (2 * length * factorial(n - length))
    require(len(cycles) == expected, f"cycle count failed n={n},k={length}")
    require(len(set(cycles)) == expected, f"cycle masks collided n={n},k={length}")
    return cycles


def triangle_data(n: int, edge_index: dict[tuple[int, int], int]) -> list[tuple[tuple[int, int, int], int]]:
    return [
        (vertices, cycle_mask(vertices, edge_index))
        for vertices in itertools.combinations(range(n), 3)
    ]


def single_edge_masks(n: int, gauge_edges: list[tuple[int, int]], edge_index: dict[tuple[int, int], int]) -> set[int]:
    masks = {1 << index for index in range(len(gauge_edges))}
    for vertex in range(1, n):
        mask = 0
        for other in range(1, n):
            if other == vertex:
                continue
            edge = (vertex, other) if vertex < other else (other, vertex)
            mask |= 1 << edge_index[edge]
        masks.add(mask)
    require(len(masks) == comb(n, 2), f"single-edge mask count failed n={n}")
    return masks


def negative_triangle_support(
    signing: int,
    triangles: list[tuple[tuple[int, int, int], int]],
) -> set[tuple[int, int, int]]:
    return {
        vertices
        for vertices, mask in triangles
        if (signing & mask).bit_count() & 1
    }


def balanced_deletions(negative_triangles: set[tuple[int, int, int]]) -> set[int]:
    require(bool(negative_triangles), "nontrivial signing had no negative triangle")
    common = set(next(iter(negative_triangles)))
    for triangle in negative_triangles:
        common.intersection_update(triangle)
    return common


def direct_negative_counts(signing: int, cycles_by_length: dict[int, list[int]]) -> dict[int, int]:
    return {
        length: sum((signing & mask).bit_count() & 1 for mask in cycles)
        for length, cycles in cycles_by_length.items()
    }


def audit_balanced_deletion_classification(
    n: int,
    signing: int,
    negative_triangles: set[tuple[int, int, int]],
    deletions: set[int],
    negative_counts: dict[int, int],
    triangles: list[tuple[tuple[int, int, int], int]],
) -> str:
    b_count = len(deletions)
    if b_count >= 2:
        u, v = sorted(deletions)[:2]
        expected_triangles = {
            tuple(sorted((u, v, x)))
            for x in range(n)
            if x not in (u, v)
        }
        require(negative_triangles == expected_triangles, f"b>=2 triangle star failed n={n},mask={signing}")
        for length, count in negative_counts.items():
            expected = factorial(n - 2) // factorial(n - length)
            require(count == expected, f"b>=2 cycle formula failed n={n},k={length},mask={signing}")
        return "b_ge_2"

    if b_count == 1:
        vertex = next(iter(deletions))
        remaining = [x for x in range(n) if x != vertex]
        root = remaining[0]
        triangle_lookup = {vertices: mask for vertices, mask in triangles}

        def is_negative_triangle(a: int, b: int, c: int) -> int:
            vertices = tuple(sorted((a, b, c)))
            return (signing & triangle_lookup[vertices]).bit_count() & 1

        potential = {root: 0}
        for x in remaining[1:]:
            potential[x] = is_negative_triangle(vertex, root, x)
        for x, y in itertools.combinations(remaining, 2):
            require(
                is_negative_triangle(vertex, x, y) == (potential[x] ^ potential[y]),
                f"b=1 cut equation failed n={n},mask={signing}",
            )
        part_size = sum(potential.values())
        complement_size = n - 1 - part_size
        require(part_size >= 2 and complement_size >= 2, f"b=1 singleton side failed n={n},mask={signing}")
        for length, count in negative_counts.items():
            expected = part_size * complement_size * factorial(n - 3) // factorial(n - length)
            require(count == expected, f"b=1 cycle formula failed n={n},k={length},mask={signing}")
        return "b1"

    return "b0"


def mixed_four_set_count(
    n: int,
    negative_triangles: set[tuple[int, int, int]],
) -> int:
    mixed = 0
    for vertices in itertools.combinations(range(n), 4):
        local_count = sum(
            tuple(triangle) in negative_triangles
            for triangle in itertools.combinations(vertices, 3)
        )
        require(local_count in (0, 2, 4), f"two-graph parity failed n={n},vertices={vertices}")
        mixed += local_count == 2
    return mixed


def d3_target(n: int) -> int:
    return (n - 2) ** 2


def d4_target(n: int) -> int:
    return (n - 2) * (n * n - 6 * n + 10)


def audit_induction_algebra() -> None:
    for n in range(5, 65):
        d3_left = n * d3_target(n - 1) - d3_target(n)
        d3_right = (n - 4) * d3_target(n) + (n - 3) * (n - 4)
        require(d3_left == d3_right, f"D3 induction identity failed n={n}")
        require(
            factorial(n) // (8 * factorial(n - 4)) == 3 * comb(n, 4),
            f"four-cycle degree identity failed n={n}",
        )
        if n >= 5:
            require(
                factorial(n) // (10 * factorial(n - 5)) == 12 * comb(n, 5),
                f"five-cycle degree identity failed n={n}",
            )
        if n >= 9:
            margin = n * d4_target(n - 1) - ((n - 4) * d4_target(n) + comb(n, 3))
            expected_margin = 5 * (n - 8) * (n - 4) * (n - 3) // 6
            require(margin == expected_margin and margin > 0, f"D4 induction margin failed n={n}")


def direct_row(n: int) -> tuple[dict[str, object], int]:
    gauge_edges, edge_index = gauge_system(n)
    dimension = len(gauge_edges)
    character_count = 1 << dimension
    cycles_by_length = {
        length: simple_cycles(n, length, edge_index)
        for length in range(3, n + 1)
    }
    triangles = triangle_data(n, edge_index)
    expected_single_edges = single_edge_masks(n, gauge_edges, edge_index)
    b_counts = {"b0": 0, "b1": 0, "b_ge_2": 0}
    b_ge_2_masks: set[int] = set()
    d3_values: dict[int, int] = {}
    d4_values: dict[int, int] = {}

    for signing in range(1, character_count):
        negative_counts = direct_negative_counts(signing, cycles_by_length)
        negative_triangles = negative_triangle_support(signing, triangles)
        deletions = balanced_deletions(negative_triangles)
        category = audit_balanced_deletion_classification(
            n,
            signing,
            negative_triangles,
            deletions,
            negative_counts,
            triangles,
        )
        b_counts[category] += 1
        if category == "b_ge_2":
            b_ge_2_masks.add(signing)

        require(
            negative_counts[4] == 2 * mixed_four_set_count(n, negative_triangles),
            f"four-set rigidity identity failed n={n},mask={signing}",
        )
        d3_values[signing] = negative_counts[3] + negative_counts[4]
        if n >= 5:
            d4_values[signing] = d3_values[signing] + negative_counts[5]

    require(b_counts == EXPECTED_B_COUNTS[n], f"balanced-deletion census changed n={n}: {b_counts}")
    require(b_ge_2_masks == expected_single_edges, f"b>=2 equality masks changed n={n}")

    d3_minimum = min(d3_values.values())
    d3_minimizers = {mask for mask, value in d3_values.items() if value == d3_minimum}
    expected_d3_masks = set(expected_single_edges)
    if n == 4:
        expected_d3_masks.add(character_count - 1)
    require((d3_minimum, len(d3_minimizers)) == EXPECTED_D3[n], f"D3 base changed n={n}")
    require(d3_minimizers == expected_d3_masks, f"D3 equality masks changed n={n}")

    d4_summary = None
    if n >= 5:
        d4_minimum = min(d4_values.values())
        d4_minimizers = {mask for mask, value in d4_values.items() if value == d4_minimum}
        require((d4_minimum, len(d4_minimizers)) == EXPECTED_D4[n], f"D4 base changed n={n}")
        require(d4_minimizers == expected_single_edges, f"D4 equality masks changed n={n}")
        d4_summary = {"minimum": d4_minimum, "minimizers": len(d4_minimizers)}

    parity_gates = (character_count - 1) * sum(len(cycles) for cycles in cycles_by_length.values())
    row = {
        "n": n,
        "dimension": dimension,
        "characters": character_count,
        "cycles_by_length": {str(length): len(cycles) for length, cycles in cycles_by_length.items()},
        "parity_gates": parity_gates,
        "balanced_deletion_counts": b_counts,
        "d3": {"minimum": d3_minimum, "minimizers": len(d3_minimizers)},
        "d4": d4_summary,
    }
    return row, parity_gates


def n8_bridge_row() -> dict[str, object]:
    n = 8
    d3_lower_from_deletions = n * d3_target(n - 1)
    d3_boundary = (n - 4) * d3_target(n) + d3_target(n)
    require(d3_lower_from_deletions > d3_boundary, "D3 n=8 bridge lost strictness")

    d4_lower_from_deletions = n * d4_target(n - 1)
    d4_upper_at_target = (n - 4) * d4_target(n) + comb(n, 3)
    require(d4_lower_from_deletions == d4_upper_at_target == 680, "D4 n=8 equality bridge changed")
    negative_five_cycles_antibalanced = factorial(n) // (10 * factorial(n - 5))
    require(negative_five_cycles_antibalanced == 672, "antibalanced n=8 five-cycle count changed")

    return {
        "n": n,
        "certification": "deletion-average bridge from exact n=7 row",
        "d3": {"minimum": d3_target(n), "minimizers": comb(n, 2)},
        "d4": {"minimum": d4_target(n), "minimizers": comb(n, 2)},
        "d4_equality_lower": d4_lower_from_deletions,
        "d4_equality_upper": d4_upper_at_target,
        "antibalanced_c5": negative_five_cycles_antibalanced,
    }


def main() -> None:
    audit_induction_algebra()
    rows: list[dict[str, object]] = []
    total_parity_gates = 0
    for n in range(4, N_DIRECT_MAX + 1):
        row, parity_gates = direct_row(n)
        rows.append(row)
        total_parity_gates += parity_gates
    rows.append(n8_bridge_row())

    require(total_parity_gates == 38606835, "direct parity-gate aggregate changed")
    payload = {"rows": rows, "total_parity_gates": total_parity_gates}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(digest == EXPECTED_SEMANTIC_SHA256, f"semantic digest changed: {digest}")
    print("THM-4083 proof-local direct classification audit")
    print("direct_n_range=4..7 certified_base_range=4..8")
    print(f"total_direct_characters={sum(1 << comb(n - 1, 2) for n in range(4, 8))} direct_cycle_parity_gates={total_parity_gates}")
    print(f"rows={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: balanced deletions, direct cycle counts, four-set rigidity, D3/D4 bases, and equality masks")


if __name__ == "__main__":
    main()
