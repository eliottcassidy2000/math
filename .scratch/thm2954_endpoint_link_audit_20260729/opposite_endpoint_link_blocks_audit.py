#!/usr/bin/env python3
"""Finite hostile/control audit for the proposed THM-2954 support theorem.

This checks only the abstract arithmetic coincidence graph

    u ~ v  iff  14*gcd(u,v) divides u+v.

It does not infer a strict-open cover transition from an endpoint
coincidence.  MISTAKE-146 shows that such a transition would need an
additional protector/overlap sidecar.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from collections import defaultdict, deque
from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def block_address(value: int) -> tuple[int, int, int]:
    a = valuation(value, 2)
    b = valuation(value, 7)
    core = value // (2**a * 7**b)
    require(core % 2 and core % 7, "core is not a 14-unit")
    return a, b, core % 7


def arithmetic_link(first: int, second: int) -> bool:
    return (
        first != second
        and (first + second) % (14 * math.gcd(first, second)) == 0
    )


def block_link(first: int, second: int) -> bool:
    a, b, residue = block_address(first)
    c, d, other = block_address(second)
    return first != second and a == c and b == d and (residue + other) % 7 == 0


def is_bipartite(vertices: tuple[int, ...], edges: set[tuple[int, int]]) -> bool:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    color: dict[int, int] = {}
    for root in vertices:
        if root in color:
            continue
        color[root] = 0
        queue = deque([root])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if neighbor not in color:
                    color[neighbor] = 1 - color[vertex]
                    queue.append(neighbor)
                elif color[neighbor] == color[vertex]:
                    return False
    return True


def cyclic_orbit(difference: int) -> set[tuple[int, int]]:
    return {
        tuple(sorted((index, (index + difference) % 7)))
        for index in range(7)
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=4_096)
    args = parser.parse_args()
    require(args.bound >= 14, "bound is too small for the positive control")

    pair_count = 0
    edge_count = 0
    block_edges: dict[tuple[int, int, int], int] = defaultdict(int)
    digest = hashlib.sha256(b"THM2954/opposite-endpoint-link-blocks/v1\n")
    for first in range(1, args.bound + 1):
        first_address = block_address(first)
        for second in range(first + 1, args.bound + 1):
            pair_count += 1
            arithmetic = arithmetic_link(first, second)
            block = block_link(first, second)
            require(
                arithmetic == block,
                f"block equivalence failed at ({first},{second})",
            )
            if arithmetic:
                edge_count += 1
                a, b, residue = first_address
                pair_residue = min(residue, (-residue) % 7)
                block_edges[(a, b, pair_residue)] += 1
                digest.update(f"{first},{second}\n".encode())

    # A physical non-equivariant support can be nonempty and bipartite.
    physical_vertices = (1, 13, 2, 4, 8, 16, 32)
    physical_edges = {
        (first, second)
        for index, first in enumerate(physical_vertices)
        for second in physical_vertices[index + 1 :]
        if arithmetic_link(first, second)
    }
    require(physical_edges == {(1, 13)}, "positive physical control changed")
    require(
        is_bipartite(physical_vertices, physical_edges),
        "arithmetic support lost bipartiteness",
    )

    # Every C7-invariant undirected graph is a union of the three difference
    # orbits.  Every nonempty one contains a seven-cycle, hence is nonbipartite.
    invariant_controls: list[tuple[int, int, bool]] = []
    for mask in range(8):
        edges: set[tuple[int, int]] = set()
        for difference in (1, 2, 3):
            if mask & (1 << (difference - 1)):
                edges.update(cyclic_orbit(difference))
        bipartite = is_bipartite(tuple(range(7)), edges)
        require(bipartite == (mask == 0), "C7 invariant control changed")
        invariant_controls.append((mask, len(edges), bipartite))

    # Qualitative sharpness.  Two edges of the transition triangle are exact
    # zero links, while the third is a positive opposite-endpoint overlap that
    # tends to zero.  Hence triangle-freeness supplies no scale-free width.
    near_link_rows: list[tuple[int, int, int, F]] = []
    previous_width: F | None = None
    for multiplier in (1, 2, 5, 20, 100, 500):
        first = 98 * multiplier + 1
        second = 98 * multiplier + 15
        third = 13
        exit_index = 6 * multiplier
        entry_index = 6 * multiplier + 1
        numerator = (
            (14 * exit_index + 1) * second
            - (14 * entry_index - 1) * first
        )
        width = F(numerator, 14 * first * second)
        require(
            arithmetic_link(first, third)
            and arithmetic_link(second, third)
            and not arithmetic_link(first, second),
            "near-link triangle type changed",
        )
        require(
            numerator == 2 and width == F(1, 7 * first * second),
            "near-link width identity changed",
        )
        require(
            previous_width is None or width < previous_width,
            "near-link widths no longer decrease",
        )
        previous_width = width
        near_link_rows.append((first, second, third, width))

    # A physical four-aligned/three-drift packet has a nontrivial type
    # partition.  No transitive C7 shift preserves its four-element side.
    aligned_types = {0, 1, 2, 3}
    require(
        all(
            {(vertex + shift) % 7 for vertex in aligned_types} != aligned_types
            for shift in range(1, 7)
        ),
        "a nontrivial C7 shift preserved the 4/3 type partition",
    )

    print("THM-2954 opposite-endpoint abstract support audit")
    print(
        f"pair_universe_1_to_{args.bound}={pair_count};"
        f"arithmetic_edges={edge_count};"
        f"nonempty_blocks={len(block_edges)}"
    )
    print(
        "block_law=(nu2_equal,nu7_equal,unit_residues_opposite_mod7);"
        "residue_pairs=((1,6),(2,5),(3,4))"
    )
    print(f"edge_digest={digest.hexdigest()}")
    print(
        f"non_equivariant_physical_control={physical_vertices};"
        f"edges={tuple(sorted(physical_edges))};bipartite=TRUE"
    )
    print(f"C7_invariant_controls={tuple(invariant_controls)}")
    print("every_nonempty_C7_invariant_support_contains_odd_cycle=PASS")
    print(f"qualitative_sharpness_family={tuple(near_link_rows)}")
    print("uniform_positive_overlap_lower_bound=REFUTED")
    print("four_aligned_three_drift_C7_type_equivariance=IMPOSSIBLE")
    print("pointwise_cover_consequence=NONE_WITHOUT_PROTECTOR_SIDECAR")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
