#!/usr/bin/env python3
"""Primary exact matching-character audit for THM-4084.

This path enumerates labelled matchings and unoriented simple cycles directly.
It checks the all-layer inclusion--exclusion formula, the D5 positivity
firewall, switching minimality, the two frustration-index-two shapes, and the
sharp n=4 D3 boundary.  All arithmetic is integral.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import comb, factorial


EXPECTED_SEMANTIC_SHA256 = "9c9b917b84055a6411abae78026ed1915d6f7a6358a399078368422f2611de9f"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def edge_system(n: int) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    edges = list(itertools.combinations(range(n), 2))
    return edges, {item: index for index, item in enumerate(edges)}


def cycle_masks(n: int, length: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    masks: list[int] = []
    for vertices in itertools.combinations(range(n), length):
        first = vertices[0]
        for tail in itertools.permutations(vertices[1:]):
            if tail[0] > tail[-1]:
                continue
            order = (first,) + tail
            mask = 0
            for index, a in enumerate(order):
                b = order[(index + 1) % length]
                mask |= 1 << edge_index[edge(a, b)]
            masks.append(mask)
    expected = factorial(n) // (2 * length * factorial(n - length))
    require(len(masks) == expected, f"cycle count changed n={n},k={length}")
    require(len(set(masks)) == expected, f"cycle collision n={n},k={length}")
    return masks


def matching_masks(n: int, size: int, edges: list[tuple[int, int]], edge_index: dict[tuple[int, int], int]) -> list[int]:
    masks: list[int] = []
    for chosen in itertools.combinations(edges, size):
        vertices = {vertex for item in chosen for vertex in item}
        if len(vertices) != 2 * size:
            continue
        mask = sum(1 << edge_index[item] for item in chosen)
        masks.append(mask)
    expected = factorial(n) // (2**size * factorial(size) * factorial(n - 2 * size))
    require(len(masks) == expected, f"matching count changed n={n},r={size}")
    return masks


def negative_count(signing: int, cycles: list[int]) -> int:
    return sum((signing & cycle).bit_count() & 1 for cycle in cycles)


def profile_formula(n: int, length: int, size: int) -> int:
    total = 0
    for chosen_edges in range(1, min(size, length // 2) + 1):
        total += (
            (-4) ** (chosen_edges - 1)
            * comb(size, chosen_edges)
            * factorial(length - chosen_edges - 1)
            * comb(n - 2 * chosen_edges, length - 2 * chosen_edges)
        )
    return total


def edge_profile(n: int, length: int) -> int:
    return factorial(n - 2) // factorial(n - length)


def d5_edge_value(n: int) -> int:
    return sum(edge_profile(n, length) for length in range(3, 7))


def d5_overlap_value(n: int) -> int:
    return 3 * n * n - 25 * n + 53


def d5_matching_value(n: int, size: int) -> int:
    return (
        size * d5_edge_value(n)
        - 4 * comb(size, 2) * d5_overlap_value(n)
        + 32 * comb(size, 3)
    )


def cut_mask(n: int, subset: tuple[int, ...], edge_index: dict[tuple[int, int], int]) -> int:
    chosen = set(subset)
    mask = 0
    for a, b in edge_index:
        if (a in chosen) != (b in chosen):
            mask |= 1 << edge_index[(a, b)]
    return mask


def canonical_matching_mask(size: int, edge_index: dict[tuple[int, int], int]) -> int:
    return sum(1 << edge_index[(2 * index, 2 * index + 1)] for index in range(size))


def audit_symbolic_range() -> tuple[int, int]:
    positivity_gates = 0
    switching_gates = 0
    for n in range(6, 129):
        edge_value = d5_edge_value(n)
        overlap = d5_overlap_value(n)
        require(
            edge_value == (n - 2) * (n**3 - 11 * n * n + 41 * n - 50),
            f"edge polynomial changed n={n}",
        )
        require(overlap == 1 + 2 * (n - 4) + 3 * (n - 4) * (n - 5), f"overlap polynomial changed n={n}")
        for size in range(1, n // 2 + 1):
            value = d5_matching_value(n, size)
            formula_sum = sum(profile_formula(n, length, size) for length in range(3, 7))
            require(value == formula_sum, f"D5 profile sum changed n={n},r={size}")
            gap = value - edge_value
            defect_numerator = 3 * edge_value - 6 * size * overlap + 16 * size * (size - 2)
            require(3 * gap == (size - 1) * defect_numerator, f"gap factorization changed n={n},r={size}")
            require(gap == 0 if size == 1 else gap > 0, f"D5 firewall failed n={n},r={size}")
            positivity_gates += 4

    for n in range(6, 13):
        edges, edge_index = edge_system(n)
        for size in range(1, n // 2 + 1):
            matching = canonical_matching_mask(size, edge_index)
            for subset_size in range(1, n):
                for subset in itertools.combinations(range(1, n), subset_size):
                    switched_weight = (matching ^ cut_mask(n, subset, edge_index)).bit_count()
                    require(switched_weight > size, f"matching was not switching-minimal n={n},r={size},S={subset}")
                    switching_gates += 1
        require(len(edges) == comb(n, 2), f"edge universe changed n={n}")
    return positivity_gates, switching_gates


def audit_direct_labelled() -> tuple[int, int, list[dict[str, int]]]:
    matching_gates = 0
    parity_gates = 0
    rows: list[dict[str, int]] = []
    for n in range(6, 10):
        edges, edge_index = edge_system(n)
        cycles = {length: cycle_masks(n, length, edge_index) for length in range(3, 7)}
        row_matchings = 0
        for size in range(1, n // 2 + 1):
            for signing in matching_masks(n, size, edges, edge_index):
                for length in range(3, 7):
                    actual = negative_count(signing, cycles[length])
                    expected = profile_formula(n, length, size)
                    require(actual == expected, f"direct profile failed n={n},r={size},k={length}")
                    matching_gates += 1
                    parity_gates += len(cycles[length])
                require(
                    sum(negative_count(signing, cycles[length]) for length in range(3, 7))
                    == d5_matching_value(n, size),
                    f"direct D5 sum failed n={n},r={size}",
                )
                matching_gates += 1
                row_matchings += 1

        adjacent = (1 << edge_index[(0, 1)]) | (1 << edge_index[(0, 2)])
        adjacent_sum = 0
        for length in range(3, 7):
            actual = negative_count(adjacent, cycles[length])
            expected = 2 * edge_profile(n, length) - 2 * factorial(n - 3) // factorial(n - length)
            require(actual == expected, f"adjacent-pair profile failed n={n},k={length}")
            adjacent_sum += actual
            matching_gates += 1
            parity_gates += len(cycles[length])
        require(adjacent_sum > d5_edge_value(n), f"adjacent frustration-two firewall failed n={n}")
        rows.append(
            {
                "n": n,
                "labelled_matchings": row_matchings,
                "edge_value": d5_edge_value(n),
                "matching_two_value": d5_matching_value(n, 2),
                "adjacent_two_value": adjacent_sum,
            }
        )
    return matching_gates, parity_gates, rows


def audit_boundary() -> tuple[int, int]:
    n = 4
    edges, edge_index = edge_system(n)
    cycles = {length: cycle_masks(n, length, edge_index) for length in (3, 4)}
    matching = canonical_matching_mask(2, edge_index)
    matching_value = sum(negative_count(matching, cycles[length]) for length in (3, 4))
    edge_value = sum(edge_profile(n, length) for length in (3, 4))
    require(matching_value == edge_value == 4, "n=4 D3 matching tie changed")
    require(len(edges) == 6, "K4 edge count changed")
    return matching_value, edge_value


def main() -> None:
    positivity_gates, switching_gates = audit_symbolic_range()
    matching_gates, parity_gates, rows = audit_direct_labelled()
    boundary_matching, boundary_edge = audit_boundary()
    payload = {
        "boundary": [boundary_matching, boundary_edge],
        "matching_gates": matching_gates,
        "parity_gates": parity_gates,
        "positivity_gates": positivity_gates,
        "rows": rows,
        "switching_gates": switching_gates,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    require(digest == EXPECTED_SEMANTIC_SHA256, f"semantic digest changed: {digest}")
    print("THM-4084 primary labelled-matching cycle audit")
    print(f"direct_n_range=6..9 symbolic_n_range=6..128")
    print(f"matching_profile_gates={matching_gates} direct_cycle_parity_gates={parity_gates}")
    print(f"positivity_factorization_gates={positivity_gates} switching_minimality_gates={switching_gates}")
    print(f"rows={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"n4_D3_boundary=matching:{boundary_matching},edge:{boundary_edge}")
    print(f"semantic_sha256={digest}")
    print("PASS: all-layer matching profile, D5 strict firewall, frustration-two shapes, and sharp D3 boundary")


if __name__ == "__main__":
    main()
