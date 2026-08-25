#!/usr/bin/env python3
"""Independent exhaustive character/FWHT audit for THM-4083.

Unlike the proof-local audit, this path does not inspect balanced deletions or
use the induction.  It independently generates cycle characters from ordered
vertex tuples and exhausts every root-gauged signing through n=8 by exact
integer Walsh transforms.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import comb, factorial

import numpy


N_MAX = 8
EXPECTED_SEMANTIC_SHA256 = "6f814b904df999b9cf432a981570ddd3a5aa30a47a38a0a64a161521c5eb64ec"
EXPECTED_D3 = {4: (4, 7), 5: (9, 10), 6: (16, 15), 7: (25, 21), 8: (36, 28)}
EXPECTED_D4 = {5: (15, 10), 6: (40, 15), 7: (85, 21), 8: (156, 28)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gauge_index(n: int) -> dict[tuple[int, int], int]:
    return {
        edge: index
        for index, edge in enumerate(itertools.combinations(range(1, n), 2))
    }


def ordered_cycle_masks(n: int, length: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    masks: list[int] = []
    for order in itertools.permutations(range(n), length):
        if order[0] != min(order) or order[1] > order[-1]:
            continue
        mask = 0
        for index, a in enumerate(order):
            b = order[(index + 1) % length]
            if a == 0 or b == 0:
                continue
            edge = (a, b) if a < b else (b, a)
            mask ^= 1 << edge_index[edge]
        masks.append(mask)
    expected = factorial(n) // (2 * length * factorial(n - length))
    require(len(masks) == expected, f"ordered cycle count failed n={n},k={length}")
    require(len(set(masks)) == expected, f"ordered cycle masks collided n={n},k={length}")
    return masks


def integer_fwht(values: numpy.ndarray) -> tuple[numpy.ndarray, int]:
    transformed = values.copy()
    width = 1
    butterflies = 0
    while width < len(transformed):
        blocks = transformed.reshape(-1, 2, width)
        left = blocks[:, 0, :].copy()
        right = blocks[:, 1, :].copy()
        blocks[:, 0, :] = left + right
        blocks[:, 1, :] = left - right
        butterflies += len(transformed) // 2
        width *= 2
    return transformed, butterflies


def single_edge_masks(n: int, edge_index: dict[tuple[int, int], int]) -> set[int]:
    masks = {1 << index for index in edge_index.values()}
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


def minimum_and_masks(values: numpy.ndarray) -> tuple[int, set[int]]:
    minimum = int(values[1:].min())
    masks = {
        int(index)
        for index in numpy.flatnonzero(values == minimum)
        if index
    }
    return minimum, masks


def row_for_n(n: int) -> tuple[dict[str, object], int, int]:
    edge_index = gauge_index(n)
    dimension = len(edge_index)
    character_count = 1 << dimension
    expected_single_edges = single_edge_masks(n, edge_index)
    cumulative = numpy.zeros(character_count, dtype=numpy.int64)
    cycles_by_length: dict[str, int] = {}
    d3_summary = None
    d4_summary = None
    total_cycles = 0
    total_butterflies = 0

    for length in range(3, min(5, n) + 1):
        masks = ordered_cycle_masks(n, length, edge_index)
        total_cycles += len(masks)
        cycles_by_length[str(length)] = len(masks)
        frequency = numpy.zeros(character_count, dtype=numpy.int64)
        frequency[numpy.array(masks, dtype=numpy.int64)] = 1
        signed_sums, butterflies = integer_fwht(frequency)
        total_butterflies += butterflies
        require(bool(numpy.all((len(masks) - signed_sums) % 2 == 0)), f"parity divisibility failed n={n},k={length}")
        negative_counts = (len(masks) - signed_sums) // 2
        require(int(negative_counts[0]) == 0, f"balanced character failed n={n},k={length}")
        cumulative += negative_counts

        if length == 4:
            minimum, minimizers = minimum_and_masks(cumulative)
            expected_masks = set(expected_single_edges)
            if n == 4:
                expected_masks.add(character_count - 1)
            require((minimum, len(minimizers)) == EXPECTED_D3[n], f"D3 exhaustive minimum failed n={n}")
            require(minimizers == expected_masks, f"D3 exhaustive equality masks failed n={n}")
            d3_summary = {"minimum": minimum, "minimizers": len(minimizers)}

        if length == 5:
            minimum, minimizers = minimum_and_masks(cumulative)
            require((minimum, len(minimizers)) == EXPECTED_D4[n], f"D4 exhaustive minimum failed n={n}")
            require(minimizers == expected_single_edges, f"D4 exhaustive equality masks failed n={n}")
            d4_summary = {"minimum": minimum, "minimizers": len(minimizers)}

    require(d3_summary is not None, f"D3 row missing n={n}")
    if n >= 5:
        require(d4_summary is not None, f"D4 row missing n={n}")
    row = {
        "n": n,
        "dimension": dimension,
        "characters": character_count,
        "cycles_by_length": cycles_by_length,
        "d3": d3_summary,
        "d4": d4_summary,
    }
    return row, total_cycles, total_butterflies


def main() -> None:
    rows: list[dict[str, object]] = []
    total_characters = 0
    total_cycles = 0
    total_butterflies = 0
    for n in range(4, N_MAX + 1):
        row, cycle_count, butterflies = row_for_n(n)
        rows.append(row)
        total_characters += int(row["characters"])
        total_cycles += cycle_count
        total_butterflies += butterflies

    require(total_characters == 2131016, "character aggregate changed")
    require(total_cycles == 1511, "cycle-vector aggregate changed")
    require(total_butterflies == 66813528, "FWHT butterfly aggregate changed")
    payload = {
        "rows": rows,
        "total_characters": total_characters,
        "total_cycles": total_cycles,
        "total_butterflies": total_butterflies,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(digest == EXPECTED_SEMANTIC_SHA256, f"semantic digest changed: {digest}")
    print("THM-4083 independent exhaustive character/FWHT audit")
    print("n_range=4..8 terminal_lengths=4,5")
    print(f"total_characters={total_characters} cycle_vectors={total_cycles} integer_butterflies={total_butterflies}")
    print(f"rows={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: exhaustive D3/D4 minima and labelled equality multiplicities through n=8")


if __name__ == "__main__":
    main()
