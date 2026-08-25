#!/usr/bin/env python3
"""Independent dual/Fourier audit for THM-4078.

Characters are gauged signed complete graphs. Exact Walsh transforms recover
all cycle-convolution eigenvalues and test the proved triangle extremum plus
the explicitly conjectural cumulative diameter-layer pattern through n=8.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import comb, factorial

import numpy


N_MAX = 8
ORBIT_N_MAX = 7
EXPECTED_ORBITS = {3: 2, 4: 3, 5: 7, 6: 16, 7: 54}
EXPECTED_PREFIX_MINIMA = {
    3: [1],
    4: [2, 4],
    5: [3, 9, 15],
    6: [4, 16, 40, 64],
    7: [5, 25, 85, 205, 325],
    8: [6, 36, 156, 516, 1236, 1956],
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gauge_edges(n: int) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    edges = list(itertools.combinations(range(1, n), 2))
    return edges, {edge: index for index, edge in enumerate(edges)}


def simple_cycle_vectors(n: int, length: int, gauge_index: dict[tuple[int, int], int]) -> list[int]:
    vectors: list[int] = []
    for vertices in itertools.combinations(range(n), length):
        first = vertices[0]
        for tail in itertools.permutations(vertices[1:]):
            if tail[0] > tail[-1]:
                continue
            order = (first,) + tail
            mask = 0
            for index in range(length):
                a, b = order[index], order[(index + 1) % length]
                if a == 0 or b == 0:
                    continue
                edge = (a, b) if a < b else (b, a)
                mask |= 1 << gauge_index[edge]
            vectors.append(mask)
    expected = factorial(n) // (2 * length * factorial(n - length))
    require(len(vectors) == expected and len(set(vectors)) == expected, f"dual cycle generation failed n={n},k={length}")
    return vectors


def walsh_transform(values: numpy.ndarray) -> numpy.ndarray:
    transformed = values.copy()
    width = 1
    size = len(transformed)
    while width < size:
        blocks = transformed.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return transformed


def single_edge_character_masks(n: int, gauge_index: dict[tuple[int, int], int]) -> set[int]:
    masks = {1 << index for index in range(len(gauge_index))}
    for vertex in range(1, n):
        mask = 0
        for other in range(1, n):
            if other == vertex:
                continue
            edge = (vertex, other) if vertex < other else (other, vertex)
            mask |= 1 << gauge_index[edge]
        masks.add(mask)
    return masks


def regauge_after_permutation(mask: int, n: int, permutation: tuple[int, ...], gauge_edges_list: list[tuple[int, int]], gauge_index: dict[tuple[int, int], int]) -> int:
    signs = [[0] * n for _ in range(n)]
    for index, (a, b) in enumerate(gauge_edges_list):
        if (mask >> index) & 1:
            x, y = permutation[a], permutation[b]
            signs[x][y] = 1
            signs[y][x] = 1
    switch = [signs[0][vertex] for vertex in range(n)]
    output = 0
    for edge, index in gauge_index.items():
        a, b = edge
        value = signs[a][b] ^ switch[a] ^ switch[b]
        if value:
            output |= 1 << index
    return output


def dual_orbit_count(n: int, gauge_edges_list: list[tuple[int, int]], gauge_index: dict[tuple[int, int], int]) -> int:
    dimension = len(gauge_edges_list)
    size = 1 << dimension
    transforms: list[list[int]] = []
    for swap_at in range(n - 1):
        permutation = list(range(n))
        permutation[swap_at], permutation[swap_at + 1] = permutation[swap_at + 1], permutation[swap_at]
        transforms.append([
            regauge_after_permutation(mask, n, tuple(permutation), gauge_edges_list, gauge_index)
            for mask in range(size)
        ])
    seen = bytearray(size)
    count = 0
    for seed in range(size):
        if seen[seed]:
            continue
        count += 1
        seen[seed] = 1
        stack = [seed]
        while stack:
            current = stack.pop()
            for transform in transforms:
                neighbor = transform[current]
                if not seen[neighbor]:
                    seen[neighbor] = 1
                    stack.append(neighbor)
    return count


def main() -> None:
    rows: list[dict[str, object]] = []
    total_characters = 0
    total_cycle_vectors = 0
    orbit_counts: dict[int, int] = {}

    for n in range(3, N_MAX + 1):
        gauge_edges_list, gauge_index = gauge_edges(n)
        dimension = len(gauge_edges_list)
        size = 1 << dimension
        total_characters += size
        single_edge_masks = single_edge_character_masks(n, gauge_index)
        cumulative = numpy.zeros(size, dtype=numpy.int64)
        prefix_minima: list[int] = []
        minimizer_counts: list[int] = []

        for length in range(3, n + 1):
            vectors = simple_cycle_vectors(n, length, gauge_index)
            total_cycle_vectors += len(vectors)
            indicator = numpy.zeros(size, dtype=numpy.int64)
            indicator[numpy.array(vectors, dtype=numpy.int64)] = 1
            eigenvalues = walsh_transform(indicator)
            negative_counts = (len(vectors) - eigenvalues) // 2
            if length == 3:
                minimum = int(negative_counts[1:].min())
                minimizers = set(int(index) for index in numpy.flatnonzero(negative_counts == minimum) if index)
                require(minimum == n - 2, f"triangle minimum failed n={n}")
                require(minimizers == single_edge_masks, f"triangle equality orbit failed n={n}")
                all_negative = size - 1
                require(int(eigenvalues[all_negative]) == -len(vectors), f"bipartite character failed n={n}")
            cumulative += negative_counts
            prefix_minimum = int(cumulative[1:].min())
            prefix_minima.append(prefix_minimum)
            minimizer_counts.append(int(numpy.count_nonzero(cumulative[1:] == prefix_minimum)))

        require(prefix_minima == EXPECTED_PREFIX_MINIMA[n], f"prefix minima changed n={n}")
        expected_single = [
            sum(factorial(n - 2) // factorial(n - length) for length in range(3, terminal + 1))
            for terminal in range(3, n + 1)
        ]
        require(prefix_minima == expected_single, f"single-edge prefix formula failed n={n}")
        if n != 4:
            require(all(count == len(single_edge_masks) for count in minimizer_counts), f"prefix minimizer count failed n={n}")
        else:
            require(minimizer_counts == [6, 7], "n=4 antibalanced tie changed")

        if n <= ORBIT_N_MAX:
            orbit_count = dual_orbit_count(n, gauge_edges_list, gauge_index)
            require(orbit_count == EXPECTED_ORBITS[n], f"dual orbit count changed n={n}")
            orbit_counts[n] = orbit_count

        rows.append(
            {
                "n": n,
                "dimension": dimension,
                "characters": size,
                "prefix_minima": prefix_minima,
                "minimizer_counts": minimizer_counts,
            }
        )

    require(total_characters == 2131018, "character aggregate changed")
    require(total_cycle_vectors == 9432, "cycle-vector aggregate changed")
    digest = hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4078 independent signed-dual/Walsh audit")
    print(f"n_range=3..{N_MAX} total_characters={total_characters} cycle_vectors={total_cycle_vectors}")
    print(f"dual_orbit_counts={json.dumps(orbit_counts, sort_keys=True, separators=(',', ':'))}")
    print(f"prefix_rows={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: dual spectrum and FINITE-EXACT prefix conjecture through n=8")


if __name__ == "__main__":
    main()
