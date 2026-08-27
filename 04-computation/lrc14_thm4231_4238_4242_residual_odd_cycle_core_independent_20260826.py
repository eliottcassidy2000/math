#!/usr/bin/env python3
"""Standard-library referee for the current LRC(14) residual graph core."""

from __future__ import annotations

import contextlib
import hashlib
import io
import runpy
import sys
from collections import defaultdict


SOURCE = "04-computation/lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py"


def load_edges() -> list[tuple[int, int]]:
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        data = runpy.run_path(SOURCE)
    edges = data["aggregate_residual"]
    assert len(edges) == 181_126
    return edges


def articulation_points(adjacency: dict[int, set[int]]) -> list[int]:
    sys.setrecursionlimit(10_000)
    discovery: dict[int, int] = {}
    low: dict[int, int] = {}
    parent: dict[int, int] = {}
    cut: set[int] = set()
    clock = 0

    def visit(vertex: int) -> None:
        nonlocal clock
        clock += 1
        discovery[vertex] = low[vertex] = clock
        children = 0
        for neighbor in sorted(adjacency[vertex]):
            if neighbor not in discovery:
                parent[neighbor] = vertex
                children += 1
                visit(neighbor)
                low[vertex] = min(low[vertex], low[neighbor])
                if vertex not in parent and children > 1:
                    cut.add(vertex)
                if vertex in parent and low[neighbor] >= discovery[vertex]:
                    cut.add(vertex)
            elif parent.get(vertex) != neighbor:
                low[vertex] = min(low[vertex], discovery[neighbor])

    for vertex in sorted(adjacency):
        if vertex not in discovery:
            visit(vertex)
    return sorted(cut)


def connected(adjacency: dict[int, set[int]]) -> bool:
    start = next(iter(adjacency))
    seen = {start}
    stack = [start]
    while stack:
        vertex = stack.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == len(adjacency)


def main() -> None:
    edges = load_edges()
    adjacency: dict[int, set[int]] = defaultdict(set)
    for left, right in edges:
        assert left < right
        adjacency[left].add(right)
        adjacency[right].add(left)

    assert len(adjacency) == 739
    assert connected(adjacency)
    leaves = sorted(vertex for vertex in adjacency if len(adjacency[vertex]) == 1)
    assert leaves == [756, 760]

    triangle_by_edge: dict[tuple[int, int], int] = {}
    triangle_sum = 0
    triangle_vertices: set[int] = set()
    for edge in edges:
        left, right = edge
        common = adjacency[left] & adjacency[right]
        triangle_by_edge[edge] = len(common)
        triangle_sum += len(common)
        if common:
            triangle_vertices.add(left)
            triangle_vertices.add(right)
            triangle_vertices.update(common)
    assert triangle_sum % 3 == 0
    triangles = triangle_sum // 3
    assert triangles == 30_912_074
    assert len(triangle_vertices) == 737

    triangle_free = [edge for edge in edges if triangle_by_edge[edge] == 0]
    assert triangle_free == [(616, 756), (616, 760)]

    # A triangle-covered edge is never a bridge. The two remaining edges end
    # in the only leaves, so the bridge set is exactly triangle_free.
    bridges = [(616, 756), (616, 760)]
    assert triangle_free == bridges
    cuts = articulation_points(adjacency)
    assert cuts == [616]

    core_vertices = set(adjacency) - set(leaves)
    core = {
        vertex: adjacency[vertex] & core_vertices
        for vertex in core_vertices
    }
    assert len(core) == 737
    assert connected(core)
    assert articulation_points(core) == []

    fixed = 50
    neighbors = adjacency[fixed]
    fixed_edges = sum(
        1 for left, right in edges if left in neighbors and right in neighbors
    )
    fixed_triangle_sum = sum(
        len(adjacency[left] & adjacency[right] & neighbors)
        for left, right in edges
        if left in neighbors and right in neighbors
    )
    assert len(neighbors) == 556
    assert fixed_edges == 150_861
    assert fixed_triangle_sum % 3 == 0
    assert fixed_triangle_sum // 3 == 26_807_003

    digest = hashlib.sha256(
        "".join(
            f"{left},{right},{triangle_by_edge[(left, right)]}\n"
            for left, right in sorted(edges)
        ).encode("ascii")
    ).hexdigest()
    assert digest == "d682ec7b70b7eecc641559e26d0afa1e23b0ede3ff90c025a637e6c02e208623"

    print("LRC14_RESIDUAL_ODD_CYCLE_CORE_INDEPENDENT")
    print(f"vertices={len(adjacency)} edges={len(edges)} connected=PASS")
    print(f"triangles={triangles} triangle_vertices={len(triangle_vertices)}")
    print(f"triangle_free_edges={triangle_free}")
    print(f"bridges={bridges} articulation_points={cuts}")
    print("biconnected_core_vertices=737 core_articulation_points=0")
    print(
        f"fixed50_degree={len(neighbors)} induced_edges={fixed_edges} "
        f"induced_triangles={fixed_triangle_sum // 3}"
    )
    print(f"edge_triangle_sha256={digest}")
    print("verdict=odd_cycles_schedule_obligations_not_safety")


if __name__ == "__main__":
    main()
