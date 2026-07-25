#!/usr/bin/env python3
"""Exact audit for THM-2281.

The theorem itself is algebraic.  This script independently checks the
direct inherited-color matching sum against the hafnian of the selected
colored pair kernel on a finite exhaustive bank, verifies the two hostile
two-vertex quotients, and checks the leaf-forcing entry probe.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import product


COLORS = ("r", "b")


@dataclass(frozen=True)
class Edge:
    u: int
    v: int
    color_u: str
    color_v: str
    weight: int

    def normalized(self) -> "Edge":
        if self.u < self.v:
            return self
        return Edge(self.v, self.u, self.color_v, self.color_u, self.weight)


def perfect_pairings(vertices: tuple[int, ...]):
    if not vertices:
        yield ()
        return
    u = vertices[0]
    for index in range(1, len(vertices)):
        v = vertices[index]
        rest = vertices[1:index] + vertices[index + 1 :]
        for tail in perfect_pairings(rest):
            yield ((u, v),) + tail


def colored_pair_kernel(
    vertices: tuple[int, ...], edges: tuple[Edge, ...]
) -> dict[tuple[int, str, int, str], int]:
    kernel = {
        (u, a, v, b): 0
        for u in vertices
        for v in vertices
        if u < v
        for a in COLORS
        for b in COLORS
    }
    for raw_edge in edges:
        edge = raw_edge.normalized()
        key = (edge.u, edge.color_u, edge.v, edge.color_v)
        kernel[key] = kernel.get(key, 0) + edge.weight
    return kernel


def hafnian_amplitude(
    vertices: tuple[int, ...],
    edges: tuple[Edge, ...],
    coloring: tuple[str, ...],
) -> int:
    kernel = colored_pair_kernel(vertices, edges)
    total = 0
    for pairing in perfect_pairings(vertices):
        term = 1
        for u, v in pairing:
            term *= kernel[(u, coloring[u], v, coloring[v])]
        total += term
    return total


def direct_amplitude(
    vertices: tuple[int, ...],
    edges: tuple[Edge, ...],
    coloring: tuple[str, ...],
) -> int:
    normalized = tuple(edge.normalized() for edge in edges)

    def recurse(unmatched: tuple[int, ...]) -> int:
        if not unmatched:
            return 1
        u = unmatched[0]
        remaining = set(unmatched[1:])
        total = 0
        for edge in normalized:
            if edge.u == u and edge.v in remaining:
                v = edge.v
                if (
                    edge.color_u == coloring[u]
                    and edge.color_v == coloring[v]
                ):
                    tail = tuple(x for x in unmatched[1:] if x != v)
                    total += edge.weight * recurse(tail)
            elif edge.v == u and edge.u in remaining:
                # This branch is retained for robustness if the vertex tuple
                # is not in increasing order.
                v = edge.u
                if (
                    edge.color_v == coloring[u]
                    and edge.color_u == coloring[v]
                ):
                    tail = tuple(x for x in unmatched[1:] if x != v)
                    total += edge.weight * recurse(tail)
        return total

    return recurse(vertices)


def response_vector(
    vertices: tuple[int, ...], edges: tuple[Edge, ...]
) -> tuple[int, ...]:
    return tuple(
        direct_amplitude(vertices, edges, coloring)
        for coloring in product(COLORS, repeat=len(vertices))
    )


def exhaustive_four_vertex_audit() -> tuple[int, str]:
    vertices = (0, 1, 2, 3)
    pairs = tuple((u, v) for u in vertices for v in vertices if u < v)
    local_options = (
        (),
        (("r", "r", 1),),
        (("b", "b", 1),),
        (("r", "b", -1),),
        (("b", "r", 2),),
        (("r", "r", 1), ("r", "r", -1)),
    )
    checked = 0
    digest = sha256()
    for option_word in product(range(len(local_options)), repeat=len(pairs)):
        edges_list: list[Edge] = []
        for (u, v), option_index in zip(pairs, option_word):
            for color_u, color_v, weight in local_options[option_index]:
                edges_list.append(Edge(u, v, color_u, color_v, weight))
        edges = tuple(edges_list)
        for coloring in product(COLORS, repeat=4):
            direct = direct_amplitude(vertices, edges, coloring)
            hafnian = hafnian_amplitude(vertices, edges, coloring)
            assert direct == hafnian
            digest.update(f"{option_word}:{coloring}:{direct}\n".encode())
            checked += 1
    return checked, digest.hexdigest()


def selector_erasure_collision():
    vertices = (0, 1)
    first = (
        Edge(0, 1, "r", "r", 1),
        Edge(0, 1, "b", "b", 1),
    )
    second = (
        Edge(0, 1, "r", "r", 1),
        Edge(0, 1, "r", "r", 1),
    )
    assert tuple(edge.weight for edge in first) == tuple(
        edge.weight for edge in second
    )
    assert response_vector(vertices, first) == (1, 0, 0, 1)
    assert response_vector(vertices, second) == (2, 0, 0, 0)
    return response_vector(vertices, first), response_vector(vertices, second)


def phase_erasure_collision():
    vertices = (0, 1)
    cancelling = (
        Edge(0, 1, "r", "r", 1),
        Edge(0, 1, "r", "r", -1),
    )
    reinforcing = (
        Edge(0, 1, "r", "r", 1),
        Edge(0, 1, "r", "r", 1),
    )
    assert tuple(
        (edge.u, edge.v, edge.color_u, edge.color_v)
        for edge in cancelling
    ) == tuple(
        (edge.u, edge.v, edge.color_u, edge.color_v)
        for edge in reinforcing
    )
    assert response_vector(vertices, cancelling) == (0, 0, 0, 0)
    assert response_vector(vertices, reinforcing) == (2, 0, 0, 0)
    return (
        response_vector(vertices, cancelling),
        response_vector(vertices, reinforcing),
    )


def leaf_probe_audit() -> tuple[int, str]:
    vertices = (0, 1, 2, 3, 4)
    base_options = (
        (),
        (("r", "r", 1),),
        (("r", "b", -2),),
        (("b", "r", 3),),
        (("b", "b", -1), ("b", "b", 4)),
    )
    checked = 0
    digest = sha256()
    for pair_index, (i, j) in enumerate(
        (u, v) for u in vertices for v in vertices if u < v
    ):
        # Vary a nontrivial original graph, including edges incident to the
        # vertices that will be leaf-forced.
        edges: list[Edge] = []
        all_pairs = [(u, v) for u in vertices for v in vertices if u < v]
        for index, (u, v) in enumerate(all_pairs):
            option = base_options[(3 * index + pair_index) % len(base_options)]
            for color_u, color_v, weight in option:
                edges.append(Edge(u, v, color_u, color_v, weight))

        target_colors = tuple(
            COLORS[(k + pair_index) % 2] for k in vertices
        )
        next_vertex = len(vertices)
        extended_vertices = list(vertices)
        context_edges: list[Edge] = []
        extended_coloring = list(target_colors)
        for k in vertices:
            if k in (i, j):
                continue
            leaf = next_vertex
            next_vertex += 1
            extended_vertices.append(leaf)
            context_edges.append(
                Edge(k, leaf, target_colors[k], COLORS[leaf % 2], 1)
            )
            extended_coloring.append(COLORS[leaf % 2])

        combined = tuple(edges + context_edges)
        amplitude = direct_amplitude(
            tuple(extended_vertices), combined, tuple(extended_coloring)
        )
        kernel = colored_pair_kernel(vertices, tuple(edges))
        expected = kernel[(i, target_colors[i], j, target_colors[j])]
        assert amplitude == expected
        digest.update(
            f"{i},{j}:{target_colors}:{amplitude}:{len(combined)}\n".encode()
        )
        checked += 1
    return checked, digest.hexdigest()


def main() -> None:
    exhaustive_count, exhaustive_digest = exhaustive_four_vertex_audit()
    selector_first, selector_second = selector_erasure_collision()
    phase_first, phase_second = phase_erasure_collision()
    leaf_count, leaf_digest = leaf_probe_audit()

    print("THM-2281 context-selected hafnian kernel exact audit")
    print(f"exhaustive direct=hafnian evaluations: {exhaustive_count}")
    print(f"exhaustive digest: {exhaustive_digest}")
    print(
        "selector-erasure response vectors (rr,rb,br,bb): "
        f"{selector_first} vs {selector_second}"
    )
    print(
        "phase-erasure response vectors (rr,rb,br,bb): "
        f"{phase_first} vs {phase_second}"
    )
    print(f"leaf-forcing entry probes: {leaf_count}")
    print(f"leaf-probe digest: {leaf_digest}")
    print("all exact assertions passed")


if __name__ == "__main__":
    main()
