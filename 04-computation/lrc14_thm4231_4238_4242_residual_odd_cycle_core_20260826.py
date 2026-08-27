#!/usr/bin/env python3
"""Exact THM-4231/4238/4242 residual odd-cycle core diagnostic."""

from __future__ import annotations

import contextlib
import hashlib
import io
import runpy
from collections import Counter

import networkx as nx


SOURCE = "04-computation/lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py"


def load_edges() -> list[tuple[int, int]]:
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        data = runpy.run_path(SOURCE)
    edges = data["aggregate_residual"]
    assert len(edges) == 181_126
    return edges


def main() -> None:
    edges = load_edges()
    graph = nx.Graph()
    graph.add_edges_from(edges)
    degrees = dict(graph.degree())
    triangles_at = nx.triangles(graph)
    triangle_count = sum(triangles_at.values()) // 3
    edge_triangle_counts = {}
    for a, b in edges:
        edge_triangle_counts[(a, b)] = len(set(graph[a]) & set(graph[b]))
    triangle_free_edges = [e for e, count in edge_triangle_counts.items() if count == 0]
    bridge_edges = [tuple(sorted(e)) for e in nx.bridges(graph)]
    leaves = sorted(v for v, degree in degrees.items() if degree == 1)
    articulation_points = sorted(nx.articulation_points(graph))
    blocks = sorted(
        (len(block), min(block), max(block))
        for block in nx.biconnected_components(graph)
    )
    components = sorted((len(c) for c in nx.connected_components(graph)), reverse=True)
    core = nx.core_number(graph)
    triangle_hist = Counter(edge_triangle_counts.values())
    triangle_counts_sorted = sorted(edge_triangle_counts.values())
    graph_digest = hashlib.sha256(
        "".join(f"{a},{b},{edge_triangle_counts[(a, b)]}\n" for a, b in sorted(edges)).encode()
    ).hexdigest()

    assert graph.number_of_nodes() == 739
    assert graph.number_of_edges() == 181_126
    assert components == [739]
    assert triangle_count == 30_912_074
    assert sum(v > 0 for v in triangles_at.values()) == 737
    assert triangle_free_edges == [(616, 756), (616, 760)]
    assert leaves == [756, 760]
    assert articulation_points == [616]
    assert bridge_edges == [(616, 756), (616, 760)]
    assert blocks == [(2, 616, 756), (2, 616, 760), (737, 1, 769)]
    assert graph_digest == "d682ec7b70b7eecc641559e26d0afa1e23b0ede3ff90c025a637e6c02e208623"
    assert max(core.values()) == 500

    print("LRC14_THM4231_4238_4242_RESIDUAL_ODD_CYCLE_CORE")
    print(f"vertices_nonisolated={graph.number_of_nodes()} edges={graph.number_of_edges()}")
    print(f"components={components}")
    print(f"degree_min={min(degrees.values())} degree_max={max(degrees.values())}")
    print(f"degree_hist_top={Counter(degrees.values()).most_common(12)}")
    print(f"triangles={triangle_count} triangle_vertices={sum(v > 0 for v in triangles_at.values())}")
    print(
        f"triangle_free_edges={len(triangle_free_edges)} "
        f"examples={triangle_free_edges[:20]}"
    )
    print(
        f"edge_triangle_min={min(edge_triangle_counts.values())} "
        f"edge_triangle_max={max(edge_triangle_counts.values())}"
    )
    print(f"bipartite={nx.is_bipartite(graph)} chordal={nx.is_chordal(graph)}")
    print(f"leaves={leaves} articulation_points={articulation_points}")
    print(f"bridges={bridge_edges}")
    print(f"biconnected_blocks={blocks}")
    print(
        "triangle_free_equals_bridges="
        f"{set(map(tuple, map(sorted, triangle_free_edges))) == set(bridge_edges)}"
    )
    print(f"edge_triangle_hist_low={sorted(triangle_hist.items())[:12]}")
    print(
        "edge_triangle_quantiles="
        f"{[(q, triangle_counts_sorted[(len(triangle_counts_sorted) - 1) * q // 100]) for q in (0, 1, 5, 25, 50, 75, 95, 99, 100)]}"
    )
    print(f"edge_triangle_sha256={graph_digest}")
    print(
        f"degeneracy={max(core.values())} "
        f"max_k_core_vertices={sum(v == max(core.values()) for v in core.values())}"
    )

    # The fixed-50 induced star/neighbourhood is a key next literal slice.
    fixed = 50
    neighbors = sorted(graph[fixed])
    induced = graph.subgraph(neighbors)
    assert len(neighbors) == 556
    assert induced.number_of_edges() == 150_861
    assert sum(nx.triangles(induced).values()) // 3 == 26_807_003
    print(
        f"fixed50_degree={len(neighbors)} neighbor_min={min(neighbors)} "
        f"neighbor_max={max(neighbors)} induced_edges={induced.number_of_edges()} "
        f"induced_triangles={sum(nx.triangles(induced).values()) // 3}"
    )
    print("verdict=proof_obligation_graph_not_safety_graph")


if __name__ == "__main__":
    main()
