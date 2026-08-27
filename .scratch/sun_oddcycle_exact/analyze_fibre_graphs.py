#!/usr/bin/env python3
"""Exact graph/hypergraph observables for the three THM-4026 hostile targets."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import sys
from collections import Counter, deque
from math import comb
from pathlib import Path


ROLES = ("w", "x", "y", "z")


def load(path: Path) -> list[tuple[int, int, int, int]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [tuple(int(row[k]) for k in ROLES) for row in csv.DictReader(handle, delimiter="\t")]
    assert len(rows) == len(set(rows))
    return rows


def shared_coordinates(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(i for i in range(4) if a[i] == b[i])


def graph(reps: list[tuple[int, ...]], predicate) -> list[set[int]]:
    adj = [set() for _ in reps]
    for i, j in itertools.combinations(range(len(reps)), 2):
        shared = shared_coordinates(reps[i], reps[j])
        if predicate(shared):
            adj[i].add(j)
            adj[j].add(i)
    return adj


def components(adj: list[set[int]]) -> list[list[int]]:
    out = []
    unseen = set(range(len(adj)))
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        todo = [root]
        comp = []
        while todo:
            v = todo.pop()
            comp.append(v)
            for u in adj[v]:
                if u in unseen:
                    unseen.remove(u)
                    todo.append(u)
        out.append(sorted(comp))
    return out


def bipartite_and_witness(adj: list[set[int]]) -> tuple[bool, list[int] | None]:
    """Return bipartiteness and one shortest odd closed-walk witness if false."""
    color: dict[int, int] = {}
    for root in range(len(adj)):
        if root in color:
            continue
        color[root] = 0
        parent = {root: -1}
        todo = deque([root])
        while todo:
            v = todo.popleft()
            for u in sorted(adj[v]):
                if u not in color:
                    color[u] = color[v] ^ 1
                    parent[u] = v
                    todo.append(u)
                elif color[u] == color[v]:
                    av = []
                    x = v
                    while x != -1:
                        av.append(x)
                        x = parent[x]
                    au = []
                    x = u
                    while x != -1:
                        au.append(x)
                        x = parent[x]
                    pos = {x: i for i, x in enumerate(av)}
                    meet = next(x for x in au if x in pos)
                    left = av[: pos[meet] + 1]
                    right = au[: au.index(meet)]
                    return False, left + list(reversed(right)) + [v]
    return True, None


def triangles(adj: list[set[int]]) -> list[tuple[int, int, int]]:
    out = []
    for i in range(len(adj)):
        for j in sorted(u for u in adj[i] if u > i):
            for k in sorted(u for u in adj[i] & adj[j] if u > j):
                out.append((i, j, k))
    return out


def graph_summary(reps: list[tuple[int, ...]], adj: list[set[int]]) -> dict:
    tris = triangles(adj)
    bip, odd = bipartite_and_witness(adj)
    comps = components(adj)
    edge_count = sum(map(len, adj)) // 2
    return {
        "vertices": len(reps),
        "edges": edge_count,
        "components": len(comps),
        "nontrivial_components": sum(len(c) > 1 for c in comps),
        "largest_component": max(map(len, comps), default=0),
        "triangles": len(tris),
        "bipartite": bip,
        "odd_cycle_vertex_indices": odd,
        "odd_cycle_representations": None if odd is None else [reps[i] for i in odd],
    }


def overlap_triangle_types(reps: list[tuple[int, ...]], tris: list[tuple[int, int, int]]) -> Counter:
    kinds = Counter()
    for i, j, k in tris:
        labels = [
            set(shared_coordinates(reps[i], reps[j])),
            set(shared_coordinates(reps[j], reps[k])),
            set(shared_coordinates(reps[k], reps[i])),
        ]
        common = labels[0] & labels[1] & labels[2]
        if common:
            kinds["common_coordinate_clique"] += 1
        elif any(not x for x in labels):
            raise AssertionError("non-edge in overlap triangle")
        else:
            # A Berge triangle exists iff three distinct role-coordinate nodes
            # can be chosen, one from each pairwise intersection.
            systems = set(itertools.product(*[sorted(x) for x in labels]))
            if any(len(set(choice)) == 3 for choice in systems):
                kinds["berge_triangle"] += 1
            else:
                kinds["overlap_triangle_no_distinct_representatives"] += 1
    return kinds


def role_multiplicities(reps: list[tuple[int, ...]]) -> dict:
    result = {}
    for role, idx in zip(ROLES, range(4)):
        counts = Counter(rep[idx] for rep in reps)
        result[role] = {
            "distinct": len(counts),
            "maximum_fibre": max(counts.values(), default=0),
            "collision_pairs": sum(v * (v - 1) // 2 for v in counts.values()),
        }
    return result


def pair_split_graph(reps: list[tuple[int, ...]]) -> dict:
    """Statistics for the intrinsic low/high pair-incidence graph.

    Low vertices are (w,x), high vertices are (y,z), and a representation is
    an edge.  Equal low sums index the complete-bipartite components.
    """
    blocks: dict[int, tuple[set[tuple[int, int]], set[tuple[int, int]]]] = {}
    for w, x, y, z in reps:
        low_sum = comb(w, 2) + comb(x, 4)
        low, high = blocks.setdefault(low_sum, (set(), set()))
        low.add((w, x))
        high.add((y, z))
    block_sizes = sorted((len(low), len(high)) for low, high in blocks.values())
    expected_edges = sum(a * b for a, b in block_sizes)
    assert expected_edges == len(reps)
    return {
        "components": len(block_sizes),
        "block_size_histogram": {
            f"K_{a},{b}": count
            for (a, b), count in sorted(Counter(block_sizes).items())
        },
        "edges": expected_edges,
        "cycle_rank": sum((a - 1) * (b - 1) for a, b in block_sizes),
        "four_cycles": sum(comb(a, 2) * comb(b, 2) for a, b in block_sizes),
        "odd_cycles": 0,
        "structural_reason": "disjoint union of complete bipartite graphs K_(L_s,H_(n-s))",
    }


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: analyze_fibre_graphs DIRECTORY", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    inputs = {
        "N-1": root / "N_minus_1.tsv",
        "N": root / "N.tsv",
        "N+1": root / "N_plus_1.tsv",
    }
    report = {"definition": {
        "overlap_graph": "representations adjacent iff they share at least one role index",
        "two_role_exchange_graph": "representations adjacent iff they share exactly two role indices",
        "incidence_graph": "bipartite graph between representations and role-index vertices",
    }, "targets": {}}
    for name, path in inputs.items():
        reps = load(path)
        overlap = graph(reps, lambda shared: len(shared) >= 1)
        exchange = graph(reps, lambda shared: len(shared) == 2)
        overlap_tris = triangles(overlap)
        report["targets"][name] = {
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "representations": len(reps),
            "role_multiplicities": role_multiplicities(reps),
            "overlap_graph": graph_summary(reps, overlap),
            "overlap_triangle_types": dict(overlap_triangle_types(reps, overlap_tris)),
            "two_role_exchange_graph": graph_summary(reps, exchange),
            "low_high_pair_incidence_graph": pair_split_graph(reps),
            "incidence_graph_bipartite_by_construction": True,
        }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
