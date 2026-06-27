#!/usr/bin/env python3
"""Small graph audit for the LRC14 Desargues/median forum lens.

The point is not to reduce LRC14 to graph theory.  It is to make the
analogy precise enough to be useful:

* median graphs model quotients whose sidecar coordinates behave like
  compatible hyperplanes, so three proof routes have a unique consensus state;
* the Desargues graph is a compact bipartite incidence graph that still fails
  medianity, so incidence/bipartiteness alone is not a safe quotient.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations_with_replacement

import networkx as nx


def girth(G: nx.Graph) -> int | None:
    best = None
    for source in G:
        dist = {source: 0}
        parent = {source: None}
        q = deque([source])
        while q:
            u = q.popleft()
            for v in G[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    q.append(v)
                elif parent[u] != v and parent[v] != u:
                    length = dist[u] + dist[v] + 1
                    best = length if best is None else min(best, length)
    return best


def interval_nodes(dists: dict, u, v) -> frozenset:
    duv = dists[u][v]
    return frozenset(x for x in dists if dists[u][x] + dists[x][v] == duv)


def median_failures(G: nx.Graph):
    dists = dict(nx.all_pairs_shortest_path_length(G))
    failures = []
    multiplicities = Counter()
    nodes = list(G.nodes())
    for a, b, c in combinations_with_replacement(nodes, 3):
        common = (
            interval_nodes(dists, a, b)
            & interval_nodes(dists, b, c)
            & interval_nodes(dists, c, a)
        )
        multiplicities[len(common)] += 1
        if len(common) != 1:
            failures.append((a, b, c, tuple(sorted(common, key=str))))
    return failures, multiplicities


def cycle_length_histogram(G: nx.Graph, max_len: int = 16) -> Counter:
    """Count simple undirected cycles up to max_len by canonical node word."""
    seen = set()
    hist = Counter()
    nodes = list(G.nodes())
    node_order = {v: i for i, v in enumerate(nodes)}

    def canon(cycle):
        idxs = [node_order[x] for x in cycle]
        rotations = []
        n = len(idxs)
        for word in (idxs, list(reversed(idxs))):
            for i in range(n):
                rotations.append(tuple(word[i:] + word[:i]))
        return min(rotations)

    def dfs(start, cur, path):
        if len(path) > max_len:
            return
        for nxt in G[cur]:
            if nxt == start and len(path) >= 3:
                key = canon(path)
                if key not in seen:
                    seen.add(key)
                    hist[len(path)] += 1
            elif nxt not in path and node_order[nxt] >= node_order[start]:
                dfs(start, nxt, path + [nxt])

    for start in nodes:
        for nbr in G[start]:
            if node_order[nbr] >= node_order[start]:
                dfs(start, nbr, [start, nbr])
    return hist


def theta_classes(G: nx.Graph):
    """Return Djokovic-Winkler relation components for a bipartite graph."""
    d = dict(nx.all_pairs_shortest_path_length(G))
    edges = [tuple(sorted(e, key=str)) for e in G.edges()]
    parent = {e: e for e in edges}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for e in edges:
        a, b = e
        for f in edges:
            x, y = f
            if d[a][x] + d[b][y] != d[a][y] + d[b][x]:
                union(e, f)

    classes = defaultdict(list)
    for e in edges:
        classes[find(e)].append(e)
    return list(classes.values())


def summarize(name: str, G: nx.Graph) -> dict:
    failures, med_mult = median_failures(G)
    cycle_hist = cycle_length_histogram(G)
    theta = theta_classes(G) if nx.is_bipartite(G) else []
    return {
        "name": name,
        "n": G.number_of_nodes(),
        "m": G.number_of_edges(),
        "degree_hist": dict(sorted(Counter(dict(G.degree()).values()).items())),
        "bipartite": nx.is_bipartite(G),
        "diameter": nx.diameter(G),
        "girth": girth(G),
        "cycle_lengths_le_16": dict(sorted(cycle_hist.items())),
        "median": not failures,
        "median_intersection_size_hist": dict(sorted(med_mult.items())),
        "median_failure_count": len(failures),
        "sample_failures": failures[:5],
        "theta_class_sizes": sorted([len(c) for c in theta]),
        "theta_class_count": len(theta),
    }


def print_summary(summary: dict) -> None:
    print(f"== {summary['name']} ==")
    for key in [
        "n",
        "m",
        "degree_hist",
        "bipartite",
        "diameter",
        "girth",
        "cycle_lengths_le_16",
        "median",
        "median_intersection_size_hist",
        "median_failure_count",
        "theta_class_count",
        "theta_class_sizes",
    ]:
        print(f"{key}: {summary[key]}")
    if summary["sample_failures"]:
        print("sample_median_failures:")
        for failure in summary["sample_failures"]:
            print(f"  {failure}")
    print()


def main() -> None:
    graphs = [
        ("Desargues graph", nx.desargues_graph()),
        ("Q4 hypercube", nx.hypercube_graph(4)),
        ("4x4 grid", nx.grid_2d_graph(4, 4)),
        ("C6 cycle", nx.cycle_graph(6)),
    ]
    print("LRC14 Desargues/median graph lens audit")
    print("Median test: every triple has exactly one common interval node.")
    print("Theta classes are Djokovic-Winkler components, used here as a sidecar sketch.")
    print()
    for name, graph in graphs:
        print_summary(summarize(name, graph))

    print("Readout")
    print("- Q4 and the 4x4 grid pass the median triple test: compatible sidecar")
    print("  hyperplanes give unique three-route consensus states.")
    print("- C6 is bipartite and small but fails medianity by triples with no")
    print("  candidate center; even-cycle structure is not enough.")
    print("- The Desargues graph is cubic, bipartite, girth 6, and has many median")
    print("  failures.  As an LRC analogy it is a compact incidence/duality warning,")
    print("  not a proof carrier by itself.")


if __name__ == "__main__":
    main()
