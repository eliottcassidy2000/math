"""Exact bounded controls for the undirected complementary traversal game.

Run: python3 05-knowledge/results/strong_blue_games_20261020f/game.py
Uses the complete NetworkX graph atlas through seven vertices, not a sample.
All checks are explicit and survive python -O. No imported research checker.
"""
from collections import Counter, deque
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import networkx as nx


def check(test, message):
    if not test:
        raise RuntimeError(message)


def rows(graph):
    return tuple(sum(1 << w for w in graph[v]) for v in range(len(graph)))


def path_endsets(adj):
    """ends[S] is the possible final vertices of an HP of the induced graph S."""
    n = len(adj)
    ends = [0] * (1 << n)
    for v in range(n):
        ends[1 << v] = 1 << v
    for mask in range(1, 1 << n):
        vertices = ends[mask]
        while vertices:
            bit = vertices & -vertices
            vertices -= bit
            available = adj[bit.bit_length() - 1] & ~mask
            while available:
                nxt = available & -available
                available -= nxt
                ends[mask | nxt] |= nxt
    return ends


def permutation_endset(adj):
    answer = 0
    for word in permutations(range(len(adj))):
        if all(adj[u] & (1 << v) for u, v in zip(word, word[1:])):
            answer |= 1 << word[-1]
    return answer


def path_cover_number(graph):
    adj = rows(nx.convert_node_labels_to_integers(graph))
    ends = path_endsets(adj)
    covers = [0] + [len(adj) + 1] * ((1 << len(adj)) - 1)
    for mask in range(1, 1 << len(adj)):
        anchor = mask & -mask
        sub = mask
        while sub:
            if sub & anchor and ends[sub]:
                covers[mask] = min(covers[mask], 1 + covers[mask ^ sub])
            sub = (sub - 1) & mask
    return covers[-1]


def shortest_spanning_walk(graph, start=None):
    n = len(graph)
    adj = rows(graph)
    starts = range(n) if start is None else [start]
    queue = deque((1 << v, v, 0) for v in starts)
    seen = {(1 << v, v) for v in starts}
    while queue:
        mask, vertex, distance = queue.popleft()
        if mask == (1 << n) - 1:
            return distance
        for nxt in graph[vertex]:
            state = (mask | (1 << nxt), nxt)
            if state not in seen:
                seen.add(state)
                queue.append((*state, distance + 1))
    return None


def connected_hostile(n):
    """Q: K_(n-2) minus 0--1, and two new leaves attached to vertex 1."""
    a = n - 2
    graph = nx.complete_graph(a)
    graph.remove_edge(0, 1)
    graph.add_edges_from([(1, a), (1, a + 1)])
    return graph


def main():
    atlas = nx.graph_atlas_g()
    check(len(atlas) == 1253, "atlas cardinality")
    totals = Counter()
    counts = {name: Counter() for name in ["HH", "WH", "HW", "WW", "HH_prescribed", "WH_prescribed", "connected_HH"]}
    witnesses = []
    independent_graphs = join_checks = tree_checks = tree_fixed_start_checks = 0
    for index, graph in enumerate(atlas):
        n = len(graph)
        if n == 0:
            continue
        totals[n] += 1
        comp = nx.complement(graph)
        endset = path_endsets(rows(graph))[-1]
        comp_ends = path_endsets(rows(comp))[-1]
        h, hc = bool(endset), bool(comp_ends)
        w, wc = nx.is_connected(graph), nx.is_connected(comp)
        tests = {"HH": not h and not hc, "WH": not w and not hc,
                 "HW": not h and not wc, "WW": not w and not wc,
                 "HH_prescribed": not h and comp_ends != (1 << n) - 1,
                 "WH_prescribed": not w and comp_ends != (1 << n) - 1,
                 "connected_HH": w and wc and not h and not hc}
        for key, value in tests.items():
            counts[key][n] += int(value)
        if tests["connected_HH"]:
            witnesses.append({"atlas_index": index, "n": n,
                              "edges": sorted(map(list, graph.edges())),
                              "graph6": nx.to_graph6_bytes(graph, header=False).decode().strip()})
        if n <= 6:
            check(endset == permutation_endset(rows(graph)), f"permutation endpoint check {index}")
            independent_graphs += 1
        # The complement components canonically supply a complete-join partition.
        if not wc:
            blocks = list(nx.connected_components(comp))
            predicted = all(path_cover_number(graph.subgraph(block)) <= n - len(block) + 1 for block in blocks)
            check(predicted == h, f"complete-join path-cover criterion {index}")
            join_checks += 1
        if w and nx.is_tree(graph):
            diameter = nx.diameter(graph)
            check(shortest_spanning_walk(graph) == 2 * (n - 1) - diameter,
                  f"tree free spanning walk {index}")
            tree_checks += 1
            for start in range(n):
                check(shortest_spanning_walk(graph, start) == 2 * (n - 1) - nx.eccentricity(graph, start),
                      f"tree fixed start {index}:{start}")
                tree_fixed_start_checks += 1
    check(dict(totals) == {1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044}, "atlas order counts")
    for key in ["HH", "WH", "HW"]:
        check(all(counts[key][n] == 0 for n in range(1, 4)), f"small threshold {key}")
        check(all(counts[key][n] > 0 for n in range(4, 8)), f"positive threshold {key}")
    check(not any(counts["WW"].values()), "no disconnected complementary pair")
    check([item["atlas_index"] for item in witnesses] == [615, 743], "connected minimal witness pair")
    check(nx.is_isomorphic(nx.complement(atlas[615]), atlas[743]), "witness complement pair")
    check(nx.is_isomorphic(connected_hostile(7), atlas[743]), "all-order family meets first witness")
    for n in range(4, 31):
        graph = nx.disjoint_union(nx.complete_graph(n - 1), nx.empty_graph(1))
        comp = nx.complement(graph)
        check(not nx.is_connected(graph), "basic witness disconnected")
        check(sum(comp.degree(v) == 1 for v in comp) == n - 1, "basic complement star")
    for n in range(7, 31):
        graph = connected_hostile(n)
        comp = nx.complement(graph)
        check(nx.is_connected(graph) and nx.is_connected(comp), "family connectivity")
        minus_cut = graph.copy()
        minus_cut.remove_node(1)
        check(nx.number_connected_components(minus_cut) == 3, "family one-cut obstruction")
        minus_cut = comp.copy()
        minus_cut.remove_nodes_from([n - 2, n - 1])
        check(nx.number_connected_components(minus_cut) == n - 3 > 3, "family two-cut obstruction")
    receipt = {
        "status": "FINITE-EXACT controls for separately proved all-order claims",
        "universe": "All 1253 NetworkX/Read-Wilson atlas graphs on 0..7 vertices; omit n=0 in the game",
        "networkx_version": nx.__version__,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "order_counts": dict(totals),
        "winning_isomorphism_classes_by_rule": {key: dict(value) for key, value in counts.items()},
        "connected_HH_witnesses": witnesses,
        "independent_permutation_endpoint_checks_through_6": independent_graphs,
        "complete_join_path_cover_checks": join_checks,
        "tree_walk_BFS_checks": tree_checks,
        "tree_fixed_start_walk_BFS_checks": tree_fixed_start_checks,
        "basic_family_orders_checked_structurally": [4, 30],
        "connected_family_orders_checked_structurally": [7, 30],
        "scope": "Unlabelled atlas counts, not counts of labelled edge colorings. P2 chooses the start in the four main rules."
    }
    target = Path(__file__).with_suffix(".json")
    target.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
