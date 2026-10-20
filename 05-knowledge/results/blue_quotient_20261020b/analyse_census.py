"""Analyse the exact n9 universe, complete blue graph, and black tree.

The complete black graph is NOT enumerated. Its connected components are
determined by an authenticated spanning tree plus exact pure-blue detection.
Run from any directory: python3 analyse_census.py
"""
from collections import Counter, deque
import csv
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import networkx as nx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def table(name):
    path = HERE / name
    stream = path.open(encoding="utf-8") if path.exists() else gzip.open(
        str(path) + ".gz", "rt", encoding="utf-8")
    with stream:
        return [{k: int(v) for k, v in row.items()}
                for row in csv.DictReader(stream, delimiter="\t")]


def rows(code, n=9):
    answer = [0] * n
    bit = n * (n - 1) // 2 - 1
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                answer[i] |= 1 << j
            else:
                answer[j] |= 1 << i
            bit -= 1
    return answer


def unpack(code):
    return tuple(code >> (4 * i) & 15 for i in range(9))


def first_paths(adj, count=4):
    answer = []

    def search(path, used):
        if len(answer) >= count:
            return
        if len(path) == 9:
            answer.append(tuple(path))
            return
        possible = ((1 << 9) - 1) if not path else adj[path[-1]] & ~used
        for v in range(9):
            if possible >> v & 1:
                search(path + [v], used | (1 << v))
    search([], 0)
    return answer


def distances(adjacency, start):
    d = {start: 0}
    todo = deque([start])
    while todo:
        v = todo.popleft()
        for w in adjacency[v]:
            if w not in d:
                d[w] = d[v] + 1
                todo.append(w)
    return d


def main():
    for name in ("n9_summary.txt", "n9_independent_audit.txt"):
        path = HERE / name
        path.write_text(path.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    classes = {c["code"]: c for c in table("n9_classes.tsv")}
    blue = table("n9_blue.tsv")
    tree = table("n9_black_tree.tsv")
    nodes = {c for c, data in classes.items() if c <= data["converse_code"]}
    sc = {c for c in nodes if classes[c]["converse_code"] == c}
    pb = {c for c in sc if classes[c]["SC_H"] // classes[c]["Aut"] == classes[c]["blue_mass"]}
    active = nodes - pb
    merged = lambda c: min(c, classes[c]["converse_code"])
    adjacency = {c: set() for c in active}
    for edge in tree:
        u, v = merged(edge["source_code"]), merged(edge["target_code"])
        adjacency[u].add(v)
        adjacency[v].add(u)
    need(len(tree) == len(active) - 1, "tree size")
    first = min(active)
    d = distances(adjacency, first)
    need(len(d) == len(active), "tree spans all active nodes")
    end_a = max(d, key=d.get)
    d_a = distances(adjacency, end_a)
    end_b = max(d_a, key=d_a.get)
    diameter = d_a[end_b]
    d_b = distances(adjacency, end_b)
    radius = min(max(d_a[v], d_b[v]) for v in active)
    centres = sorted(v for v in active if max(d_a[v], d_b[v]) == radius)
    hub = max(active, key=lambda v: (len(adjacency[v]), -v))
    hub_distances = distances(adjacency, hub)
    blue_graph = nx.Graph()
    blue_graph.add_nodes_from(sc)
    blue_mass = Counter()
    quotient = nx.MultiGraph()
    quotient.add_nodes_from(["active", *sorted(pb)])
    pb_neighbors = {c: [] for c in pb}
    blue_multiplicities = {}
    for edge in blue:
        a, b, multiplicity = (edge[k] for k in ("source_code", "target_code", "multiplicity"))
        blue_graph.add_edge(a, b)
        blue_mass[a] += multiplicity
        blue_mass[b] += multiplicity
        blue_multiplicities[tuple(sorted((a, b)))] = multiplicity
        qa, qb = (a if a in pb else "active"), (b if b in pb else "active")
        for _ in range(multiplicity):
            quotient.add_edge(qa, qb, original=(a, b))
        for v, w in ((a, b), (b, a)):
            if v in pb:
                need(w not in pb, "pure-blue neighbour not mixed")
                pb_neighbors[v].append((w, multiplicity))
    need(all(blue_mass[c] == classes[c]["blue_mass"] for c in sc), "blue degree law")
    actual_bridges = []
    for a, b in nx.bridges(quotient):
        data = next(iter(quotient[a][b].values()))
        actual_bridges.append(tuple(sorted(data["original"])))
    blue_only_bridges = [tuple(sorted(e)) for e in nx.bridges(blue_graph)
                         if blue_multiplicities[tuple(sorted(e))] == 1]
    need(len(actual_bridges) == 5, "n9 full bridge count")
    need(all(a in pb or b in pb for a, b in actual_bridges), "bridge is not pure-blue incident")
    for a, b in actual_bridges:
        leaf = a if a in pb else b
        need(classes[leaf]["blue_mass"] == 1, "bridge at trivalent PB")
    rare_edges = tree[len(active) - 237:]
    need(len(rare_edges) == 236, "second-round witness edges")
    rank_hist, rare_classes, rare_details = Counter(), {}, []
    for edge in rare_edges:
        source = edge["source_code"]
        adj = rows(source)
        path = unpack(edge["path_order"])
        paths = first_paths(adj)
        need(path in paths, "rare witness not among first four paths")
        rank = paths.index(path) + 1
        need(rank in (2, 3, 4), "unexpected second round rank")
        rank_hist[rank] += 1
        rare_classes[source] = max(rare_classes.get(source, 0), rank)
        graph = nx.DiGraph()
        graph.add_nodes_from(range(9))
        graph.add_edges_from((i, j) for i in range(9) for j in range(9) if adj[i] >> j & 1)
        components = sorted((len(s) for s in nx.strongly_connected_components(graph)), reverse=True)
        score = [x.bit_count() for x in adj]
        rare_details.append(dict(source_code=source, target_code=edge["target_code"],
                                 path=path, rank=rank, Aut=classes[source]["Aut"],
                                 self_converse=source in sc, scores=sorted(score),
                                 scc_sizes=components,
                                 triangle_change=score[path[0]] - score[path[-1]] - 1))
    delta_hist = Counter()
    for edge in tree:
        score = [x.bit_count() for x in rows(edge["source_code"])]
        path = unpack(edge["path_order"])
        delta_hist[score[path[0]] - score[path[-1]] - 1] += 1
    report = dict(
        universe=dict(n=9, ordinary_classes=len(classes), SC_classes=len(sc),
                      merged_nodes=len(nodes), all_blue_tilings=65536, all_blue_lines=32768,
                      complete_black_edge_census=False),
        connectivity=dict(black_active_nodes=len(active), authenticated_black_tree_edges=len(tree),
                          black_components=1 + len(pb), full_components=1,
                          blue_touched_components=nx.number_connected_components(blue_graph)),
        quotient=dict(nodes=1 + len(pb), loop_lines=quotient.number_of_edges("active", "active"),
                      pure_blue_incident_lines=sum(classes[c]["blue_mass"] for c in pb),
                      full_metagraph_bridges=sorted(actual_bridges), non_leaf_bridges=[],
                      blue_only_bridges=sorted(blue_only_bridges),
                      blue_only_nonleaf_bridges=sorted(
                          e for e in blue_only_bridges if not any(x in pb for x in e)),
                      pure_blue_neighbors={str(c): sorted(pb_neighbors[c]) for c in sorted(pb)}),
        tree_geometry=dict(diameter=diameter, diameter_endpoints=[end_a, end_b],
                           radius=radius, centres=centres, max_degree=len(adjacency[hub]),
                           max_degree_vertex=hub, max_distance_from_hub=max(hub_distances.values()),
                           mean_distance_from_hub=str(Fraction(sum(hub_distances.values()), len(active))),
                           degree_histogram=dict(sorted(Counter(map(len, adjacency.values())).items()))),
        second_round=dict(edge_count=len(rare_edges), source_classes=len(rare_classes),
                          path_rank_histogram=dict(sorted(rank_hist.items())),
                          self_converse_source_classes=sum(c in sc for c in rare_classes),
                          max_rank_by_source_histogram=dict(sorted(Counter(rare_classes.values()).items())),
                          detail_file="n9_second_round.json"),
        tree_triangle_change_histogram=dict(sorted(delta_hist.items())),
    )
    (HERE / "n9_analysis.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    (HERE / "n9_second_round.json").write_text(json.dumps(rare_details, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: report[k] for k in ("universe", "connectivity", "tree_geometry", "second_round")}, indent=2))
    print("Full graph bridges:", actual_bridges)
    print("Blue graph bridges:", len(blue_only_bridges))
    print("ALL ANALYSIS CHECKS PASSED")
    for name in ("n9_classes.tsv", "n9_blue.tsv", "n9_black_tree.tsv"):
        raw = HERE / name
        if raw.exists():
            raw_bytes = raw.read_bytes()
            (HERE / (name + ".gz")).write_bytes(gzip.compress(raw_bytes, mtime=0))
    hashed = ["census.cpp", "audit.cpp", "analyse_census.py", "replay_census.py",
              "census.md", ".gitignore", "n9_analysis.json",
              "n9_second_round.json", "n9_summary.txt", "n9_independent_audit.txt",
              "n9_classes.tsv.gz", "n9_blue.tsv.gz", "n9_black_tree.tsv.gz"]
    manifest = {
        "owned_sha256": {
            name: hashlib.sha256((HERE / name).read_bytes() if name.endswith(".gz")
                                 else (HERE / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
            for name in hashed},
        "inherited_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
            for name in (
                "04-computation/experiments/tiling_line_metagraph_census_20261020.cpp",
                "04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv",
            )
        },
        "compiler": "g++.exe (GCC) 15.2.0",
        "flags": "-O3 -std=c++17",
        "encoding": "Text hashes use canonical LF; raw archived Windows TSVs retain CRLF; deterministic gzip mtime=0",
    }
    (HERE / "census_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
