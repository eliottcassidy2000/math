"""Exact bridge and abstract wiring audit for the tournament line metagraph.

Universe: all committed full class/line TSV censuses at n=3,...,7.
The n=8 line TSV is not committed, so this script does not claim to replay it.
Run: python3 04-computation/boundary_support_20261020_metagraph.py

No assertions: the audit is unchanged under python -O.
"""
from collections import Counter, defaultdict
import csv
import importlib.util
from itertools import combinations
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation" / "experiments"


def need(test, label):
    if not test:
        raise RuntimeError(label)


def load(n):
    prefix = SOURCE / f"tiling_line_metagraph_census_20261020"
    with open(str(prefix) + f"_classes_n{n}.tsv", encoding="utf-8") as stream:
        classes = {}
        for row in csv.DictReader(stream, delimiter="\t"):
            item = dict(row)
            for key in ("class", "H", "Aut", "tc", "bm", "SC", "sigmaClass"):
                item[key] = int(item[key])
            classes[item["class"]] = item
    with open(str(prefix) + f"_lines_n{n}.tsv", encoding="utf-8") as stream:
        lines = [tuple(int(row[key]) for key in
                       ("classA", "classB", "blue", "black"))
                 for row in csv.DictReader(stream, delimiter="\t")]
    return classes, lines


def degree_polynomial(n):
    """Conditional interpretation: degree counts, loops counted twice."""
    with open(SOURCE / f"tiling_line_metagraph_census_20261020_classes_n{n}.tsv",
              encoding="utf-8") as stream:
        classes = list(csv.DictReader(stream, delimiter="\t"))
    coefficients = Counter()
    sc_count = sum(int(row["SC"]) for row in classes)
    for row in classes:
        if int(row["class"]) > int(row["sigmaClass"]):
            continue
        degree = int(row["tc"]) * (1 if int(row["SC"]) else 2)
        coefficients[degree] += 1
    need(sum(coefficients.values()) == (len(classes) + sc_count) // 2,
         "degree polynomial value at 1")
    need(sum(k * value for k, value in coefficients.items())
         == 2 ** ((n - 1) * (n - 2) // 2), "degree polynomial derivative at 1")
    need(sum((-1) ** k * value for k, value in coefficients.items())
         == (len(classes) - 3 * sc_count) // 2, "degree polynomial at -1")
    need(all(k % 4 for k in coefficients), "no degrees divisible by four")
    print(f"D_{n}: M_{n}={max(coefficients)}, coefficients={dict(sorted(coefficients.items()))}")
    maximal = max(coefficients)
    maximizers = [(int(row["class"]), int(row["H"]), int(row["Aut"]), int(row["SC"]))
                  for row in classes
                  if int(row["class"]) <= int(row["sigmaClass"])
                  and int(row["tc"]) * (1 if int(row["SC"]) else 2) == maximal]
    print(f"  H_max={max(int(row['H']) for row in classes)}; "
          f"M-maximizers (class,H,Aut,SC)={maximizers}")


def merge(classes, lines):
    mid = {c: min(c, row["sigmaClass"]) for c, row in classes.items()}
    blue, black = Counter(), Counter()
    for a, b, bl, bk in lines:
        pair = tuple(sorted((mid[a], mid[b])))
        blue[pair] += bl
        black[pair] += bk
    return set(mid.values()), +blue, +black


def full_graph(nodes, blue, black):
    graph = nx.MultiGraph()
    graph.add_nodes_from(nodes)
    for colour, edges in (("blue", blue), ("black", black)):
        for (a, b), multiplicity in edges.items():
            for _ in range(multiplicity):
                graph.add_edge(a, b, colour=colour)
    return graph


def bridges_by_deletion(graph):
    """Independent component-count test, only possible singleton edge pairs."""
    base = nx.number_connected_components(graph)
    answer = set()
    for a, b in nx.Graph(graph).edges:
        if a == b or graph.number_of_edges(a, b) != 1:
            continue
        key, data = next(iter(graph[a][b].items()))
        saved = dict(data)
        graph.remove_edge(a, b, key)
        if nx.number_connected_components(graph) > base:
            answer.add(tuple(sorted((a, b))))
        graph.add_edge(a, b, key=key, **saved)
    return answer


def black_quotient(nodes, blue, black):
    black_graph = nx.Graph()
    black_graph.add_nodes_from(nodes)
    black_graph.add_edges_from(black)
    components = list(nx.connected_components(black_graph))
    owner = {v: i for i, part in enumerate(components) for v in part}
    quotient = nx.MultiGraph()
    quotient.add_nodes_from(range(len(components)))
    for (a, b), multiplicity in blue.items():
        for _ in range(multiplicity):
            quotient.add_edge(owner[a], owner[b], original=(a, b))
    lifted = set()
    for a, b in nx.bridges(quotient):
        data = next(iter(quotient[a][b].values()))
        lifted.add(tuple(sorted(data["original"])))
    return components, quotient, lifted


def involution_realization(nodes, sc_nodes, blue, black):
    """Build commuting sigma,tau on flags; black inputs are actual multiplicities."""
    need(set(sc_nodes) <= set(nodes), "self-converse node domain")
    states, sigma, tau, labels, colours = [], {}, {}, {}, {}

    def put(kind, endpoints):
        start = len(states)
        if kind == "blue":
            need(all(v in sc_nodes for v in endpoints), "blue endpoint type")
            ids = (start, start + 1)
            states.extend(ids)
            for side in (0, 1):
                x = ids[side]
                sigma[x] = x
                tau[x] = ids[1 - side]
                labels[x] = (endpoints[side], 0)
                colours[x] = kind
        else:
            ids = {(side, sheet): start + 2 * side + sheet
                   for side in (0, 1) for sheet in (0, 1)}
            states.extend(ids.values())
            for (side, sheet), x in ids.items():
                sigma[x] = ids[(side, 1 - sheet)]
                tau[x] = ids[(1 - side, sheet)]
                v = endpoints[side]
                labels[x] = (v, 0 if v in sc_nodes else sheet)
                colours[x] = kind

    for pair, multiplicity in blue.items():
        for _ in range(multiplicity):
            put("blue", pair)
    for pair, multiplicity in black.items():
        need(multiplicity % 2 == 0, "even black multiplicity")
        for _ in range(multiplicity // 2):
            put("black", pair)
    for x in states:
        need(sigma[sigma[x]] == x and tau[tau[x]] == x,
             "involution laws")
        need(sigma[tau[x]] == tau[sigma[x]], "commutation")
        need(tau[x] != x and sigma[tau[x]] != x, "free tau and sigma*tau")
        need((sigma[x] == x) == (colours[x] == "blue"), "blue fixed locus")
        v, sheet = labels[x]
        target = (v, 0 if v in sc_nodes else 1 - sheet)
        need(labels[sigma[x]] == target, "converse label law")
    rebuilt_blue, rebuilt_black = Counter(), Counter()
    for x in states:
        if x < tau[x]:
            pair = tuple(sorted((labels[x][0], labels[tau[x]][0])))
            (rebuilt_blue if colours[x] == "blue" else rebuilt_black)[pair] += 1
    need(+rebuilt_blue == +blue and +rebuilt_black == +black,
         "exact reconstructed wiring")
    return Counter(labels.values()), len(states)


def census(n):
    classes, lines = load(n)
    nodes, blue, black = merge(classes, lines)
    need(all(v % 2 == 0 for v in black.values()), "sigma-twin multiplicities")
    graph = full_graph(nodes, blue, black)
    need(nx.is_connected(graph), "actual census connectivity")
    connected_parts = list(nx.connected_components(graph))
    bipartite_parts = sum(nx.is_bipartite(graph.subgraph(part))
                          for part in connected_parts)
    defect_rank = len(graph) - len(connected_parts) - bipartite_parts
    direct = {tuple(sorted(edge)) for edge in nx.bridges(graph)}
    independently = bridges_by_deletion(graph)
    need(direct == independently, "bridge deletion cross-check")
    components, quotient, lifted = black_quotient(nodes, blue, black)
    need(direct == lifted, "black contraction bridge equivalence")
    sc_nodes = {v for v in nodes if classes[v]["SC"]}
    fibres, total = involution_realization(nodes, sc_nodes, blue, black)
    for v in nodes:
        if v in sc_nodes:
            need(fibres[v, 0] == classes[v]["tc"], "SC fibre size")
        else:
            need(fibres[v, 0] == fibres[v, 1] == classes[v]["tc"], "NS fibre sizes")
    need(total == 2 ** ((n - 1) * (n - 2) // 2), "total tiling mass")
    nonleaf = []
    for a, b in sorted(direct):
        need(blue[a, b] == 1 and black[a, b] == 0, "bridge singleton blue")
        if graph.degree(a) != 1 and graph.degree(b) != 1:
            nonleaf.append((a, b))
        cut = graph.copy()
        cut.remove_edge(a, b)
        for part in nx.connected_components(cut):
            need(len(part & sc_nodes) % 2 == 1, "odd SC bridge shore")
    print(f"n={n}: nodes={len(nodes)}, line_instances={graph.number_of_edges()}, "
          f"black_components={len(components)}, bridges={len(direct)}, "
          f"nonleaf_bridges={len(nonleaf)}")
    print(f"  fullsupport_components={len(connected_parts)}, "
          f"bipartite_components={bipartite_parts}, "
          f"proved_defect_rank={defect_rank}")
    for a, b in sorted(direct):
        print("  bridge", (a, b), "degrees", (graph.degree(a), graph.degree(b)),
              "H", (classes[a]["H"], classes[b]["H"]),
              "fibres", (classes[a]["tc"], classes[b]["tc"]))
    return classes, graph


def independent_small_census(n):
    """Recompute tiling classes with the separately authored Python producer."""
    name = ROOT / "05-knowledge/results/metagraph_growth_20261020/enumeration.py"
    spec = importlib.util.spec_from_file_location("inherited_enumeration", name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    m = (n - 1) * (n - 2) // 2
    canonical = {}
    blue, black = Counter(), Counter()
    for t in range(1 << m):
        rows = module.tiling(n, t)
        c = module.canonical(rows)
        dual = module.canonical(module.converse(rows))
        canonical[t] = min(c, dual)
    for t in range(1 << m):
        mate = t ^ ((1 << m) - 1)
        if t < mate:
            pair = tuple(sorted((canonical[t], canonical[mate])))
            (blue if module.sigma(n, t) == t else black)[pair] += 1
    graph = full_graph(set(canonical.values()), blue, black)
    if n <= 5:
        from sympy import Matrix
        labels = sorted(set(canonical.values()))
        index = {c: i for i, c in enumerate(labels)}
        indicator = Matrix([[int(index[canonical[t]] == j) for j in range(len(labels))]
                            for t in range(1 << m)])
        translated = Matrix([[int(index[canonical[t ^ ((1 << m) - 1)]] == j)
                              for j in range(len(labels))] for t in range(1 << m)])
        degree = indicator.T * indicator
        incidence = indicator.T * translated
        need(incidence == incidence.T, "symmetric contraction kernel")
        residual = translated - indicator * degree.inv() * incidence
        defect = degree - incidence * degree.inv() * incidence
        need(residual.T * residual == defect, "exact quotient-loss Gram identity")
        parts = list(nx.connected_components(graph))
        expected_rank = len(graph) - len(parts) - sum(
            nx.is_bipartite(graph.subgraph(part)) for part in parts)
        need(defect.rank() == expected_rank, "rational rank versus even-walk rank theorem")
        print(f"  kernel n={n}: nodes={len(labels)}, "
              f"continuation_defect_rank={defect.rank()}")
    return (len(graph), graph.number_of_edges(), len(list(nx.bridges(graph))),
            Counter(dict(graph.degree()).values()))


def hostile():
    from sympy import Matrix, diag
    disconnected = nx.MultiGraph()
    disconnected.add_edges_from([(0, 1), (2, 2), (3, 4), (4, 5), (3, 5)])
    adjacency = Matrix([[disconnected.number_of_edges(i, j) * (2 if i == j else 1)
                         for j in range(6)] for i in range(6)])
    degree = diag(*(disconnected.degree(i) for i in range(6)))
    parts = list(nx.connected_components(disconnected))
    bipartite = sum(nx.is_bipartite(disconnected.subgraph(part)) for part in parts)
    need(len(parts) == 3 and bipartite == 1, "disconnected rank control types")
    need((degree - adjacency * degree.inv() * adjacency).rank() == 6 - 3 - 1,
         "disconnected mixed-bipartiteness rank control")
    print("Rank hostile/control: disjoint K2 + looped singleton + triangle: "
          "r=6,c=3,b=1, exact rational defect rank=2.")

    # All six blue degrees are odd. The central blue bridge has degrees 3,3.
    # Parallel black pairs inside the two shores are not bridges.
    nodes = set(range(6))
    blue = Counter({(0, 2): 1, (1, 2): 1, (2, 3): 1, (3, 4): 1, (3, 5): 1})
    black = Counter({(0, 1): 2, (4, 5): 2})
    graph = full_graph(nodes, blue, black)
    direct = {tuple(sorted(edge)) for edge in nx.bridges(graph)}
    need(direct == {(2, 3)}, "non-leaf bridge hostile")
    components, _, lifted = black_quotient(nodes, blue, black)
    need(lifted == direct, "hostile contraction")
    involution_realization(nodes, nodes, blue, black)
    need(all(graph.degree(v) == 3 for v in nodes), "cubic hostile")
    print("Abstract hostile: 6 SC nodes, all full multidegrees 3, "
          "unique bridge (2,3), both endpoint degrees 3; V4-realizable.")

    # Even black multiplicities, all fibre parity/mass constraints at n=4,
    # but impossible actual tournament wiring: transitive node gets a self-loop.
    blue = Counter({(0, 0): 1, (0, 1): 1})
    black = Counter({(1, 2): 2})
    fibres, total = involution_realization({0, 1, 2}, {0, 1}, blue, black)
    need(total == 8 and fibres == Counter({(0, 0): 3, (1, 0): 3,
                                         (2, 0): 1, (2, 1): 1}),
         "abstract n4 mass/parity hostile")
    print("Abstract n=4 mass hostile: SC fibres 3,3; NS sheet fibres 1,1; "
          "4 blue and 4 black flags. V4-realizable; impossible tournament "
          "fibres because actual SC fibres are 1,5.")

    classes, lines = load(5)
    nodes, blue, black = merge(classes, lines)
    actual_fibres, _ = involution_realization(
        nodes, {v for v in nodes if classes[v]["SC"]}, blue, black)
    need(blue[0, 11] == blue[2, 6] == 1, "identified n5 switch edges")
    blue[0, 11] -= 1
    blue[2, 6] -= 1
    blue[2, 11] += 1
    blue[0, 6] += 1
    altered_fibres, total = involution_realization(
        nodes, {v for v in nodes if classes[v]["SC"]}, +blue, black)
    need(altered_fibres == actual_fibres and total == 64,
         "all actual ordinary fibre sizes retained")
    need(classes[11]["H"] == 1 and classes[0]["H"] == 9
         and classes[2]["H"] == 13, "fixed transitive partner obstruction")
    print("Stronger n=5 hostile: blue switch (0,11),(2,6) -> (2,11),(0,6) "
          "preserves every actual fibre and blue count, all black multiplicities, "
          "and V4 realization; impossible labelled tournament wiring because "
          "the transitive node now points to H=13 instead of its forced H=9 partner.")


def main():
    print("STATUS: exact n=3..7 committed-census bridge audit; abstract V4 construction.")
    print("Inherited full line data: THM-4731 TSV producer; independent small path-cube producer.")
    for n in range(3, 8):
        _, graph = census(n)
        if n <= 6:
            observed = independent_small_census(n)
            expected = (len(graph), graph.number_of_edges(), len(list(nx.bridges(graph))),
                        Counter(dict(graph.degree()).values()))
            need(observed == expected, "independent path-cube census")
    print("Declared interpretation D_n(z)=sum z^degree, M_n=degree(D_n):")
    for n in range(3, 9):
        degree_polynomial(n)
    hostile()
    print("PASS: direct lowlink, independent edge deletion, black contraction, "
          "actual fibre reconstruction, small independent path-cube controls.")


if __name__ == "__main__":
    main()
