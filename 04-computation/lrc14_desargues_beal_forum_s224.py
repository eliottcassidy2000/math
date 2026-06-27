"""
lrc14_desargues_beal_forum_s224.py

Small exact scout for S224:

1. Build the Desargues graph as the bipartite double cover of the Petersen
   graph: two copies of the 10 two-subsets of {0,1,2,3,4}, with incidence by
   disjointness.  This is the Levi graph of the Desargues configuration.
2. Check graph invariants that matter for LRC proof-carrier use: cubic,
   bipartite, no rectangles, girth 6, diameter, automorphism orbits, and
   cycle counts.  The point is not to use a graph scalar, but to expose the
   first non-rectangle incidence residue after S217 rectangle/hourglass data.
3. Run a tiny Beal-style perfect-power scout.  The computation is not evidence
   for Beal's conjecture; it is a sanity check for the common-owner metaphor:
   every exact triple-power equality found in the box has a shared base factor.
4. Emit a proof-carrier tournament whose vertices are sidecar families, not
   runners.  This obeys the repo's Tournament Analysis default.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from typing import Dict, Iterable, List, Tuple

import networkx as nx


Vertex = Tuple[str, Tuple[int, int]]


def two_subsets(n: int = 5) -> List[Tuple[int, int]]:
    return list(combinations(range(n), 2))


def desargues_graph() -> nx.Graph:
    graph = nx.Graph()
    pairs = two_subsets()
    for side in ("P", "L"):
        for pair in pairs:
            graph.add_node((side, pair), side=side, pair=pair)
    for p in pairs:
        for line in pairs:
            if set(p).isdisjoint(line):
                graph.add_edge(("P", p), ("L", line))
    return graph


def generalized_petersen_10_3() -> nx.Graph:
    graph = nx.Graph()
    for i in range(10):
        graph.add_node(("outer", i))
        graph.add_node(("inner", i))
    for i in range(10):
        graph.add_edge(("outer", i), ("outer", (i + 1) % 10))
        graph.add_edge(("outer", i), ("inner", i))
        graph.add_edge(("inner", i), ("inner", (i + 3) % 10))
    return graph


def girth(graph: nx.Graph) -> int:
    best = None
    for start in graph.nodes:
        dist = {start: 0}
        parent = {start: None}
        q = deque([start])
        while q:
            u = q.popleft()
            for v in graph.neighbors(u):
                if v not in dist:
                    parent[v] = u
                    dist[v] = dist[u] + 1
                    q.append(v)
                elif parent[u] != v and parent[v] != u:
                    length = dist[u] + dist[v] + 1
                    best = length if best is None else min(best, length)
    if best is None:
        raise ValueError("acyclic graph has no girth")
    return best


def canonical_cycle(cycle: Tuple[Vertex, ...]) -> Tuple[Vertex, ...]:
    seq = list(cycle)
    rotations = []
    for orient in (seq, list(reversed(seq))):
        for i in range(len(seq)):
            rotations.append(tuple(orient[i:] + orient[:i]))
    return min(rotations, key=repr)


def cycle_counts(graph: nx.Graph, max_len: int = 10) -> Dict[int, int]:
    cycles = set()
    nodes = list(graph.nodes)

    def dfs(start: Vertex, current: Vertex, path: List[Vertex]) -> None:
        if len(path) > max_len:
            return
        for nxt in graph.neighbors(current):
            if nxt == start and len(path) >= 3:
                cycles.add(canonical_cycle(tuple(path)))
            elif nxt not in path and repr(nxt) >= repr(start):
                dfs(start, nxt, path + [nxt])

    for node in nodes:
        dfs(node, node, [node])
    return dict(sorted(Counter(len(c) for c in cycles).items()))


def automorphism_summary(graph: nx.Graph) -> Dict[str, object]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph, graph)
    autos = list(matcher.isomorphisms_iter())

    def orbit_of(items: Iterable[object], image_fn) -> List[int]:
        unseen = set(items)
        sizes = []
        while unseen:
            seed = next(iter(unseen))
            orb = {image_fn(auto, seed) for auto in autos}
            sizes.append(len(orb))
            unseen -= orb
        return sorted(sizes)

    nodes = list(graph.nodes)
    edges = [tuple(sorted(edge, key=repr)) for edge in graph.edges]

    def edge_image(auto, edge):
        return tuple(sorted((auto[edge[0]], auto[edge[1]]), key=repr))

    return {
        "automorphism_count": len(autos),
        "vertex_orbit_sizes": orbit_of(nodes, lambda auto, node: auto[node]),
        "edge_orbit_sizes": orbit_of(edges, edge_image),
    }


def beal_scout(base_bound: int = 80, exponent_bound: int = 7) -> Dict[str, object]:
    powers: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for base in range(2, base_bound + 1):
        for exp in range(3, exponent_bound + 1):
            powers[base**exp].append((base, exp))

    values = sorted(powers)
    hits = []
    primitive_hits = []
    for i, va in enumerate(values):
        for vb in values[i:]:
            s = va + vb
            if s not in powers:
                continue
            for a, x in powers[va]:
                for b, y in powers[vb]:
                    for c, z in powers[s]:
                        common = gcd(gcd(a, b), c)
                        row = {
                            "A": a,
                            "x": x,
                            "B": b,
                            "y": y,
                            "C": c,
                            "z": z,
                            "gcd": common,
                        }
                        hits.append(row)
                        if common == 1:
                            primitive_hits.append(row)
    hits.sort(key=lambda r: (r["C"] ** r["z"], r["A"], r["B"], r["C"], r["x"], r["y"], r["z"]))
    return {
        "base_bound": base_bound,
        "exponent_bound": exponent_bound,
        "hit_count": len(hits),
        "primitive_hit_count": len(primitive_hits),
        "first_hits": hits[:12],
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    features: Tuple[int, ...]


def tournament_fingerprint(carriers: List[Carrier]) -> Dict[str, object]:
    # Larger feature vector is better; tie path is the listed order.
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    for i, ci in enumerate(carriers):
        for j, cj in enumerate(carriers):
            if i == j:
                continue
            if ci.features > cj.features:
                adj[i][j] = True
            elif ci.features == cj.features:
                adj[i][j] = i < j
            else:
                adj[j][i] = True
    score_hist = Counter(sum(row) for row in adj)
    cycles3 = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cycles3 += 1

    digraph = nx.DiGraph()
    digraph.add_nodes_from(range(n))
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                digraph.add_edge(i, j)
    scc_sizes = sorted((len(c) for c in nx.strongly_connected_components(digraph)), reverse=True)

    # Count Hamiltonian paths by DP over subsets; n is tiny.
    dp = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    h_paths = sum(dp.get((full, last), 0) for last in range(n))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles3,
        "scc_sizes": scc_sizes,
        "hamiltonian_path_count": h_paths,
        "tie_path": [c.name for c in carriers],
    }


def main() -> None:
    graph = desargues_graph()
    gp = generalized_petersen_10_3()
    auto = automorphism_summary(graph)
    cycles = cycle_counts(graph, 10)
    beal = beal_scout()
    carriers = [
        Carrier("labelled_packet_sheaf", (9, 9, 9, 9, 9, 9, 8, 9)),
        Carrier("observer_cut_orbit_ledger", (8, 8, 8, 8, 8, 7, 6, 8)),
        Carrier("desargues_girth6_incidence_residue", (7, 7, 7, 9, 9, 8, 5, 8)),
        Carrier("beal_common_owner_gate", (7, 7, 8, 8, 6, 9, 5, 7)),
        Carrier("endpoint_owner_strip", (8, 8, 6, 9, 7, 6, 7, 8)),
        Carrier("residual_capacitor_min_cut", (7, 8, 7, 7, 8, 6, 6, 7)),
        Carrier("haar_zeta_cocycle", (7, 7, 6, 6, 8, 5, 9, 7)),
        Carrier("fejer_interval_certificate", (8, 7, 6, 5, 5, 4, 9, 6)),
        Carrier("raw_desargues_scalar", (2, 2, 2, 2, 3, 2, 2, 3)),
        Carrier("raw_beal_scalar", (2, 2, 4, 3, 2, 4, 2, 3)),
    ]
    fp = tournament_fingerprint(carriers)

    print("S224 Desargues/Beal LRC14 forum scout")
    print("=" * 72)
    print("Desargues graph model")
    print(f"nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")
    print(f"degree_hist={dict(sorted(Counter(dict(graph.degree()).values()).items()))}")
    print(f"bipartite={nx.is_bipartite(graph)}")
    print(f"girth={girth(graph)} diameter={nx.diameter(graph)}")
    print(f"isomorphic_to_generalized_petersen_G_10_3={nx.is_isomorphic(graph, gp)}")
    print(f"cycle_counts_len_<=10={cycles}")
    print(f"automorphism_count={auto['automorphism_count']}")
    print(f"vertex_orbit_sizes={auto['vertex_orbit_sizes']}")
    print(f"edge_orbit_sizes={auto['edge_orbit_sizes']}")
    print()

    print("Beal-style perfect-power scout")
    print(f"base_bound={beal['base_bound']} exponent_bound={beal['exponent_bound']}")
    print(f"hit_count={beal['hit_count']} primitive_hit_count={beal['primitive_hit_count']}")
    print("first_hits:")
    for row in beal["first_hits"]:
        print(
            "  "
            f"{row['A']}^{row['x']} + {row['B']}^{row['y']} = "
            f"{row['C']}^{row['z']}  gcd={row['gcd']}"
        )
    print()

    print("Tournament Analysis")
    print("vertices=proof carriers, not runners")
    print(
        "pairwise_observable="
        "retained boundary/open status, route, exact scale, owner incidence, "
        "topology, arithmetic common-owner gate, harmonic certificate, "
        "visible automorphism orbit"
    )
    print("switch=lexicographic retained-payload vector; tie path listed below")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("tie_path=" + " > ".join(fp["tie_path"]))
    print()

    print("Assumption challenge")
    print(
        "Candidate vertices considered: runners, graph vertices, Desargues "
        "point-line flags, exponent triples, packet sidecar fields, and proof "
        "obligations.  The useful quotient uses sidecar/proof obligations: it "
        "preserves boundary/open status and route schedulability, while raw "
        "Desargues or Beal scalars destroy owner names, exact speeds, and the "
        "observer-cut payload."
    )
    print()
    print("Proof-use summary")
    print(
        "S217 rectangle/hourglass residues are 4-cycle/coboundary checks.  "
        "The Desargues graph has no 4-cycles and first exposes girth-6 "
        "incidence residue, so it is a candidate next sidecar for a residual "
        "that is invisible to rectangle tests.  Beal's common-factor pattern "
        "suggests the matching arithmetic gate: a primitive three-channel "
        "power collision should force a common owner/prime/packet coordinate "
        "or become named F7/THM-572 debt."
    )


if __name__ == "__main__":
    main()
