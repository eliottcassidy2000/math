#!/usr/bin/env python3
"""S132: Farey mutations meet graph carriers and second-moment proof roles.

This extends HYP-2931.  The four mutated Farey payloads

    p + q, p*q, q^p, p^q

are not replacements for the binding denominator q in

    M(S)-1/14 = (14p-q)/(14q).

Instead they are side channels.  This script asks where those channels should
land in the existing LRC14 carrier stack:

* q and p+q: additive Farey / n+2 recursion lane;
* p*q: multiplicative coimage / product-ledger lane;
* q^p and p^q: magnitude-leak stress tests;
* octahedron L(K4): support-six repeated-packet current carrier;
* Clebsch / halved 5-cube: folded residual-mask covariance carrier;
* Paley-Zygmund / second moment: existence lower-bound tool, but too lossy
  unless upgraded to the HYP-2823 degree-4 moment-region target.

The script is deliberately finite and exact where possible.  It produces a
role tournament on proof carriers, not a proof of LRC14.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import comb, log


THR = F(1, 14)


def popcount(x: int) -> int:
    return x.bit_count()


Graph = dict[object, set[object]]


def add_edge(graph: Graph, a: object, b: object) -> None:
    graph.setdefault(a, set()).add(b)
    graph.setdefault(b, set()).add(a)


def folded_five_cube_clebsch() -> Graph:
    """Clebsch graph as the folded 5-cube."""
    verts = sorted({min(x, x ^ 31) for x in range(32)})
    graph: Graph = {v: set() for v in verts}
    for a, b in combinations(verts, 2):
        d = popcount(a ^ b)
        if min(d, 5 - d) == 1:
            add_edge(graph, a, b)
    return graph


def halved_cube(dim: int = 5) -> Graph:
    """Halved cube on even-parity words, adjacent at Hamming distance 2."""
    verts = [x for x in range(1 << dim) if popcount(x) % 2 == 0]
    graph: Graph = {v: set() for v in verts}
    for a, b in combinations(verts, 2):
        if popcount(a ^ b) == 2:
            add_edge(graph, a, b)
    return graph


def octahedron_line_k4() -> Graph:
    """Octahedron as the line graph L(K4)."""
    k4_edges = tuple(combinations(range(4), 2))
    graph: Graph = {edge: set() for edge in k4_edges}
    for a, b in combinations(k4_edges, 2):
        if set(a) & set(b):
            add_edge(graph, a, b)
    return graph


def num_edges(graph: Graph) -> int:
    return sum(len(neigh) for neigh in graph.values()) // 2


def connected_components(graph: Graph) -> list[set[object]]:
    unseen = set(graph)
    comps: list[set[object]] = []
    while unseen:
        start = unseen.pop()
        comp = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in graph[cur]:
                if nxt in unseen:
                    unseen.remove(nxt)
                    comp.add(nxt)
                    stack.append(nxt)
        comps.append(comp)
    return comps


def triangle_count(graph: Graph) -> int:
    count = 0
    nodes = list(graph)
    for a, b, c in combinations(nodes, 3):
        if b in graph[a] and c in graph[a] and c in graph[b]:
            count += 1
    return count


def complement_graph(graph: Graph) -> Graph:
    nodes = list(graph)
    comp: Graph = {v: set() for v in nodes}
    for a, b in combinations(nodes, 2):
        if b not in graph[a]:
            add_edge(comp, a, b)
    return comp


def graph_stats(graph: Graph) -> dict[str, object]:
    comps = len(connected_components(graph))
    degree_hist = Counter(len(neigh) for neigh in graph.values())
    triangles = triangle_count(graph)
    cycle_rank = num_edges(graph) - len(graph) + comps
    return {
        "vertices": len(graph),
        "edges": num_edges(graph),
        "degree_hist": dict(sorted(degree_hist.items())),
        "triangles": triangles,
        "cycle_rank": cycle_rank,
    }


def common_neighbor_hist(graph: Graph) -> tuple[Counter[int], Counter[int]]:
    adj = Counter()
    non = Counter()
    for a, b in combinations(graph, 2):
        c = len(graph[a] & graph[b])
        if b in graph[a]:
            adj[c] += 1
        else:
            non[c] += 1
    return adj, non


def closed_neighborhood_design(graph: Graph) -> dict[str, object]:
    blocks = {v: {v, *graph[v]} for v in graph}
    point_rep = Counter()
    pair_rep = Counter()
    for block in blocks.values():
        for v in block:
            point_rep[v] += 1
        for a, b in combinations(sorted(block), 2):
            pair_rep[(a, b)] += 1
    return {
        "blocks": len(blocks),
        "block_size_hist": dict(Counter(len(b) for b in blocks.values())),
        "point_rep_hist": dict(Counter(point_rep.values())),
        "pair_rep_hist": dict(Counter(pair_rep.values())),
    }


def halved_cube_is_complement_clebsch() -> bool:
    """Direct even-word model for complement(Clebsch).

    In dimension 5, each folded-cube complement pair has exactly one even word.
    Two even words have Hamming distance 0, 2, or 4.  Folded-cube adjacency is
    min(d,5-d)=1, which on even words is exactly d=4.  Therefore the complement
    has edges exactly at d=2, the halved-cube rule.
    """
    half = halved_cube(5)
    cleb = folded_five_cube_clebsch()
    comp = complement_graph(cleb)
    # Re-key complement(Clebsch) by the even representative of each folded pair.
    rekeyed: Graph = {}
    for v, neigh in comp.items():
        even_v = v if popcount(v) % 2 == 0 else v ^ 31
        rekeyed.setdefault(even_v, set())
        for w in neigh:
            even_w = w if popcount(w) % 2 == 0 else w ^ 31
            rekeyed[even_v].add(even_w)
    return all(half[v] == rekeyed[v] for v in half)


def unit_excess_chain(limit: int = 10) -> list[dict[str, object]]:
    rows = []
    for p in range(1, limit + 1):
        q = 14 * p - 1
        rows.append(
            {
                "p": p,
                "value": F(p, q),
                "q": q,
                "sum": p + q,
                "product": p * q,
                "log_qpow": p * log(q),
                "log_ppow": q * log(p) if p > 1 else 0.0,
            }
        )
    return rows


def finite_differences(values: list[int]) -> tuple[list[int], list[int]]:
    first = [b - a for a, b in zip(values, values[1:])]
    second = [b - a for a, b in zip(first, first[1:])]
    return first, second


def binomial_miss_model(k: int) -> dict[str, F]:
    """Decorrelated six-sector empty model used only to calibrate PZ loss.

    Each inner sector is independently empty with probability a=(6/7)^k.  This
    is not the exact LRC row law.  It gives a transparent Paley-Zygmund scale:
    for N=sum 1[sector empty], P(N>0) is the event lower-bounded by second
    moments.
    """
    a = F(6, 7) ** k
    en = 6 * a
    en2 = en + 2 * comb(6, 2) * a * a
    pz = en * en / en2 if en2 else F(0)
    exact = 1 - (1 - a) ** 6
    return {"a": a, "E": en, "E2": en2, "PZ": pz, "exact_union": exact}


@dataclass(frozen=True)
class Carrier:
    name: str
    role: str
    theorem_scale: int
    label_retention: int
    lrc_relevance: int
    scalarization_risk: int

    @property
    def score(self) -> tuple[int, int, int, int]:
        # Higher theorem/retention/relevance is better; lower scalarization
        # risk is better, so use its negative in the role tournament.
        return (
            self.theorem_scale,
            self.label_retention,
            self.lrc_relevance,
            -self.scalarization_risk,
        )


def carrier_list() -> list[Carrier]:
    return [
        Carrier("q_binding_scale", "theorem coordinate", 4, 4, 4, 0),
        Carrier("p_plus_q_additive", "n+2 / Stern-Brocot ledger", 3, 3, 4, 1),
        Carrier("p_times_q_product", "multiplicative coimage side channel", 2, 4, 3, 2),
        Carrier("octahedron_LK4", "support-six current / face-curl carrier", 2, 4, 4, 1),
        Carrier("Clebsch_halfcube", "folded residual covariance carrier", 1, 4, 3, 2),
        Carrier("PZ_second_moment", "existence lower-bound / moment gateway", 1, 2, 3, 3),
        Carrier("power_payloads", "magnitude-leak stress tests", 0, 2, 2, 4),
    ]


def role_tournament(carriers: list[Carrier]) -> int:
    """Orient carrier i->j if i has lexicographically better role score."""
    mask = 0
    bit = 0
    for i in range(len(carriers)):
        for j in range(i + 1, len(carriers)):
            if carriers[i].score >= carriers[j].score:
                mask |= 1 << bit
            bit += 1
    return mask


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i < j:
        bit = i * n - i * (i + 1) // 2 + (j - i - 1)
        return bool(mask & (1 << bit))
    return not edge(mask, n, j, i)


def tournament_fingerprint(mask: int, n: int) -> dict[str, object]:
    scores = [sum(1 for j in range(n) if i != j and edge(mask, n, i, j)) for i in range(n)]
    c3 = 0
    for tri in combinations(range(n), 3):
        local = [sum(1 for j in tri if i != j and edge(mask, n, i, j)) for i in tri]
        if sorted(local) == [1, 1, 1]:
            c3 += 1
    return {"score_hist": dict(Counter(scores)), "c3": c3}


def hamiltonian_order(carriers: list[Carrier]) -> list[str]:
    return [c.name for c in sorted(carriers, key=lambda c: c.score, reverse=True)]


def print_chain() -> None:
    print("[1] Unit-excess Farey child chain")
    rows = unit_excess_chain(10)
    print("  p/q = p/(14p-1); q moves by +14, matching the apex-14 n+2 lane.")
    print("  p   value   q   p+q   p*q   log(q^p)  log(p^q)")
    for r in rows:
        print(
            f"  {r['p']:2d} {str(r['value']):>7s} {r['q']:3d} "
            f"{r['sum']:5d} {r['product']:5d} "
            f"{r['log_qpow']:9.3f} {r['log_ppow']:9.3f}"
        )
    q_vals = [int(r["q"]) for r in rows]
    s_vals = [int(r["sum"]) for r in rows]
    prod_vals = [int(r["product"]) for r in rows]
    print(f"  Delta q: {finite_differences(q_vals)[0]}")
    print(f"  Delta(p+q): {finite_differences(s_vals)[0]}")
    print(f"  Delta product: {finite_differences(prod_vals)[0]}")
    print(f"  Delta^2 product: {finite_differences(prod_vals)[1]}")
    print("  readout: p+q is additive/linear, p*q is a quadratic area ledger,")
    print("           power payloads are scale amplifiers rather than local addresses.")
    print()


def print_graphs() -> None:
    print("[2] Graph carriers requested by the prompt")
    octa = octahedron_line_k4()
    cleb = folded_five_cube_clebsch()
    half = halved_cube(5)
    adj, non = common_neighbor_hist(cleb)
    print(f"  octahedron L(K4): {graph_stats(octa)}")
    print("    readout: six packet vertices are K4 edges; cycle rank 7 is the apex-7 curl module.")
    print(f"  Clebsch folded 5-cube: {graph_stats(cleb)}")
    print(f"    common-neighbor hist adjacent={dict(adj)} nonadjacent={dict(non)}")
    print(f"    closed-neighborhood design={closed_neighborhood_design(cleb)}")
    print(f"  halved 5-cube: {graph_stats(half)}")
    print(
        "    identical to complement(Clebsch) under even representatives: "
        f"{halved_cube_is_complement_clebsch()}"
    )
    print("  readout: Clebsch/half-cube is the cut/covariance side; octahedron is")
    print("           the support-six curl side.  Neither should be collapsed to a")
    print("           raw runner tournament.")
    print()


def print_pz() -> None:
    print("[3] Paley-Zygmund / second-moment calibration")
    print("  Toy model: six independent empty-sector indicators with a=(6/7)^k.")
    print("  PZ lower-bounds P(N>0) by E[N]^2/E[N^2].  It is an existence tool,")
    print("  but the repo's live wide cap needs HYP-2823's degree-4 moment target.")
    print("  k   a          E[N]       exact P(N>0)  PZ bound   loss")
    for k in range(8, 13):
        row = binomial_miss_model(k)
        loss = row["exact_union"] - row["PZ"]
        print(
            f"  {k:2d} {float(row['a']):10.6f} {float(row['E']):10.6f} "
            f"{float(row['exact_union']):13.6f} {float(row['PZ']):10.6f} "
            f"{float(loss):8.6f}"
        )
    print("  readout: second moment can guarantee nonzero witness mass, but it loses")
    print("           too much for tight LRC14 cap comparison unless upgraded to the")
    print("           labelled factorial-moment/Delsarte lane.")
    print()


def print_tournament() -> None:
    print("[4] Tournament Analysis on seven proof carriers")
    carriers = carrier_list()
    mask = role_tournament(carriers)
    fp = tournament_fingerprint(mask, len(carriers))
    print("  vertices: proof carriers, not runners")
    print("  observable: (theorem scale, label retention, LRC relevance, scalarization risk)")
    print("  switch/gauge: lexicographically better role score wins")
    print("  tie Hamiltonian path: listed carrier order")
    for c in carriers:
        print(f"    {c.name:22s} score={c.score} role={c.role}")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']}")
    print("  Hamiltonian role order:")
    print("    " + " > ".join(hamiltonian_order(carriers)))
    print()


def main() -> None:
    print("S132 LRC14 FAREY GRAPH / PZ CARRIER SYNTHESIS")
    print("=" * 78)
    print("[Assumption challenge]")
    print("  considered vertices: runners, Farey fractions, fraction payloads,")
    print("    residual masks, K4 edges, half-cube cuts, sector-empty events, and")
    print("    proof obligations.")
    print("  chosen quotient: a role tournament on proof carriers, backed by exact")
    print("    graph and sequence fingerprints.")
    print("  preserved predicate: which carrier keeps the data needed for the")
    print("    LRC14 binding-scale / state-lift proof target.")
    print("  destroyed predicate: raw row-specific witnesses; those remain in")
    print("    HYP-2930/HYP-2931 and the exact-period packet ledgers.")
    print()
    print_chain()
    print_graphs()
    print_pz()
    print_tournament()
    print("[Proof readout]")
    print("  HYP-2931's denominator verdict survives: q is the theorem coordinate.")
    print("  The n+2 analogy is literal on unit-excess Farey children: q -> q+14.")
    print("  The n*2/product analogy should be a side-channel area/coimage ledger,")
    print("  not an ordering replacement for q.")
    print("  Octahedral and Clebsch/half-cube carriers explain where labelled packet")
    print("  covariance and curl data should go after the fraction address is kept.")
    print("  Paley-Zygmund is useful only as an existence gateway; for the cap route")
    print("  it must be upgraded to HYP-2823's degree-4 factorial moment inequality.")
    print("  Next theorem target: a reduced |M14|<=6 atom should carry the labelled")
    print("  address (e,q,p+q,p*q) plus Clebsch/octahedral packet data into either")
    print("  HYP-2930's non-tight Farey class or HYP-2908's forbidden H=7 state lift.")


if __name__ == "__main__":
    main()
