#!/usr/bin/env python3
"""Median-sidecar audit for the LRC14 proof-route spine.

This is a proof-interface scout, not a proof of LRC14.  It turns the S233
Desargues-median warning and the forum owner/root sidecar post into an
executable checklist:

* raw even-cycle/incidence shadows can have empty median centers;
* adding the correct owner/root/sidecar coordinate should replace an empty or
  ambiguous center by a unique proof-state center;
* the sidecar-ranking tournament uses proof obligations as vertices, not
  runners or raw graph nodes.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations


def cycle_graph(n: int) -> dict[str, set[str]]:
    return {
        str(i): {str((i - 1) % n), str((i + 1) % n)}
        for i in range(n)
    }


def hypercube(dim: int) -> dict[str, set[str]]:
    vertices = [format(i, f"0{dim}b") for i in range(1 << dim)]
    graph = {v: set() for v in vertices}
    for v in vertices:
        for bit in range(dim):
            w = list(v)
            w[bit] = "1" if w[bit] == "0" else "0"
            graph[v].add("".join(w))
    return graph


def all_pairs_dist(graph: dict[str, set[str]]) -> dict[str, dict[str, int]]:
    dists: dict[str, dict[str, int]] = {}
    for source in graph:
        dist = {source: 0}
        q = deque([source])
        while q:
            cur = q.popleft()
            for nxt in graph[cur]:
                if nxt not in dist:
                    dist[nxt] = dist[cur] + 1
                    q.append(nxt)
        dists[source] = dist
    return dists


def interval(dists: dict[str, dict[str, int]], a: str, b: str) -> set[str]:
    dab = dists[a][b]
    return {x for x in dists if dists[a][x] + dists[x][b] == dab}


def median_centers(graph: dict[str, set[str]], a: str, b: str, c: str) -> set[str]:
    dists = all_pairs_dist(graph)
    return interval(dists, a, b) & interval(dists, b, c) & interval(dists, c, a)


@dataclass(frozen=True)
class RouteCase:
    case_id: str
    coarse_fiber: str
    route_triple: tuple[str, str, str]
    coarse_shadow: str
    before_status: str
    after_status: str
    first_missing_sidecar: str
    added_sidecars: tuple[str, ...]
    proof_readout: str


CASES = [
    RouteCase(
        "q23_endpoint_owner_center",
        "B18Z6 q=23 petal/covering diagonal",
        ("exact_M_diagonal", "haar_zeta_square", "endpoint_route"),
        "same coarse endpoint word and M=2/23",
        "empty",
        "unique",
        "endpoint_owner_strip",
        ("exact_M_zeta", "endpoint_owner_strip", "safe_component_body"),
        "Exact M and zeta diagnose the diagonal, but the owner strip names the center.",
    ),
    RouteCase(
        "a000568_rootless_center",
        "first rooted-perspective failure U(6)-P(5)=8",
        ("node_root", "directed_edge_sector", "cycle_conflict"),
        "node perspective cache sees only P(5)=48",
        "empty",
        "unique",
        "rootless_cycle_object",
        ("rootless_cycle_object", "ordered_pair_sector_deck", "observer_cut_orbit"),
        "The [3,3] Burnside term has no fixed vertex, so the root object must change.",
    ),
    RouteCase(
        "desargues_beal_owner_center",
        "post-rectangle Desargues/Beal finalizer",
        ("incidence_route", "owner_route", "arithmetic_collision_route"),
        "girth-six incidence address without common-owner gate",
        "empty",
        "unique",
        "desargues_girth6_residue",
        ("desargues_girth6_residue", "beal_common_owner_gate", "endpoint_owner_strip"),
        "A girth-six address is useful only when paired with owner/factor payload.",
    ),
    RouteCase(
        "fejer_haar_ramanujan_center",
        "harmonic certificate handoff fiber",
        ("fejer_interval", "haar_zeta", "ramanujan_period"),
        "small negative/scalar certificate without value-origin type",
        "multiple",
        "unique",
        "value_origin_type",
        ("value_origin_type", "primitive_period_deck", "exact_M_zeta"),
        "Multiple centers mean the same number is playing incompatible proof roles.",
    ),
    RouteCase(
        "observer_deletion_rectangle_center",
        "observer/deletion/rectangle transport fiber",
        ("observer_cut", "deletion_parent", "rectangle_hourglass"),
        "observer quotient before the cut payload orbit is named",
        "empty",
        "unique",
        "observer_cut_orbit",
        ("observer_cut_orbit", "deletion_parent_fiber", "rectangle_hourglass_residue"),
        "The center is the orbit of the cut payload, not the raw rooted class.",
    ),
    RouteCase(
        "pair_good_barcode_center",
        "pair-good decoy / barcode / normal-fan fiber",
        ("pair_good_generator", "barcode_owner", "normal_fan_support"),
        "raw pair-good count and barcode shadow",
        "empty",
        "unique",
        "active_owner_barcode_support",
        ("active_owner_barcode_support", "endpoint_owner_strip", "residual_capacitor_cut"),
        "The proof object is the active owner/barcode tooth, not the decoy count.",
    ),
]


SIDECAR_TIE_PATH = [
    "endpoint_owner_strip",
    "observer_cut_orbit",
    "rootless_cycle_object",
    "active_owner_barcode_support",
    "desargues_girth6_residue",
    "beal_common_owner_gate",
    "value_origin_type",
    "residual_capacitor_cut",
    "rectangle_hourglass_residue",
    "primitive_period_deck",
    "exact_M_zeta",
    "raw_scalar_count",
]


SIDECAR_REPAIRS = {
    sidecar: {
        case.case_id
        for case in CASES
        if sidecar in case.added_sidecars or sidecar == case.first_missing_sidecar
    }
    for sidecar in SIDECAR_TIE_PATH
}

SIDECAR_PRESERVES = {
    "endpoint_owner_strip": {"owner", "route", "boundary_status"},
    "observer_cut_orbit": {"observer", "route", "boundary_status"},
    "rootless_cycle_object": {"root_object", "extension_legality"},
    "active_owner_barcode_support": {"owner", "topology", "route"},
    "desargues_girth6_residue": {"incidence", "topology"},
    "beal_common_owner_gate": {"owner", "arithmetic_collision"},
    "value_origin_type": {"certificate_role", "small_value_origin"},
    "residual_capacitor_cut": {"route", "first_cut"},
    "rectangle_hourglass_residue": {"cycle_payload", "transport"},
    "primitive_period_deck": {"period", "route"},
    "exact_M_zeta": {"exact_scale", "haar_cocycle"},
    "raw_scalar_count": set(),
}

SIDECAR_DEBT = {
    "endpoint_owner_strip": 1,
    "observer_cut_orbit": 1,
    "rootless_cycle_object": 2,
    "active_owner_barcode_support": 2,
    "desargues_girth6_residue": 2,
    "beal_common_owner_gate": 2,
    "value_origin_type": 1,
    "residual_capacitor_cut": 2,
    "rectangle_hourglass_residue": 2,
    "primitive_period_deck": 1,
    "exact_M_zeta": 1,
    "raw_scalar_count": 5,
}


def sidecar_rank(sidecar: str) -> tuple[int, int, int, int]:
    return (
        len(SIDECAR_REPAIRS[sidecar]),
        len(SIDECAR_PRESERVES[sidecar]),
        -SIDECAR_DEBT[sidecar],
        -SIDECAR_TIE_PATH.index(sidecar),
    )


def tournament_edges(vertices: list[str]) -> dict[str, set[str]]:
    edges = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner, loser = (a, b) if sidecar_rank(a) >= sidecar_rank(b) else (b, a)
        edges[winner].add(loser)
    return edges


def directed_3cycles(edges: dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(edges, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            count += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            count += 1
    return count


def scc_sizes(edges: dict[str, set[str]]) -> list[int]:
    vertices = list(edges)
    rev = {v: set() for v in vertices}
    for u, outs in edges.items():
        for v in outs:
            rev[v].add(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs1(v: str) -> None:
        seen.add(v)
        for w in edges[v]:
            if w not in seen:
                dfs1(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs1(v)

    assigned: set[str] = set()
    sizes = []

    for v in reversed(order):
        if v in assigned:
            continue
        comp: set[str] = set()

        def dfs2_assigned(u: str) -> None:
            assigned.add(u)
            comp.add(u)
            for w in rev[u]:
                if w not in assigned:
                    dfs2_assigned(w)

        dfs2_assigned(v)
        sizes.append(len(comp))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: dict[str, set[str]]) -> int:
    vertices = list(edges)
    n = len(vertices)
    index = {v: i for i, v in enumerate(vertices)}
    dp = defaultdict(int)
    for v in vertices:
        dp[(1 << index[v], index[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if ways == 0:
                continue
            last_v = vertices[last]
            for nxt_v in edges[last_v]:
                nxt = index[nxt_v]
                if mask & (1 << nxt) == 0:
                    dp[(mask | (1 << nxt), nxt)] += ways
    full = (1 << n) - 1
    return sum(dp[(full, last)] for last in range(n))


def main() -> None:
    print("LRC14 median owner/root sidecar audit")
    print("Scope: proof-interface checklist; not a proof of LRC14.")
    print()

    c6_centers = median_centers(cycle_graph(6), "0", "2", "4")
    q3_centers = median_centers(hypercube(3), "000", "110", "101")
    print("Toy median lift")
    print(f"  C6 route triple (0,2,4) centers: {sorted(c6_centers)}")
    print(f"  Q3 sidecar triple (000,110,101) centers: {sorted(q3_centers)}")
    print()

    before_hist = Counter(case.before_status for case in CASES)
    after_hist = Counter(case.after_status for case in CASES)
    print("Route-triple audit table")
    for case in CASES:
        print(f"- {case.case_id}")
        print(f"  coarse_fiber: {case.coarse_fiber}")
        print(f"  route_triple: {case.route_triple}")
        print(f"  coarse_shadow: {case.coarse_shadow}")
        print(f"  median_center_status: {case.before_status} -> {case.after_status}")
        print(f"  first_missing_sidecar: {case.first_missing_sidecar}")
        print(f"  sidecars_added: {case.added_sidecars}")
        print(f"  readout: {case.proof_readout}")
    print()
    print(f"before_status_hist: {dict(sorted(before_hist.items()))}")
    print(f"after_status_hist: {dict(sorted(after_hist.items()))}")
    print()

    vertices = SIDECAR_TIE_PATH[:]
    edges = tournament_edges(vertices)
    score_hist = Counter(len(edges[v]) for v in vertices)
    print("Sidecar tournament")
    print("  vertices are sidecar/proof-obligation objects, not runners.")
    print("  pairwise observable: repaired route triples, preserved predicates, debt cost.")
    print(f"  score_hist: {dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles: {directed_3cycles(edges)}")
    print(f"  scc_sizes: {scc_sizes(edges)}")
    print(f"  hamiltonian_path_count: {hamiltonian_path_count(edges)}")
    print("  ranked_path:")
    for sidecar in sorted(vertices, key=sidecar_rank, reverse=True):
        print(
            "    "
            f"{sidecar}: repairs={len(SIDECAR_REPAIRS[sidecar])}, "
            f"preserves={len(SIDECAR_PRESERVES[sidecar])}, "
            f"debt={SIDECAR_DEBT[sidecar]}"
        )
    print()

    print("Assumption challenge")
    print("  considered_vertices: runners, gaps, circle sections, section boundaries,")
    print("    wall crossings, residues, cover arcs, Fourier modes, matroid circuits,")
    print("    proof obligations, owner objects, root objects")
    print("  chosen_vertices: proof obligations plus owner/root/sidecar objects")
    print("  preserved_lrc_predicate: boundary/open status and route schedulability")
    print("  destroyed_information: median center when owner, root, value-origin,")
    print("    observer-cut, rectangle/hourglass, or cocycle payload is forgotten")
    print()

    print("Readout")
    print("- The useful next table is not a larger scalar count.  It is the")
    print("  case_id -> first_missing_sidecar ledger above, run on actual HYP-2963")
    print("  coarse fibers.")
    print("- q=23/B18Z6 points first at endpoint_owner_strip.")
    print("- A000568/P(5)->U(6) points first at rootless_cycle_object.")
    print("- Desargues/Beal points first at girth-six incidence plus common owner.")
    print("- Multiple centers should be treated as value-origin ambiguity before")
    print("  they are treated as new residual debt.")


if __name__ == "__main__":
    main()
