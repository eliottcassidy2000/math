#!/usr/bin/env python3
"""
unit_distance_n21_graph6_endpoint_audit_codex_s129.py

S129 addendum to HYP-2620.

Audit the five exact n=21, 57-edge graph6 cores from S614 at the level of
Hamiltonian endpoint masks.  The goal is to decide whether the remaining
unit-distance n=22 proof can be attacked by traceability/unit-spine arguments,
or whether those arguments are already too strong and the proof must move to
unit-cocyclic geometry.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import unit_distance_n22_tournament_lrc_s614 as s614  # noqa: E402


U_EXACT_20 = 54
U_EXACT_21 = 57


def bit_index(bit: int) -> int:
    return bit.bit_length() - 1


def graph_bits(adj: list[list[int]]) -> list[int]:
    bit_adj: list[int] = []
    for i, row in enumerate(adj):
        mask = 0
        for j, value in enumerate(row):
            if i != j and value:
                mask |= 1 << j
        bit_adj.append(mask)
    return bit_adj


def edge_count(adj: list[int]) -> int:
    return sum(bits.bit_count() for bits in adj) // 2


def endpoint_dp(adj: list[int]) -> list[int]:
    """For every mask, bitset of endpoints of Hamiltonian paths on mask."""
    n = len(adj)
    ends = [0] * (1 << n)
    for v in range(n):
        ends[1 << v] = 1 << v
    for mask in range(1 << n):
        e = ends[mask]
        while e:
            bit = e & -e
            e ^= bit
            v = bit_index(bit)
            nxts = adj[v] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                ends[mask | nb] |= nb
    return ends


def connected_components_after_delete(adj: list[int], drop: int) -> int:
    n = len(adj)
    remaining = ((1 << n) - 1) ^ (1 << drop)
    components = 0
    while remaining:
        components += 1
        bit = remaining & -remaining
        remaining ^= bit
        stack = [bit_index(bit)]
        while stack:
            v = stack.pop()
            nxts = adj[v] & remaining
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                remaining ^= nb
                stack.append(bit_index(nb))
    return components


def count_full_endpoint_pairs(adj: list[int], ends: list[int]) -> int:
    """Count unordered endpoint pairs that occur for some Hamiltonian path."""
    n = len(adj)
    full = (1 << n) - 1
    pairs: set[tuple[int, int]] = set()
    # For each possible final endpoint b, peel one edge to learn possible a.
    for b in range(n):
        if not ((ends[full] >> b) & 1):
            continue
        prev_mask = full ^ (1 << b)
        candidates = ends[prev_mask] & adj[b]
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            a = bit_index(bit)
            pairs.add(tuple(sorted((a, b))))
    return len(pairs)


def directed_triangles(out: list[int]) -> int:
    total = 0
    for i, j, k in combinations(range(len(out)), 3):
        if (out[i] >> j) & 1 and (out[j] >> k) & 1 and (out[k] >> i) & 1:
            total += 1
        if (out[i] >> k) & 1 and (out[k] >> j) & 1 and (out[j] >> i) & 1:
            total += 1
    return total


def scc_sizes(out: list[int]) -> tuple[int, ...]:
    n = len(out)
    rev = [0] * n
    for i, bits in enumerate(out):
        while bits:
            bit = bits & -bits
            bits ^= bit
            rev[bit_index(bit)] |= 1 << i

    def reach(start: int, graph: list[int]) -> int:
        seen = 1 << start
        stack = [start]
        while stack:
            v = stack.pop()
            nxts = graph[v] & ~seen
            while nxts:
                bit = nxts & -nxts
                nxts ^= bit
                seen |= bit
                stack.append(bit_index(bit))
        return seen

    remaining = (1 << n) - 1
    sizes: list[int] = []
    while remaining:
        bit = remaining & -remaining
        v = bit_index(bit)
        comp = reach(v, out) & reach(v, rev)
        remaining &= ~comp
        sizes.append(comp.bit_count())
    return tuple(sorted(sizes, reverse=True))


def hp_count_tournament(out: list[int]) -> int:
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last, value in enumerate(dp[mask]):
            if not value:
                continue
            nxts = out[last] & ~mask
            while nxts:
                bit = nxts & -nxts
                nxts ^= bit
                dp[mask | bit][bit_index(bit)] += value
    return sum(dp[-1])


def route_tournament() -> list[int]:
    routes = (
        ("unit-cocyclic 4-set geometry", (5, 5, 5, 3, -1)),
        ("M_L exact extension census", (5, 4, 4, 5, -1)),
        ("endpoint universality audit", (3, 2, 5, 5, -1)),
        ("totally-unfaithful obstruction library", (5, 5, 3, 3, -2)),
        ("degree-deletion core ledger", (3, 3, 4, 5, -1)),
        ("raw graph traceability", (1, 1, 2, 5, -3)),
        ("raw edge count", (1, 1, 1, 5, -4)),
    )
    out = [0] * len(routes)
    for i, j in combinations(range(len(routes)), 2):
        ai = routes[i][1]
        aj = routes[j][1]
        iv = sum(x > y for x, y in zip(ai, aj))
        jv = sum(y > x for x, y in zip(ai, aj))
        if iv > jv or (iv == jv and i < j):
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    print("TOURNAMENT ANALYSIS OVER PROOF ROUTES")
    print("Pairwise observable: which route can still affect the u(22)=60/61 bit after endpoint universality?")
    print("Switch/gauge: proof power, geometry retained, n22 specificity, computability, and risk.")
    order = sorted(range(len(routes)), key=lambda idx: (-out[idx].bit_count(), routes[idx][0]))
    for idx in order:
        print(f"  score={out[idx].bit_count()} {routes[idx][0]}")
    print(f"score_hist={dict(sorted(Counter(bits.bit_count() for bits in out).items()))}")
    print(f"directed_3cycles={directed_triangles(out)} scc={scc_sizes(out)} H={hp_count_tournament(out)}")
    print()
    return out


def audit_core(index: int, code: str) -> dict[str, object]:
    adj = graph_bits(s614.graph6_decode(code))
    n = len(adj)
    full = (1 << n) - 1
    edges = edge_count(adj)
    if edges != U_EXACT_21:
        raise AssertionError((index, edges))
    ends = endpoint_dp(adj)
    full_endpoints = ends[full]
    degrees = [bits.bit_count() for bits in adj]

    deletion_traceable = 0
    endpoint_ears = 0
    all_incident_edges_endpoint_compatible = 0
    exact_edge_deletion_ears = 0
    endpoint_options: list[int] = []
    branch_cuts = 0
    max_components = 0
    delete_edge_deck: Counter[int] = Counter()

    for v in range(n):
        submask = full ^ (1 << v)
        sub_endpoints = ends[submask]
        options = adj[v] & sub_endpoints
        endpoint_options.append(options.bit_count())
        deletion_traceable += int(sub_endpoints != 0)
        endpoint_ears += int(options != 0)
        all_incident_edges_endpoint_compatible += int(options == adj[v])
        exact_edge_deletion_ears += int(options != 0 and edges - degrees[v] == U_EXACT_20)
        comps = connected_components_after_delete(adj, v)
        branch_cuts += int(comps > 2)
        max_components = max(max_components, comps)
        delete_edge_deck[edges - degrees[v]] += 1

    # Graph-only one-vertex extensions over a full traceable core: if every
    # core vertex can be a Hamiltonian endpoint, then any positive-degree new
    # vertex can be attached as a new endpoint.  Degree-4 is the relevant
    # 61-edge extension over a 57-edge core.
    degree4_sets = 0
    degree4_endpoint_compatible = 0
    degree5_sets = 0
    degree5_endpoint_compatible = 0
    for subset in combinations(range(n), 4):
        degree4_sets += 1
        mask = sum(1 << v for v in subset)
        degree4_endpoint_compatible += int((mask & full_endpoints) != 0)
    for subset in combinations(range(n), 5):
        degree5_sets += 1
        mask = sum(1 << v for v in subset)
        degree5_endpoint_compatible += int((mask & full_endpoints) != 0)

    return {
        "core": index,
        "endpoints": full_endpoints.bit_count(),
        "endpoint_pair_lower_bound": count_full_endpoint_pairs(adj, ends),
        "deletion_traceable": deletion_traceable,
        "endpoint_ears": endpoint_ears,
        "all_incident_edges_endpoint_compatible": all_incident_edges_endpoint_compatible,
        "exact_edge_deletion_ears": exact_edge_deletion_ears,
        "endpoint_options_hist": dict(sorted(Counter(endpoint_options).items())),
        "degree_hist": dict(sorted(Counter(degrees).items())),
        "delete_edge_deck": dict(sorted(delete_edge_deck.items())),
        "branch_cuts": branch_cuts,
        "max_components_after_delete": max_components,
        "degree4_extension_sets": degree4_sets,
        "degree4_endpoint_compatible": degree4_endpoint_compatible,
        "degree5_extension_sets": degree5_sets,
        "degree5_endpoint_compatible": degree5_endpoint_compatible,
    }


def main() -> None:
    print("S129 EXACT n=21 GRAPH6 ENDPOINT-UNIVERSALITY AUDIT")
    print("====================================================")
    print()
    print("Question: can the remaining u(22) proof be finished by showing that")
    print("57-edge n=21 cores lose their unit Hamiltonian spine under one-vertex")
    print("extension?  The graph-only answer below is no: the exact cores are")
    print("endpoint-universal.")
    print()

    rows = [audit_core(i, code) for i, code in enumerate(s614.N21_GRAPH6, 1)]
    print("core | endpoints | endpoint-pair lower | deletions traceable | ears | all incident | exact ears | degree hist | endpoint-options hist | branch cuts")
    print("--- | --- | --- | --- | --- | --- | --- | --- | --- | ---")
    for row in rows:
        print(
            f"{row['core']} | {row['endpoints']}/21 | {row['endpoint_pair_lower_bound']} | "
            f"{row['deletion_traceable']}/21 | {row['endpoint_ears']}/21 | "
            f"{row['all_incident_edges_endpoint_compatible']}/21 | "
            f"{row['exact_edge_deletion_ears']} | {row['degree_hist']} | "
            f"{row['endpoint_options_hist']} | {row['branch_cuts']}, maxcomp={row['max_components_after_delete']}"
        )
    print()

    print("GRAPH-ONLY EXTENSION READOUT")
    print("----------------------------")
    for row in rows:
        print(
            f"core {row['core']}: degree-4 neighbor sets endpoint-compatible "
            f"{row['degree4_endpoint_compatible']}/{row['degree4_extension_sets']}; "
            f"degree-5 sets {row['degree5_endpoint_compatible']}/{row['degree5_extension_sets']}"
        )
    print()
    print("Interpretation:")
    print("- In all five exact graph6 cores, every vertex can be a Hamiltonian endpoint.")
    print("- For every existing vertex v, every incident edge can serve as the reattachment")
    print("  edge after deleting v; endpoint-options hist equals degree hist.")
    print("- Therefore any positive-degree graph-only one-vertex extension of a 57-edge")
    print("  core remains traceable by attaching the new vertex to a core endpoint.")
    print("- So traceability/unit-spine arguments cannot prove u(22)<=60 by themselves.")
    print("- The proof must forbid a realizable unit-cocyclic degree-4 neighbor set over")
    print("  the five 57-edge cores, or else use the totally-unfaithful obstruction library.")
    print()

    route_tournament()
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertices considered: graph vertices, deletion ears, endpoint masks,")
    print("4-subsets of core vertices, candidate extension centres, unit-cocyclic sets,")
    print("and totally-unfaithful obstruction packets.  The graph quotient preserves")
    print("traceability and endpoint universality but destroys planar realization and")
    print("circle-centre/cocyclicity constraints.  That destroyed geometry is exactly")
    print("where the remaining u(22) proof has to live.")


if __name__ == "__main__":
    main()
