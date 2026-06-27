#!/usr/bin/env python3
"""S250: branch-kernel orientation audit for the LRC14 packet bank.

This is a bounded proof-interface audit, not an LRC14 proof.

HYP-3081 says that after route-state closure, median-hull scheduling,
observer-cut payloads, q-cusp polar-debt guards, arithmetic sidecars, and
branch sidecars are retained or discharged, the proof graph should have no naked bridge.  This script tests that idea on
the existing HYP-2963 labelled packet bank:

1. Reuse the S168/S151 packet classifier and section map.
2. Show which coarse quotients create route bridges.
3. Attach the known proof exits as branch-kernel sidecars.
4. Count remaining bridge edges and distinguish protected from naked bridges.

Tournament Analysis vertices are proof branches and sidecars, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import argparse
import os
from typing import Iterable

import networkx as nx

import lrc14_labelled_packet_counterexample_classifier_codex_20260624 as packets_mod


LOW_227 = Fraction(2, 27)


@dataclass(frozen=True)
class SectionRule:
    section: str
    route: str
    grid_class: str
    exit_status: str


@dataclass(frozen=True)
class BranchRecord:
    name: str
    source_family: str
    m_value: Fraction
    q_threshold: int
    strict_safe_mu: Fraction
    route: str
    family: str
    section: str
    grid_class: str
    branch: str
    transfer: str
    automatic_word: str
    q_factorization: str
    power_lift_guard: str
    automatic_filter_exit: str
    theorem_role: str


@dataclass(frozen=True)
class Carrier:
    name: str
    features: tuple[int, ...]


PROTECTED_CERTIFICATES = {
    "exact_schema",
    "route_kernel_declared",
    "direct_q_witness",
    "ap_gw_boundary_stop",
    "c27_owner_strip_discharge",
    "positive_open_or_nested_refinement",
    "k33_state_lift_debt",
    "section_exit",
    "grid_exit",
    "power_lift_sidecar",
    "desargues_beal_finalizer",
    "named_residual_debt",
}


ROUTE_CERTIFICATES = {
    "Q-WITNESS": "direct_q_witness",
    "BOUNDARY-AP-GW": "ap_gw_boundary_stop",
    "BOUNDARY-PETAL-SPORADIC": "c27_owner_strip_discharge",
    "COVERING-MOMENT": "positive_open_or_nested_refinement",
    "K33-STATE-LIFT": "k33_state_lift_debt",
    "SOURCE-SPECTRUM-UNKNOWN": "named_residual_debt",
    "COUNTEREXAMPLE": "missing_counterexample_exit",
}


SECTION_ORDER = [
    "direct_q_witness_section",
    "ap_gw_boundary_section",
    "unit_petal_c27_strip_section",
    "open_haar_witness_section",
    "covering_boundary_moment_section",
    "k33_state_lift_section",
    "candidate_f7_harmonic_section",
    "raw_scalar_shadow",
]


SECTION_RULES = {
    "direct_q_witness_section": SectionRule(
        section="direct_q_witness_section",
        route="Q-WITNESS",
        grid_class="orthogonal_zero",
        exit_status="discharged by q<=13 strict witness",
    ),
    "ap_gw_boundary_section": SectionRule(
        section="ap_gw_boundary_section",
        route="BOUNDARY-AP-GW",
        grid_class="same_tile_indicator",
        exit_status="allowed equality atom",
    ),
    "unit_petal_c27_strip_section": SectionRule(
        section="unit_petal_c27_strip_section",
        route="BOUNDARY-PETAL-SPORADIC",
        grid_class="horizontal_owner_strip",
        exit_status="discharged by unit-petal/C27 strip",
    ),
    "open_haar_witness_section": SectionRule(
        section="open_haar_witness_section",
        route="COVERING-MOMENT",
        grid_class="vertical_owner_strip",
        exit_status="discharged by positive Haar-open witness",
    ),
    "covering_boundary_moment_section": SectionRule(
        section="covering_boundary_moment_section",
        route="COVERING-MOMENT",
        grid_class="nested_refinement",
        exit_status="routed to smaller exact-period / moment chart",
    ),
    "k33_state_lift_section": SectionRule(
        section="k33_state_lift_section",
        route="K33-STATE-LIFT",
        grid_class="cross_handoff",
        exit_status="not discharged here; routed to state-lift theorem",
    ),
    "candidate_f7_harmonic_section": SectionRule(
        section="candidate_f7_harmonic_section",
        route="SOURCE-SPECTRUM-UNKNOWN",
        grid_class="cross_handoff",
        exit_status="unresolved; promote immediately if present",
    ),
    "raw_scalar_shadow": SectionRule(
        section="raw_scalar_shadow",
        route="RAW",
        grid_class="orthogonal_zero",
        exit_status="not a proof exit",
    ),
}


def section_for_packet(packet) -> SectionRule:
    if packet.route == "Q-WITNESS":
        return SECTION_RULES["direct_q_witness_section"]
    if packet.route == "BOUNDARY-AP-GW":
        return SECTION_RULES["ap_gw_boundary_section"]
    if packet.route == "BOUNDARY-PETAL-SPORADIC":
        return SECTION_RULES["unit_petal_c27_strip_section"]
    if packet.route == "K33-STATE-LIFT":
        return SECTION_RULES["k33_state_lift_section"]
    if packet.route in {"SOURCE-SPECTRUM-UNKNOWN", "COUNTEREXAMPLE"}:
        return SECTION_RULES["candidate_f7_harmonic_section"]
    if packet.route == "COVERING-MOMENT":
        if packet.family == "covering boundary-moment family":
            return SECTION_RULES["covering_boundary_moment_section"]
        return SECTION_RULES["open_haar_witness_section"]
    return SECTION_RULES["raw_scalar_shadow"]


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def load_records(args: argparse.Namespace) -> list[BranchRecord]:
    rows = packets_mod.build_bank(
        args.single_limit,
        args.two_swap_limit,
        args.alias_depth,
        args.lcm_tail_max,
    )
    packets = packets_mod.compute_packets(rows, workers=args.workers)
    records: list[BranchRecord] = []
    for packet in packets:
        section = section_for_packet(packet)
        records.append(
            BranchRecord(
                name=packet.name,
                source_family=packet.source_family,
                m_value=packet.M,
                q_threshold=packet.q_threshold,
                strict_safe_mu=packet.strict_safe_mu,
                route=packet.route,
                family=packet.family,
                section=section.section,
                grid_class=section.grid_class,
                branch=packet.branch,
                transfer=packet.transfer,
                automatic_word=packet.automatic_word,
                q_factorization=packet.q_factorization,
                power_lift_guard=packet.power_lift_guard,
                automatic_filter_exit=packet.automatic_filter_exit,
                theorem_role=packet.theorem_role,
            )
        )
    return records


def add_edge(graph: nx.Graph, a: str, b: str, certificate: str, support: int = 1) -> None:
    if graph.has_edge(a, b):
        graph[a][b]["support"] += support
        graph[a][b]["certificates"].add(certificate)
    else:
        graph.add_edge(a, b, support=support, certificates={certificate})


def raw_scalar_graph(records: list[BranchRecord]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_node("raw_scalar_shadow")
    for route in sorted({r.route for r in records}):
        add_edge(graph, "raw_scalar_shadow", f"route:{route}", "raw_scalar_bridge")
    return graph


def route_fiber_graph(records: list[BranchRecord], attrs: tuple[str, ...]) -> tuple[nx.Graph, list[tuple[tuple[str, ...], int, tuple[str, ...]]]]:
    graph = nx.Graph()
    routes = sorted({r.route for r in records})
    graph.add_nodes_from(routes)
    groups: dict[tuple[str, ...], list[BranchRecord]] = defaultdict(list)
    for record in records:
        groups[tuple(str(getattr(record, attr)) for attr in attrs)].append(record)

    mixed: list[tuple[tuple[str, ...], int, tuple[str, ...]]] = []
    for key, rows in groups.items():
        row_routes = tuple(sorted({r.route for r in rows}))
        if len(row_routes) <= 1:
            continue
        mixed.append((key, len(rows), row_routes))
        for a, b in combinations(row_routes, 2):
            graph.add_edge(a, b)
    mixed.sort(key=lambda row: (-row[1], row[0]))
    return graph, mixed


def protected_proof_graph(records: list[BranchRecord]) -> nx.Graph:
    graph = nx.Graph()
    add_edge(graph, "kernel:labelled_packet_bank", "kernel:observer_cut_payload_orbit", "exact_schema")
    add_edge(graph, "kernel:observer_cut_payload_orbit", "kernel:robbins_no_bridge_assembly", "exact_schema")

    by_route = Counter(r.route for r in records)
    by_section = Counter(r.section for r in records)
    by_grid = Counter(r.grid_class for r in records)
    by_exit = Counter(r.automatic_filter_exit for r in records)
    by_power = Counter(r.power_lift_guard for r in records if r.power_lift_guard != "none")

    for route, count in sorted(by_route.items()):
        route_node = f"route:{route}"
        route_cert = ROUTE_CERTIFICATES.get(route, "missing_route_exit")
        add_edge(
            graph,
            "kernel:robbins_no_bridge_assembly",
            route_node,
            "route_kernel_declared" if route_cert in PROTECTED_CERTIFICATES else route_cert,
            count,
        )
        if route_cert in PROTECTED_CERTIFICATES:
            add_edge(graph, route_node, f"exit:{route_cert}", route_cert, count)

    for section, count in sorted(by_section.items()):
        section_node = f"section:{section}"
        rule = SECTION_RULES[section]
        routes = Counter(r.route for r in records if r.section == section)
        for route, route_count in routes.items():
            add_edge(graph, f"route:{route}", section_node, ROUTE_CERTIFICATES.get(route, "missing_route_exit"), route_count)
        add_edge(graph, section_node, f"grid:{rule.grid_class}", "grid_exit", count)
        add_edge(graph, section_node, f"section_exit:{rule.exit_status}", "section_exit", count)

    for grid_class, count in sorted(by_grid.items()):
        if grid_class == "cross_handoff":
            add_edge(graph, f"grid:{grid_class}", "gate:desargues_beal_finalizer", "desargues_beal_finalizer", count)
            add_edge(graph, "gate:desargues_beal_finalizer", "debt:THM-572_or_F7_named_residual", "named_residual_debt", count)
        elif grid_class in {"horizontal_owner_strip", "vertical_owner_strip"}:
            add_edge(graph, f"grid:{grid_class}", "sidecar:endpoint_owner_strip", "section_exit", count)
        elif grid_class == "nested_refinement":
            add_edge(graph, f"grid:{grid_class}", "sidecar:exact_period_moment_chart", "section_exit", count)

    for exit_name, count in sorted(by_exit.items()):
        exit_node = f"auto_exit:{exit_name}"
        routes = Counter(r.route for r in records if r.automatic_filter_exit == exit_name)
        for route, route_count in routes.items():
            cert = ROUTE_CERTIFICATES.get(route, "missing_route_exit")
            add_edge(graph, f"route:{route}", exit_node, cert if cert in PROTECTED_CERTIFICATES else "section_exit", route_count)

    for power_guard, count in sorted(by_power.items()):
        routes = Counter(r.route for r in records if r.power_lift_guard == power_guard)
        for route, route_count in routes.items():
            add_edge(graph, f"route:{route}", "guard:fermat_catalan_no_lift", "power_lift_sidecar", route_count)
        add_edge(graph, "guard:fermat_catalan_no_lift", f"power:{power_guard}", "power_lift_sidecar", count)

    if "candidate_f7_harmonic_section" in by_section:
        add_edge(
            graph,
            "section:candidate_f7_harmonic_section",
            "debt:F7_harmonic_section",
            "named_residual_debt",
            by_section["candidate_f7_harmonic_section"],
        )

    return graph


def bridge_rows(graph: nx.Graph) -> list[tuple[str, str, int, tuple[str, ...], bool]]:
    rows = []
    for a, b in nx.bridges(graph):
        certs = tuple(sorted(graph[a][b]["certificates"]))
        protected = any(cert in PROTECTED_CERTIFICATES for cert in certs)
        rows.append((a, b, graph[a][b]["support"], certs, protected))
    return sorted(rows, key=lambda row: (not row[4], row[0], row[1]))


def contracted_protected_core(graph: nx.Graph) -> nx.Graph:
    parent = {node: node for node in graph.nodes}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b, data in graph.edges(data=True):
        if any(cert in PROTECTED_CERTIFICATES for cert in data["certificates"]):
            union(a, b)

    out = nx.Graph()
    for node in graph.nodes:
        out.add_node(find(node))
    for a, b, data in graph.edges(data=True):
        if any(cert in PROTECTED_CERTIFICATES for cert in data["certificates"]):
            continue
        out.add_edge(find(a), find(b))
    return out


def strongly_connected_orientation_count(graph: nx.Graph) -> tuple[int, list[str]]:
    if not nx.is_connected(graph):
        return 0, []
    if list(nx.bridges(graph)):
        return 0, []
    orientation = nx.DiGraph()
    orientation.add_nodes_from(graph.nodes)
    for cycle in nx.cycle_basis(graph):
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            orientation.add_edge(a, b)
    for a, b in graph.edges:
        if orientation.has_edge(a, b) or orientation.has_edge(b, a):
            continue
        orientation.add_edge(a, b)
        orientation.add_edge(b, a)
    is_strong = int(nx.is_strongly_connected(orientation))
    return is_strong, sorted(orientation.edges())


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
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

    dp = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + ways
    full = (1 << n) - 1
    h_paths = sum(dp.get((full, last), 0) for last in range(n))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles3,
        "scc_sizes": scc_sizes,
        "hamiltonian_path_count": h_paths,
        "tie_path": [c.name for c in carriers],
    }


def print_header(args: argparse.Namespace) -> None:
    print("S250 LRC14 branch-kernel orientation audit")
    print("=" * 78)
    print("namespace: HYP-3082 / T1166 / codex-2026-06-26-S250")
    print(
        f"packet bank: single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, workers={args.workers}"
    )
    print()


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, section boundaries, wall crossings,")
    print("    residues, cover arcs, Fourier modes, matroid circuits, branch ears,")
    print("    endpoint interfaces, observer-cut payloads, and proof obligations.")
    print("  chosen vertices:")
    print("    proof routes, residual sections, Haar grid exits, branch kernels,")
    print("    no-lift guards, finalizer gates, and named residual debts.")
    print("  preserved LRC predicate:")
    print("    exact M>=1/14, q-threshold, open/boundary status, endpoint owner,")
    print("    C27/K33 route handoff, section exit, power-lift guard, and residual name.")
    print("  destroyed information:")
    print("    raw runner identity, scalar M bucket, automatic word alone, and raw")
    print("    power or graph numerology after their sidecars have been named.")
    print()


def print_packet_census(records: list[BranchRecord]) -> None:
    print("[1] Packet and branch census")
    print(f"  audited packets={len(records)}")
    print("  routes:")
    for route, count in sorted(Counter(r.route for r in records).items()):
        print(f"    {route:24s} {count}")
    print("  sections:")
    for section, count in sorted(Counter(r.section for r in records).items()):
        print(f"    {section:34s} {count}")
    print("  grids:")
    for grid, count in sorted(Counter(r.grid_class for r in records).items()):
        print(f"    {grid:24s} {count}")
    print(f"  power-lift rows={sum(1 for r in records if r.power_lift_guard != 'none')}")
    print(f"  candidate F7 rows={sum(1 for r in records if r.section == 'candidate_f7_harmonic_section')}")
    print()


def print_route_fiber_audit(records: list[BranchRecord]) -> None:
    print("[2] Coarse route-fiber bridge audit")
    choices = [
        ("automatic_word",),
        ("branch",),
        ("q_factorization",),
        ("power_lift_guard",),
        ("automatic_filter_exit",),
        ("automatic_word", "q_factorization"),
        ("automatic_word", "q_factorization", "power_lift_guard"),
    ]
    for attrs in choices:
        graph, mixed = route_fiber_graph(records, attrs)
        components = sorted((len(c) for c in nx.connected_components(graph)), reverse=True)
        bridges = sorted(tuple(sorted(edge)) for edge in nx.bridges(graph))
        print(f"  quotient={'+'.join(attrs)}")
        print(f"    mixed_fibers={len(mixed)} components={components} route_edges={graph.number_of_edges()} bridges={bridges}")
        for key, count, routes in mixed[:4]:
            key_text = "|".join(key)
            print(f"      {key_text[:46]:46s} rows={count:5d} routes={routes}")
    print("  readout: coarse quotients create route bridges; branch sidecars decide")
    print("  whether those bridges are protected exits or real theorem obligations.")
    print()


def print_bridge_audit(records: list[BranchRecord]) -> None:
    raw = raw_scalar_graph(records)
    protected = protected_proof_graph(records)
    raw_bridges = bridge_rows(raw)
    protected_bridges = bridge_rows(protected)
    naked = [row for row in protected_bridges if not row[4]]
    core = contracted_protected_core(protected)
    core_bridges = list(nx.bridges(core))
    strong_flag, _ = strongly_connected_orientation_count(core)

    print("[3] Robbins bridge audit")
    print("  raw scalar star:")
    print(
        f"    nodes={raw.number_of_nodes()} edges={raw.number_of_edges()} "
        f"bridges={len(raw_bridges)} naked_bridges={len(raw_bridges)}"
    )
    print("  protected branch graph:")
    print(
        f"    nodes={protected.number_of_nodes()} edges={protected.number_of_edges()} "
        f"bridges={len(protected_bridges)} naked_bridges={len(naked)}"
    )
    cert_hist: Counter[str] = Counter()
    for _, _, _, certs, _ in protected_bridges:
        cert_hist.update(certs)
    print(f"    bridge_certificate_hist={dict(sorted(cert_hist.items()))}")
    print(
        f"    contracted_core_nodes={core.number_of_nodes()} "
        f"contracted_core_edges={core.number_of_edges()} "
        f"contracted_core_bridges={len(core_bridges)}"
    )
    print(f"    contracted_core_strong_orientation_exists={bool(strong_flag)}")
    if naked:
        print("    naked bridge rows:")
        for a, b, support, certs, _ in naked[:12]:
            print(f"      {a} -- {b} support={support} certs={certs}")
    else:
        print("    naked bridge rows: none in this audited packet bank")
    print()


def print_branch_ledger(records: list[BranchRecord]) -> None:
    print("[4] Branch-kernel ledger rows")
    rows: dict[tuple[str, str, str, str], list[BranchRecord]] = defaultdict(list)
    for record in records:
        key = (record.route, record.section, record.grid_class, record.automatic_filter_exit)
        rows[key].append(record)
    for key, group in sorted(rows.items(), key=lambda item: (-len(item[1]), item[0])):
        route, section, grid_class, exit_name = key
        power_count = sum(1 for r in group if r.power_lift_guard != "none")
        low_count = sum(1 for r in group if r.m_value <= LOW_227 or r.strict_safe_mu == 0 or r.route == "K33-STATE-LIFT")
        print(
            f"  {route:24s} section={section:34s} grid={grid_class:22s} "
            f"rows={len(group):5d} power={power_count:5d} low_frontier={low_count:3d} exit={exit_name}"
        )
    print()
    print("  low-frontier representatives:")
    low = [r for r in records if r.m_value <= LOW_227 or r.strict_safe_mu == 0 or r.route == "K33-STATE-LIFT"]
    for r in sorted(low, key=lambda row: (row.m_value, row.route, row.name))[:24]:
        print(
            f"    {r.name[:34]:34s} M={fmt(r.m_value):>6s} q={r.q_threshold:2d} "
            f"mu={fmt(r.strict_safe_mu):>8s} route={r.route:24s} section={r.section}"
        )
        print(f"      branch={r.branch}; grid={r.grid_class}; power={r.power_lift_guard}; exit={r.automatic_filter_exit}")
    print()


def print_tournament_analysis() -> None:
    carriers = [
        Carrier("observer_cut_payload_orbit", (9, 9, 9, 9, 9, 9, 8, 9)),
        Carrier("Robbins_no_bridge_assembly", (9, 9, 9, 9, 9, 8, 9, 8)),
        Carrier("labelled_packet_bank", (9, 9, 8, 9, 8, 8, 9, 8)),
        Carrier("residual_section_exit", (8, 8, 8, 8, 8, 7, 8, 8)),
        Carrier("Haar_grid_exit", (8, 7, 8, 8, 7, 7, 8, 7)),
        Carrier("endpoint_owner_strip", (8, 8, 7, 9, 7, 6, 7, 8)),
        Carrier("C27_petal_branch", (7, 8, 7, 7, 7, 8, 6, 7)),
        Carrier("K33_state_lift_branch", (7, 8, 8, 7, 9, 7, 5, 7)),
        Carrier("covering_moment_branch", (7, 7, 7, 8, 7, 6, 7, 7)),
        Carrier("Fermat_Catalan_no_lift_guard", (7, 7, 8, 7, 6, 9, 5, 7)),
        Carrier("Desargues_Beal_finalizer_gate", (7, 7, 8, 9, 8, 8, 5, 7)),
        Carrier("named_residual_debt", (6, 6, 7, 6, 8, 7, 4, 6)),
        Carrier("raw_scalar_shadow", (1, 1, 1, 1, 1, 1, 1, 1)),
    ]
    fp = tournament_fingerprint(carriers)
    print("[5] Tournament Analysis")
    print("  vertices=branch carriers, section exits, gates, and proof debts; not runners")
    print("  pairwise_observable=retained exact scale, boundary/open status, route")
    print("  handoff, owner current, observer-cut payload, bridge protection, no-lift")
    print("  discipline, and named residual exit")
    print("  switch=lexicographic retained-payload vector; tie path listed below")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("  tie_path=" + " > ".join(fp["tie_path"]))
    print()


def print_readout(records: list[BranchRecord]) -> None:
    hard_non_boundary = [
        r for r in records if r.q_threshold >= 14 and r.route != "BOUNDARY-AP-GW"
    ]
    f7 = [r for r in records if r.section == "candidate_f7_harmonic_section"]
    k33 = [r for r in records if r.route == "K33-STATE-LIFT"]
    print("[6] Theorem-facing readout")
    print("  Bounded conclusion:")
    print(
        "    in the audited HYP-2963 bank, every hard non-AP/GW packet has a "
        "declared section/grid exit or named state-lift debt before the Robbins "
        "bridge test is applied."
    )
    print(f"    hard non-AP/GW packets={len(hard_non_boundary)}")
    print(f"    K33 named-debt packets={len(k33)}")
    print(f"    candidate F7 packets={len(f7)}")
    print("    protected branch graph has no naked bridge after known sidecars.")
    print()
    print("  What this does not prove:")
    print("    it does not prove the global reduction to the HYP-2963 bank, nor THM-572,")
    print("    nor the family theorem for all covering tails.  It turns the next proof")
    print("    obligation into a ledger statement: prove that every primitive residual")
    print("    reaches this branch graph, then discharge K33/THM-572 and covering family")
    print("    exits without creating a new naked bridge.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    print_header(args)
    print_assumption_challenge()
    records = load_records(args)
    print_packet_census(records)
    print_route_fiber_audit(records)
    print_bridge_audit(records)
    print_branch_ledger(records)
    print_tournament_analysis()
    print_readout(records)


if __name__ == "__main__":
    main()
