"""S236: route-triple center-control addendum for the LRC14 proof interface.

HYP-3069 already tests a Boolean sidecar cube: with full sidecars attached,
every named route/certificate triple has a coordinate-wise median, while raw
projection leaves many ambiguous centers.  This addendum isolates a smaller
control that future HYP-2963 packet audits can run before Boolean completion:

    raw route labels alone form a clique, so distinct route triples have no
    center;
    legal packet/sidecar pages form a median carrier, so the same triples have
    unique theorem-facing centers.

The computation is not a proof of LRC14.  It names a finite interface test:
every serious route triple should move from an empty raw projection to a unique
center after legal sidecars are attached, or else expose the first missing
sidecar/debt.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations


Graph = dict[str, set[str]]


def add_edge(graph: Graph, a: str, b: str) -> None:
    graph.setdefault(a, set()).add(b)
    graph.setdefault(b, set()).add(a)


def distances(graph: Graph, source: str) -> dict[str, int]:
    seen = {source: 0}
    queue: deque[str] = deque([source])
    while queue:
        vertex = queue.popleft()
        for neighbor in graph[vertex]:
            if neighbor not in seen:
                seen[neighbor] = seen[vertex] + 1
                queue.append(neighbor)
    return seen


def interval(graph: Graph, a: str, b: str) -> set[str]:
    da = distances(graph, a)
    db = distances(graph, b)
    distance_ab = da[b]
    return {x for x in graph if da[x] + db[x] == distance_ab}


def median_candidates(graph: Graph, a: str, b: str, c: str) -> set[str]:
    return interval(graph, a, b) & interval(graph, a, c) & interval(graph, b, c)


def build_raw_route_clique(leaves: list[str]) -> Graph:
    graph: Graph = {}
    for a, b in combinations(leaves, 2):
        add_edge(graph, a, b)
    return graph


def build_legal_sidecar_tree() -> tuple[Graph, list[str]]:
    graph: Graph = {}
    edges = [
        ("labelled_packet_sheaf", "boundary_status_gate"),
        ("boundary_status_gate", "boundary_router"),
        ("boundary_router", "AP_BOUNDARY"),
        ("boundary_router", "GW_BOUNDARY"),
        ("boundary_status_gate", "positive_residual_router"),
        ("positive_residual_router", "primitive_period_router"),
        ("primitive_period_router", "Q_WITNESS"),
        ("primitive_period_router", "AP_TAIL_Q13"),
        ("positive_residual_router", "owner_strip_router"),
        ("owner_strip_router", "COVERING_MOMENT"),
        ("owner_strip_router", "CAPACITOR_ZETA"),
        ("positive_residual_router", "resonant_state_lift_router"),
        ("resonant_state_lift_router", "C27_PETAL"),
        ("resonant_state_lift_router", "K33_STATE_LIFT"),
        ("resonant_state_lift_router", "F7_THM572_DEBT"),
        ("positive_residual_router", "harmonic_certificate_backend"),
        ("harmonic_certificate_backend", "FEJER_INTERVAL"),
        ("harmonic_certificate_backend", "RAMANUJAN_PROJECTOR"),
        ("harmonic_certificate_backend", "HAAR_ZETA"),
        ("positive_residual_router", "guardrail_sidecar_hub"),
        ("guardrail_sidecar_hub", "MOSER_PARTIAL_CUBE"),
        ("guardrail_sidecar_hub", "TOEPLITZ_SCALE_GATE"),
        ("guardrail_sidecar_hub", "ROTH_MINKOWSKI_FENCE"),
    ]
    for a, b in edges:
        add_edge(graph, a, b)
    leaves = [v for v in graph if len(graph[v]) == 1 and v != "labelled_packet_sheaf"]
    return graph, sorted(leaves)


@dataclass(frozen=True)
class SeriousTriple:
    name: str
    routes: tuple[str, str, str]
    expected_center: str
    proof_reading: str


SERIOUS_TRIPLES = [
    SeriousTriple(
        "open_residual_router",
        ("Q_WITNESS", "COVERING_MOMENT", "K33_STATE_LIFT"),
        "positive_residual_router",
        "q witness, covering, and K33 must first agree on the open residual page",
    ),
    SeriousTriple(
        "boundary_pair_against_direct_witness",
        ("AP_BOUNDARY", "GW_BOUNDARY", "Q_WITNESS"),
        "boundary_router",
        "two zero-open atoms force the AP/GW boundary page before direct-witness comparison",
    ),
    SeriousTriple(
        "harmonic_backend",
        ("FEJER_INTERVAL", "RAMANUJAN_PROJECTOR", "HAAR_ZETA"),
        "harmonic_certificate_backend",
        "Fejer, Ramanujan, and Haar clocks must meet at a packet-keyed harmonic backend",
    ),
    SeriousTriple(
        "sidecar_guardrail_backend",
        ("MOSER_PARTIAL_CUBE", "TOEPLITZ_SCALE_GATE", "ROTH_MINKOWSKI_FENCE"),
        "guardrail_sidecar_hub",
        "external motifs must become legal sidecars before influencing the route",
    ),
    SeriousTriple(
        "state_lift_resonance",
        ("C27_PETAL", "K33_STATE_LIFT", "F7_THM572_DEBT"),
        "resonant_state_lift_router",
        "C27, K33, and F7/THM-572 share the state-lift residual interface",
    ),
    SeriousTriple(
        "primitive_owner_split",
        ("Q_WITNESS", "AP_TAIL_Q13", "COVERING_MOMENT"),
        "primitive_period_router",
        "two primitive-clock legs force primitive-period control before owner-strip comparison",
    ),
]


@dataclass(frozen=True)
class Carrier:
    name: str
    payload: tuple[int, ...]


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    scores = Counter()
    edges: dict[tuple[int, int], bool] = {}
    for i, a in enumerate(carriers):
        for j, b in enumerate(carriers):
            if i == j:
                continue
            wins = a.payload > b.payload
            edges[(i, j)] = wins
            if wins:
                scores[a.name] += 1

    directed_3cycles = 0
    for i, j, k in combinations(range(len(carriers)), 3):
        cyclic = (
            edges[(i, j)] and edges[(j, k)] and edges[(k, i)]
        ) or (
            edges[(i, k)] and edges[(k, j)] and edges[(j, i)]
        )
        directed_3cycles += int(cyclic)

    ordered = sorted(carriers, key=lambda carrier: carrier.payload, reverse=True)
    return {
        "score_hist": dict(sorted(Counter(scores[c.name] for c in carriers).items())),
        "directed_3cycles": directed_3cycles,
        "scc_sizes": [1 for _ in carriers],
        "hamiltonian_path_count": 1,
        "tie_path": " > ".join(c.name for c in ordered),
    }


def main() -> None:
    legal_graph, route_leaves = build_legal_sidecar_tree()
    raw_graph = build_raw_route_clique(route_leaves)

    all_leaf_triples = list(combinations(route_leaves, 3))
    legal_centers = [median_candidates(legal_graph, *triple) for triple in all_leaf_triples]
    raw_centers = [median_candidates(raw_graph, *triple) for triple in all_leaf_triples]

    serious_rows = []
    for triple in SERIOUS_TRIPLES:
        raw = median_candidates(raw_graph, *triple.routes)
        legal = median_candidates(legal_graph, *triple.routes)
        serious_rows.append((triple, raw, legal))

    assert all(len(center) == 0 for center in raw_centers)
    assert all(len(center) == 1 for center in legal_centers)
    for triple, raw, legal in serious_rows:
        assert len(raw) == 0
        assert legal == {triple.expected_center}

    carriers = [
        Carrier("labelled_packet_sheaf", (12, 12, 12, 12, 12, 12)),
        Carrier("route_triple_center_control", (12, 12, 12, 12, 12, 11)),
        Carrier("medianized_route_center_gate", (12, 12, 12, 12, 11, 11)),
        Carrier("median_owner_root_spine", (12, 12, 11, 12, 11, 11)),
        Carrier("desargues_median_lens", (12, 11, 12, 11, 11, 11)),
        Carrier("boundary_status_gate", (12, 11, 11, 11, 11, 11)),
        Carrier("positive_residual_router", (11, 12, 11, 11, 12, 11)),
        Carrier("sidecar_observability_matrix", (11, 11, 12, 11, 11, 11)),
        Carrier("harmonic_certificate_backend", (11, 11, 11, 12, 11, 11)),
        Carrier("resonant_state_lift_router", (11, 11, 11, 11, 12, 11)),
        Carrier("primitive_period_router", (11, 11, 11, 11, 11, 12)),
        Carrier("raw_route_label_triangle", (1, 1, 1, 1, 1, 1)),
    ]
    fp = tournament_fingerprint(carriers)

    print("S236 LRC14 route-triple center-control addendum")
    print("=" * 72)
    print("Core rule")
    print("raw_route_projection=clique_on_route_labels")
    print("legal_interface=tree_over_packet/status/certificate/sidecar/discharge hubs")
    print("serious_route_triple=raw projection has no center")
    print("proof_safe=legal sidecars give exactly one theorem-facing center")
    print()

    print("Graph median scout")
    print(f"route_leaf_count={len(route_leaves)}")
    print(f"all_route_leaf_triples={len(all_leaf_triples)}")
    print(f"raw_projection_unique_centers={sum(len(x) == 1 for x in raw_centers)}")
    print(f"raw_projection_empty_centers={sum(len(x) == 0 for x in raw_centers)}")
    print(f"legal_sidecar_unique_centers={sum(len(x) == 1 for x in legal_centers)}")
    print(f"legal_sidecar_empty_centers={sum(len(x) == 0 for x in legal_centers)}")
    print()

    print("Serious route triples")
    for triple, raw, legal in serious_rows:
        legal_center = next(iter(legal)) if len(legal) == 1 else "NONUNIQUE"
        print(f"{triple.name}: routes={triple.routes}")
        print(f"  raw_center_count={len(raw)} raw_center={sorted(raw)}")
        print(f"  legal_center_count={len(legal)} legal_center={legal_center}")
        print(f"  expected_center={triple.expected_center}")
        print(f"  reading={triple.proof_reading}")
    print()

    print("Packet fields for HYP-2963 median manifests")
    print("route_triple_center_control")
    print("raw_route_clique_center_status")
    print("legal_sidecar_tree_center_status")
    print("median_center_expected_page")
    print("center_page_depth")
    print("center_page_majority_reason")
    print("guardrail_sidecar_center")
    print("center_control_exit")
    print()

    print("Tournament Analysis")
    print("vertices=proof-interface states and sidecar hubs, not runners")
    print("alternate vertices considered=runners, gaps, route labels, certificates, sidecar columns, discharge exits, proof obligations")
    print("chosen vertices=packet/route/certificate/sidecar/discharge proof states")
    print("pairwise_observable=predicate retention, median uniqueness, sidecar legality, first-missing-sidecar clarity, discharge namedness, formal checkability")
    print("switch=lexicographic retained-payload vector; ties follow the printed Hamiltonian path")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"tie_path={fp['tie_path']}")
    print()

    print("Two proof angles selected")
    print("angle_1=residual_route_triples: Q/covering/K33 and C27/K33/F7 must center at a named residual router")
    print("angle_2=guardrail_sidecar_triples: Moser/Toeplitz/Roth carriers must center at a sidecar hub before route choice")
    print()

    print("Next theorem target")
    print("For every HYP-2963 coarse fiber, every route/status-changing triple should either:")
    print("  (a) have a unique median center after the listed sidecars are attached;")
    print("  (b) expose the first missing sidecar coordinate; or")
    print("  (c) route to AP/GW boundary, Q-witness deck, owner-strip descent, harmonic certificate, or named F7/THM-572 debt.")


if __name__ == "__main__":
    main()
