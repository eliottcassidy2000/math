#!/usr/bin/env python3
"""LRC14 colored-extension carrier scout.

This is a synthesis scout, not a theorem prover.  It connects older coloring
work (HYP-2594/2595, HYP-2247, HYP-2250) with the current gate/floor frontier
(HYP-3438, HYP-3455, HYP-3456, HYP-3457, HYP-3458, HYP-3459, HYP-3460) and
the observer-extension guardrails (HYP-3056, HYP-3133, HYP-3134).

Tournament Analysis:
  vertices: proof carriers and quotient shadows;
  pairwise observable: retained LRC predicate payload under coloring plus
    extension/gluing;
  switch: higher weighted carrier score, tie-broken by proof-facing order;
  tie Hamiltonian path: declared by the chosen proof-facing order.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


WEIGHTS = {
    "predicate": 5,
    "color_payload": 4,
    "extension_payload": 5,
    "endpoint_owner": 4,
    "finite_audit": 3,
    "frontier_fit": 5,
    "scalar_guardrail": 4,
    "discharge": 5,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    refs: tuple[str, ...]
    scores: dict[str, int]
    preserves: str
    destroys: str
    proof_pull: str
    challenged_assumption: str

    def weighted_score(self) -> int:
        return sum(self.scores[key] * WEIGHTS[key] for key in WEIGHTS)

    def score_tuple(self) -> tuple[int, ...]:
        return tuple(self.scores[key] for key in WEIGHTS)


def carrier_bank() -> list[Carrier]:
    return [
        Carrier(
            name="random031_mirror_colored_extension_clause",
            refs=("HYP-3455", "HYP-3438", "HYP-3453"),
            scores={
                "predicate": 5,
                "color_payload": 4,
                "extension_payload": 5,
                "endpoint_owner": 5,
                "finite_audit": 5,
                "frontier_fit": 5,
                "scalar_guardrail": 5,
                "discharge": 5,
            },
            preserves=(
                "finite two-color gate-gluing predicate with branch mask, "
                "B-S-B-S-B word, owner deltas, mirror cut, and low-rank escape"
            ),
            destroys="unrestricted row family generality",
            proof_pull=(
                "prove the seven-owner mirror pair cannot extend to a full "
                "cover without one of the named HYP-3455 exits"
            ),
            challenged_assumption=(
                "a rank-6 noncanonical obstruction is an uncontrolled family"
            ),
        ),
        Carrier(
            name="AP84_colored_endpoint_floor_packet",
            refs=(
                "HYP-3454",
                "HYP-3456",
                "HYP-3457",
                "HYP-3458",
                "HYP-3459",
                "HYP-3460",
                "HYP-3431",
            ),
            scores={
                "predicate": 5,
                "color_payload": 4,
                "extension_payload": 5,
                "endpoint_owner": 5,
                "finite_audit": 5,
                "frontier_fit": 5,
                "scalar_guardrail": 5,
                "discharge": 4,
            },
            preserves=(
                "AP84 low-corridor branch-union predicate, endpoint wall "
                "labels, moving high-grid count, mod-35 clock, endpoint-rank "
                "subcolor, phase-branch pullback, finite transient packets"
            ),
            destroys="non-AP wall geometry and arbitrary owner adjacency",
            proof_pull=(
                "splice HYP-3454/HYP-3456/HYP-3457/HYP-3458/HYP-3459/HYP-3460 "
                "into HYP-3439 after the HYP-3431 fixed-corridor identity is "
                "imported"
            ),
            challenged_assumption=(
                "AP84 still needs sampled transients or an unnamed mod-35 law"
            ),
        ),
        Carrier(
            name="survivor_gate_colored_payload",
            refs=("HYP-3438", "HYP-3436", "HYP-3425"),
            scores={
                "predicate": 5,
                "color_payload": 4,
                "extension_payload": 4,
                "endpoint_owner": 5,
                "finite_audit": 5,
                "frontier_fit": 4,
                "scalar_guardrail": 5,
                "discharge": 4,
            },
            preserves=(
                "mixed E_safe component gate word plus branch mask, endpoint "
                "walls, parent even walls, minimal B0/B1 owner covers, and "
                "legal sidecar route"
            ),
            destroys="raw runner order and scalar survivor mass",
            proof_pull=(
                "promote every survivor gate to a colored extension-orbit row "
                "before quotienting by word type"
            ),
            challenged_assumption="gate word type alone is a safe quotient",
        ),
        Carrier(
            name="observer_cut_orbit_ledger",
            refs=("HYP-3056", "HYP-3054", "LTT-101"),
            scores={
                "predicate": 5,
                "color_payload": 3,
                "extension_payload": 5,
                "endpoint_owner": 4,
                "finite_audit": 3,
                "frontier_fit": 4,
                "scalar_guardrail": 5,
                "discharge": 5,
            },
            preserves=(
                "quotient action, observer address, cut-payload orbit, changed "
                "LRC predicate, separating sidecar, discharge mode, residual "
                "debt"
            ),
            destroys="nothing by itself; it forces every forgetting step to name its loss",
            proof_pull=(
                "make gate extension legal only after reconstruction, exactness, "
                "dual annihilation, descent, boundary stop, or named debt"
            ),
            challenged_assumption="payload retention can be stated as a scalar",
        ),
        Carrier(
            name="colored_resonance_half_boundary",
            refs=("HYP-2595", "HYP-2594", "LTI-048"),
            scores={
                "predicate": 4,
                "color_payload": 5,
                "extension_payload": 3,
                "endpoint_owner": 3,
                "finite_audit": 4,
                "frontier_fit": 4,
                "scalar_guardrail": 5,
                "discharge": 4,
            },
            preserves=(
                "mod-14 color, half-boundary endpoint correction, Fourier "
                "resonance conditions, and O(k+c_GP) discrepancy target"
            ),
            destroys="local gate-gluing topology unless endpoint payload is added",
            proof_pull=(
                "turn branch-mask and endpoint-wall events into a color-resonant "
                "boundary vector instead of charging raw component count"
            ),
            challenged_assumption="micro-component count K measures true deficit",
        ),
        Carrier(
            name="minimal_two_color_bad_core",
            refs=("HYP-3436", "HYP-3425", "HYP-3426"),
            scores={
                "predicate": 5,
                "color_payload": 4,
                "extension_payload": 3,
                "endpoint_owner": 4,
                "finite_audit": 5,
                "frontier_fit": 4,
                "scalar_guardrail": 4,
                "discharge": 4,
            },
            preserves=(
                "E_safe minus B0_odd cap B1_odd, local minimal B0/B1 owner "
                "covers, and branch-union identity"
            ),
            destroys="long-range concatenation data between local bad blocks",
            proof_pull=(
                "use two-color local covers as the base sheaf whose extension "
                "obstructions are measured by survivor-gate deltas"
            ),
            challenged_assumption="local two-color bad blocks can concatenate freely",
        ),
        Carrier(
            name="A000568_edge_envelope_controlled_forgetting",
            refs=("HYP-3134", "HYP-3133", "LTT-159", "LTT-160"),
            scores={
                "predicate": 4,
                "color_payload": 2,
                "extension_payload": 5,
                "endpoint_owner": 4,
                "finite_audit": 4,
                "frontier_fit": 3,
                "scalar_guardrail": 5,
                "discharge": 4,
            },
            preserves=(
                "middle global-consistency quotient between sector shadow and "
                "paired endpoint-deletion child deck"
            ),
            destroys="direct LRC geometry until a packet schema attaches it",
            proof_pull=(
                "treat a raw gate word like an A000568 middle shadow and keep "
                "paired endpoint/owner child data until a legal quotient is named"
            ),
            challenged_assumption="extension counts are proof certificates",
        ),
        Carrier(
            name="PH_bad_coloring_outer_extension_rank",
            refs=("HYP-2247", "HYP-2243"),
            scores={
                "predicate": 3,
                "color_payload": 4,
                "extension_payload": 4,
                "endpoint_owner": 2,
                "finite_audit": 4,
                "frontier_fit": 2,
                "scalar_guardrail": 4,
                "discharge": 3,
            },
            preserves=(
                "coloring plus coherent bad-child extension profile under outer "
                "extension"
            ),
            destroys="metric LRC endpoint geometry",
            proof_pull=(
                "reuse the warning: side-choice colors need extension rank, not "
                "only color labels"
            ),
            challenged_assumption="a color histogram predicts extendability",
        ),
        Carrier(
            name="metagraph_GF2_color_boundary",
            refs=("HYP-2250", "metagraph_blue_black_parity_s675"),
            scores={
                "predicate": 3,
                "color_payload": 5,
                "extension_payload": 3,
                "endpoint_owner": 2,
                "finite_audit": 4,
                "frontier_fit": 2,
                "scalar_guardrail": 4,
                "discharge": 3,
            },
            preserves=(
                "color layer as a GF(2) boundary/1-chain rather than a visible "
                "edge color"
            ),
            destroys="the LRC branch/endpoint semantics until reattached",
            proof_pull=(
                "encode branch-mask imbalance and owner deltas as a boundary "
                "vector whose zero/nonzero status controls gluing"
            ),
            challenged_assumption="visible blue/black colors are the invariant",
        ),
        Carrier(
            name="phase_color_reservoir",
            refs=("HYP-2593", "lrc14_phase_color_reservoir_codex"),
            scores={
                "predicate": 3,
                "color_payload": 5,
                "extension_payload": 2,
                "endpoint_owner": 2,
                "finite_audit": 4,
                "frontier_fit": 3,
                "scalar_guardrail": 4,
                "discharge": 3,
            },
            preserves=(
                "CRT color reservoir C_b(E), positive Sigma floor, and actual "
                "grid-hit identity at q=14V"
            ),
            destroys="component-cover extension topology",
            proof_pull=(
                "supply the color layer under a future finite V<711 placement "
                "check after HYP-2595 is formalized"
            ),
            challenged_assumption="mass floor alone chooses the covering-floor exit",
        ),
        Carrier(
            name="raw_gate_word_shadow",
            refs=("HYP-3438",),
            scores={
                "predicate": 2,
                "color_payload": 2,
                "extension_payload": 2,
                "endpoint_owner": 1,
                "finite_audit": 5,
                "frontier_fit": 3,
                "scalar_guardrail": 2,
                "discharge": 2,
            },
            preserves="the B/S word skeleton of a mixed component",
            destroys="branch mask, owner delta, endpoint wall, and extension status",
            proof_pull="use only as a shadow field inside the richer gate orbit",
            challenged_assumption="the word B-S-B-S-B is the obstruction by itself",
        ),
        Carrier(
            name="raw_component_count",
            refs=("HYP-2594", "HYP-3436"),
            scores={
                "predicate": 1,
                "color_payload": 1,
                "extension_payload": 1,
                "endpoint_owner": 1,
                "finite_audit": 4,
                "frontier_fit": 1,
                "scalar_guardrail": 1,
                "discharge": 1,
            },
            preserves="a size/count proxy",
            destroys="color resonance, boundary cancellation, and extension payload",
            proof_pull="demote to a pessimistic guardrail or finite audit statistic",
            challenged_assumption="large K or many components means proof hardness",
        ),
    ]


def beats(a: Carrier, b: Carrier, tie_order: dict[str, int]) -> bool:
    a_weight = a.weighted_score()
    b_weight = b.weighted_score()
    if a_weight != b_weight:
        return a_weight > b_weight
    if a.score_tuple() != b.score_tuple():
        return a.score_tuple() > b.score_tuple()
    return tie_order[a.name] < tie_order[b.name]


def tournament_edges(carriers: list[Carrier]) -> dict[str, set[str]]:
    tie_order = {carrier.name: idx for idx, carrier in enumerate(carriers)}
    edges = {carrier.name: set() for carrier in carriers}
    for a, b in combinations(carriers, 2):
        if beats(a, b, tie_order):
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    return edges


def count_directed_3cycles(edges: dict[str, set[str]]) -> int:
    total = 0
    names = list(edges)
    for a, b, c in combinations(names, 3):
        ab = b in edges[a]
        bc = c in edges[b]
        ca = a in edges[c]
        ac = c in edges[a]
        cb = b in edges[c]
        ba = a in edges[b]
        if (ab and bc and ca) or (ac and cb and ba):
            total += 1
    return total


def strongly_connected_components(edges: dict[str, set[str]]) -> list[list[str]]:
    names = list(edges)
    reverse = {name: set() for name in names}
    for src, targets in edges.items():
        for target in targets:
            reverse[target].add(src)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in edges[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    components: list[list[str]] = []
    seen.clear()

    def rdfs(v: str, component: list[str]) -> None:
        seen.add(v)
        component.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, component)

    for name in reversed(order):
        if name not in seen:
            component: list[str] = []
            rdfs(name, component)
            components.append(sorted(component))
    return components


def hamiltonian_path(carriers: list[Carrier], edges: dict[str, set[str]]) -> list[str]:
    path: list[str] = []
    for vertex in sorted(carriers, key=lambda c: c.weighted_score(), reverse=True):
        inserted = False
        for idx in range(len(path) + 1):
            candidate = path[:idx] + [vertex.name] + path[idx:]
            if all(candidate[i + 1] in edges[candidate[i]] for i in range(len(candidate) - 1)):
                path = candidate
                inserted = True
                break
        if not inserted:
            raise RuntimeError(f"failed to insert {vertex.name}")
    return path


def format_table(rows: Iterable[tuple[str, str, str]]) -> str:
    widths = [0, 0, 0]
    material = list(rows)
    for row in material:
        for idx, item in enumerate(row):
            widths[idx] = max(widths[idx], len(item))
    return "\n".join(
        f"{row[0]:<{widths[0]}}  {row[1]:<{widths[1]}}  {row[2]}"
        for row in material
    )


def main() -> None:
    carriers = carrier_bank()
    edges = tournament_edges(carriers)
    outdegrees = {name: len(targets) for name, targets in edges.items()}
    score_hist: dict[int, int] = {}
    for score in outdegrees.values():
        score_hist[score] = score_hist.get(score, 0) + 1
    path = hamiltonian_path(carriers, edges)
    sccs = strongly_connected_components(edges)

    print("LRC14 colored-extension gate carrier scout")
    print("source=codex-2026-06-29")
    print("claim_status=EXPLORATORY_SYNTHESIS_NOT_PROOF")
    print()

    print("Tournament Analysis")
    print("vertices=proof carriers and quotient shadows")
    print(
        "pairwise_observable=predicate+color_payload+extension_payload+"
        "endpoint_owner+finite_audit+frontier_fit+scalar_guardrail+discharge"
    )
    print("switch=higher weighted score; tie=score tuple then proof-facing order")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={count_directed_3cycles(edges)}")
    print(f"scc_sizes={sorted((len(component) for component in sccs), reverse=True)}")
    print("hamiltonian_path=")
    for item in path:
        print(f"  {item}")
    print()

    ranked = sorted(carriers, key=lambda c: (-c.weighted_score(), c.name))
    rows = [
        ("rank", "carrier", "weighted_score refs"),
        ("----", "-------", "-------------------"),
    ]
    for idx, carrier in enumerate(ranked, 1):
        rows.append(
            (
                str(idx),
                carrier.name,
                f"{carrier.weighted_score()} {'/'.join(carrier.refs)}",
            )
        )
    print("Carrier Ranking")
    print(format_table(rows))
    print()

    print("Top Carrier Payloads")
    for carrier in ranked[:6]:
        print(f"- {carrier.name}")
        print(f"  preserves: {carrier.preserves}")
        print(f"  destroys: {carrier.destroys}")
        print(f"  proof_pull: {carrier.proof_pull}")
    print()

    print("Proof-Facing Obligations")
    obligations = [
        (
            "O1_color_boundary_vector",
            "Promote branch masks, endpoint-wall labels, and owner deltas to a "
            "colored boundary vector.  HYP-2595 says raw component K is the "
            "wrong charge; surviving color-compatible resonance is the charge.",
        ),
        (
            "O2_gate_extension_orbit",
            "Define a gate-extension orbit modulo row/mirror automorphisms with "
            "fields: gate word, branch mask, endpoint walls, minimal B0/B1 "
            "owner covers, cover-delta vector, and low-rank escape availability.",
        ),
        (
            "O3_random031_finite_clause",
            "Prove the HYP-3455 mirror B-S-B-S-B pair on owner union "
            "(23,45,93,113,147,169,173) cannot glue to a full cover unless a "
            "rank<=2 escape, endpoint spine, owner-current imbalance, two-adic "
            "descent, or SPEC/Rprime debt fires.",
        ),
        (
            "O4_AP84_floor_splice",
            "Use HYP-3456 and HYP-3457 as closed endpoint-color sidecars: the "
            "remaining AP84 task is the HYP-3431 fixed-corridor identity plus "
            "the splice into HYP-3439, not another sampled clock.",
        ),
        (
            "O5_controlled_forgetting",
            "Before forgetting to a raw word, raw A000568 shadow, or raw count, "
            "record the observer-cut orbit and discharge mode from HYP-3056.",
        ),
        (
            "O6_finite_colored_placement",
            "If HYP-2595 is formalized, combine the O(k+c_GP) deficit target "
            "with the HYP-2593 Sigma floor and exact V<711 checking; then feed "
            "the surviving color residue into the gate-extension ledger.",
        ),
    ]
    for name, text in obligations:
        print(f"- {name}: {text}")
    print()

    print("Assumption Challenge")
    print(
        "alternate_vertices_considered=runners,gaps,fixed sections,section "
        "boundaries,wall crossings,residues,cover arcs,endpoint walls,gate "
        "words,owner-delta events,observer cuts,color boundary vectors,A000568 "
        "shadows,AP84 floor packets,proof obligations"
    )
    print(
        "chosen_quotient=colored extension-orbit rows attached to survivor gates "
        "and AP/random031 sidecars"
    )
    print(
        "preserved_predicate=two-color covering-floor gluing predicate plus "
        "existence of a legal LRC14 escape/discharge"
    )
    print(
        "destroyed_information=arbitrary runner order, scalar survivor mass, raw "
        "component count, and untyped extension counts"
    )
    print(
        "challenged_assumption=colors are labels and extensions are counts; here "
        "colors are boundary charges and extensions are gluing orbits"
    )


if __name__ == "__main__":
    main()
