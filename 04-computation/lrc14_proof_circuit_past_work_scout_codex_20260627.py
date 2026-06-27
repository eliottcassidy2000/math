#!/usr/bin/env python3
"""HYP-3117 scout: past-work proof-circuit recompilation for LRC14.

This is a synthesis scout, not a proof.  It searches the existing research
surface for circuit-complexity-shaped proof carriers and recompiles them as a
uniform proof-state circuit for the LRC14 frontier.

Tournament Analysis uses proof gates and sidecars as vertices, not runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


FEATURE_WEIGHTS = {
    "lrc_predicate": 5,
    "uniform_circuit": 5,
    "exact_measure": 4,
    "sidecar_retention": 4,
    "terminal_exit": 5,
    "lean_facing": 3,
    "route_purity": 3,
    "residual_separation": 3,
    "missing_input_detector": 3,
    "counterexample_filter": 2,
    "raw_scalar_penalty": -7,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    source: str
    role: str
    preserves: str
    destroys: str
    missing: str
    features: dict[str, int]

    def score(self) -> int:
        return sum(FEATURE_WEIGHTS[key] * value for key, value in self.features.items())


CARRIERS = [
    Carrier(
        "lean_frontier_obligation_bus",
        "HYP-3107",
        "Ten named proof obligations and the terminal LRC14 interface.",
        "LRC14Statement routing through finite-address or observer-gluing certificates.",
        "Raw runner labels and any scalar shortcut that bypasses the proof fields.",
        "Replace Prop fields by theorem-producing packets.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "terminal_exit": 1,
            "lean_facing": 1,
            "route_purity": 1,
            "missing_input_detector": 1,
            "sidecar_retention": 1,
        },
    ),
    Carrier(
        "first_obstruction_cocycle_gate",
        "HYP-3102, HYP-2997",
        "Turns illegal forgetting into a first obstruction cochain.",
        "Observer/cut payload legality and generated certificate-cycle exits.",
        "Packet identity after the obstruction class has been generated.",
        "Actual residual cochains and a stability proof for the basis.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "sidecar_retention": 1,
            "terminal_exit": 1,
            "route_purity": 1,
            "residual_separation": 1,
            "missing_input_detector": 1,
            "lean_facing": 1,
        },
    ),
    Carrier(
        "labelled_packet_decision_tree",
        "HYP-2961, HYP-2962, HYP-2963",
        "Deterministic counterexample-family decision tree and packet emitter.",
        "Strict-counterexample predicate, qdiv, open-front, packet family, and state-lift flags.",
        "Fine row identity once converted to a fixed-margin labelled packet.",
        "Uniform theorem that every primitive residual emits a legal packet.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "residual_separation": 1,
            "counterexample_filter": 1,
            "missing_input_detector": 1,
        },
    ),
    Carrier(
        "phi_gap_output_wire",
        "HYP-2112",
        "Exact circuit-to-gap functional G(v)=Phi(C).",
        "Positive gap measure and ker Phi as the tight/worry set.",
        "Endpoint-owner history after the exact gap value is formed.",
        "Proof that Phi>0 on every residual multiple-of-n configuration.",
        {
            "lrc_predicate": 1,
            "exact_measure": 1,
            "sidecar_retention": 1,
            "terminal_exit": 1,
            "route_purity": 1,
            "residual_separation": 1,
        },
    ),
    Carrier(
        "endpoint_cover_circuit_gate",
        "HYP-2108",
        "Exact endpoint-cover positivity circuit P(S).",
        "Simultaneous resonance obstruction for all component midpoints.",
        "Which local owner argument proved a positive P term.",
        "A global simultaneous-resonance infeasibility theorem for the residual.",
        {
            "lrc_predicate": 1,
            "exact_measure": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "residual_separation": 1,
            "counterexample_filter": 1,
        },
    ),
    Carrier(
        "fold_gate_virtual_sum",
        "HYP-2114, HYP-2115",
        "3-term fold gate plus hidden virtual-sum nodes.",
        "The delta-clock/fold boundary and why 4-term energy is only a shadow.",
        "Translation-invariant additive energy when hidden sums are collapsed.",
        "Fold pressure theorem with virtual-sum nodes retained.",
        {
            "lrc_predicate": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "residual_separation": 1,
            "missing_input_detector": 1,
        },
    ),
    Carrier(
        "lee_yang_ear_payload_gate",
        "HYP-3109, HYP-3112",
        "Whole PGF root curve plus exact one-runner ear payload.",
        "Root stratum, root-motion reconstruction, and ear payload A_t(E,a).",
        "Raw scalar p0 and raw row identity after the root packet is formed.",
        "Classify domain-wall/root-collision ears as legal exits or forbidden walls.",
        {
            "lrc_predicate": 1,
            "exact_measure": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "residual_separation": 1,
        },
    ),
    Carrier(
        "minkowski_relation_wall_gate",
        "HYP-3062, HYP-3111, HYP-3115",
        "Relation-lattice short-wall and convex-body sidecar.",
        "Relation lattice, covolume/body id, low-height wall deletion, and root-stratum compatibility.",
        "Generic geometry slogans and raw Bravais concentration.",
        "LRC-native convex body whose inequalities are packet predicates.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "missing_input_detector": 1,
        },
    ),
    Carrier(
        "automaton_branch_program_guard",
        "HYP-3008..HYP-3016",
        "Finite automaton or branching-program quotient with magnitude cocycle.",
        "Residue-language state only when magnitude/Farey data separates mixed fibers.",
        "Exact magnitude inside a mixed residue-language fiber if no cocycle is retained.",
        "Fiber-purity theorem or magnitude-cocycle discharge.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "sidecar_retention": 1,
            "route_purity": 1,
            "missing_input_detector": 1,
            "counterexample_filter": 1,
        },
    ),
    Carrier(
        "mciq_uniformity_guard",
        "HYP-3111, HYP-3115",
        "Circuit-complexity audit: fixed inputs, fan-in, depth, and uniformity.",
        "Proof-circuit basis and missing-input vector for any shortcut.",
        "Finite-bank fitted thresholds after their sidecars are compiled.",
        "Train/test and residual-bank circuit with fixed legal inputs.",
        {
            "lrc_predicate": 1,
            "uniform_circuit": 1,
            "sidecar_retention": 1,
            "missing_input_detector": 1,
            "counterexample_filter": 1,
        },
    ),
    Carrier(
        "raw_scalar_maximizer_rule",
        "HYP-3115 warning",
        "One-literal finite-bank classifier such as apex7_error<=5.",
        "Only a bounded-bank signal that suggests a sidecar to test.",
        "Uniformity, basis declaration, root/observer payloads, and proof exits.",
        "All theorem-grade circuit fields.",
        {
            "exact_measure": 1,
            "raw_scalar_penalty": 1,
        },
    ),
]


TIE_PATH = [
    "lean_frontier_obligation_bus",
    "first_obstruction_cocycle_gate",
    "labelled_packet_decision_tree",
    "phi_gap_output_wire",
    "endpoint_cover_circuit_gate",
    "fold_gate_virtual_sum",
    "lee_yang_ear_payload_gate",
    "minkowski_relation_wall_gate",
    "automaton_branch_program_guard",
    "mciq_uniformity_guard",
    "raw_scalar_maximizer_rule",
]


NEW_SIGNALS = [
    "proof_circuit_input_basis_id",
    "proof_circuit_missing_input_vector",
    "uniformity_guard_status",
    "endpoint_cover_P_gate",
    "Phi_gap_output_wire",
    "fold_gate_depth",
    "hidden_virtual_sum_count",
    "automaton_fiber_mixing_bit",
    "magnitude_cocycle_height",
    "first_obstruction_class",
    "lee_yang_ear_payload_mean_level",
    "root_motion_reconstruction_status",
    "relation_wall_class",
    "sidecar_fanin_profile",
    "minimal_certificate_depth",
    "gate_route_purity",
    "terminal_exit_kind",
]


def orient(a: Carrier, b: Carrier) -> tuple[Carrier, Carrier]:
    if a.score() > b.score():
        return a, b
    if b.score() > a.score():
        return b, a
    rank = {name: i for i, name in enumerate(TIE_PATH)}
    return (a, b) if rank[a.name] < rank[b.name] else (b, a)


def build_edges() -> dict[str, set[str]]:
    edges = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        winner, loser = orient(a, b)
        edges[winner.name].add(loser.name)
    return edges


def directed_3cycles(edges: dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(edges, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            count += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            count += 1
    return count


def sccs(edges: dict[str, set[str]]) -> list[list[str]]:
    names = list(edges)
    reverse = {name: set() for name in names}
    for src, dsts in edges.items():
        for dst in dsts:
            reverse[dst].add(src)

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

    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, comp)

    for name in reversed(order):
        if name not in seen:
            comp: list[str] = []
            rdfs(name, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(edges: dict[str, set[str]]) -> tuple[int, list[str]]:
    names = list(edges)
    incoming_count = {name: 0 for name in names}
    for src in names:
        for dst in edges[src]:
            incoming_count[dst] += 1

    count = 0
    first_path: list[str] = []

    def backtrack(path: list[str], used: set[str]) -> None:
        nonlocal count, first_path
        if len(path) == len(names):
            count += 1
            if not first_path:
                first_path = path[:]
            return
        for candidate in names:
            if candidate in used:
                continue
            if path and candidate not in edges[path[-1]]:
                continue
            path.append(candidate)
            used.add(candidate)
            backtrack(path, used)
            used.remove(candidate)
            path.pop()

    for start in sorted(names, key=lambda name: (-incoming_count[name], name)):
        backtrack([start], {start})
    return count, first_path


def print_gate_map() -> None:
    print("MAP 1: past-work gates recompiled as proof-circuit carriers")
    print("-" * 72)
    for carrier in sorted(CARRIERS, key=lambda c: (-c.score(), c.name)):
        print(f"{carrier.name} [{carrier.source}] score={carrier.score()}")
        print(f"  role: {carrier.role}")
        print(f"  preserves: {carrier.preserves}")
        print(f"  destroys: {carrier.destroys}")
        print(f"  missing: {carrier.missing}")
    print()


def print_compiler_map() -> None:
    print("MAP 2: residual-row proof-circuit compiler sketch")
    print("-" * 72)
    steps = [
        ("primitive_residual_row", "emit labelled packet with qdiv, Haar front, owner skeleton, K33/C27, source family"),
        ("labelled_packet_decision_tree", "route to a safe family, live family, or named sporadic/state-lift debt"),
        ("endpoint_cover_circuit_gate", "test simultaneous endpoint resonance by P(S)"),
        ("phi_gap_output_wire", "compute exact positive gap Phi(C), or identify ker Phi worry set"),
        ("fold_gate_virtual_sum", "retain 3-term folds and hidden virtual-sum nodes"),
        ("automaton_branch_program_guard", "reject mixed finite-state fibers unless magnitude cocycle is retained"),
        ("lee_yang_ear_payload_gate", "attach whole PGF root structure and one-runner ear payload"),
        ("minkowski_relation_wall_gate", "declare relation body, short walls, and low-height wall deletion"),
        ("first_obstruction_cocycle_gate", "turn any illegal quotient into a generated obstruction class"),
        ("lean_frontier_obligation_bus", "supply finite-address or observer-gluing terminal certificate"),
    ]
    for left, right in steps:
        print(f"{left} -> {right}")
    print()


def print_tournament() -> None:
    edges = build_edges()
    scores = {carrier.name: carrier.score() for carrier in CARRIERS}
    score_hist = Counter(scores.values())
    path_count, first_path = hamiltonian_path_count(edges)

    print("TOURNAMENT ANALYSIS")
    print("-" * 72)
    print("vertices = proof gates/sidecars, not runners or arcs")
    print(
        "pairwise observable = which carrier preserves the LRC predicate, "
        "retains a legal sidecar, supplies a uniform proof-circuit input, "
        "or exposes a missing-input vector before scalarization"
    )
    print("switch/gauge = larger weighted proof-readiness score; ties follow the declared compiler path")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycles(edges)}")
    print(f"scc_sizes={sorted(len(comp) for comp in sccs(edges))}")
    print(f"hamiltonian_path_count={path_count}")
    print("priority_path=" + " -> ".join(first_path))
    print()


def print_signals() -> None:
    print("NEW SIGNALS TO MEASURE")
    print("-" * 72)
    for signal in NEW_SIGNALS:
        print(signal)
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("-" * 72)
    considered = [
        "runners",
        "gaps",
        "fixed circle sections",
        "section boundaries",
        "wall-crossing events",
        "residues",
        "cover arcs",
        "Fourier modes",
        "matroid circuits",
        "proof obligations",
        "Boolean gates",
        "automaton states",
        "cocycle basis atoms",
        "Lee-Yang root-motion ears",
    ]
    print("candidate vertices considered=" + ", ".join(considered))
    print("chosen vertices=proof gates and sidecar obligations")
    print(
        "preserved predicate=residual packet progress to LRC14Statement through "
        "positive Phi/P, finite-address packets, observer-gluing certificates, "
        "or named state-lift debt"
    )
    print(
        "destroyed data=raw row identity, runner order, raw time, some owner "
        "history, and scalar maximizer labels unless a sidecar retains them"
    )
    print()


def main() -> None:
    print("HYP-3117 LRC14 proof-circuit past-work scout")
    print("=" * 72)
    print_gate_map()
    print_compiler_map()
    print_tournament()
    print_signals()
    print_assumption_challenge()
    print("SUMMARY")
    print("-" * 72)
    print(
        "Circuit complexity is best used here as a uniform proof compiler: "
        "fixed inputs, declared fan-in/depth, missing-input vectors, and "
        "terminal proof exits.  The old endpoint/Phi/additive-circuit work "
        "supplies real gates, while the newer packet, automaton, cocycle, "
        "Lee-Yang, and Lean frontiers supply sidecars and legality checks."
    )


if __name__ == "__main__":
    main()
