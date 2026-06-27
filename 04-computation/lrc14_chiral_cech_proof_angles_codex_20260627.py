#!/usr/bin/env python3
"""HYP-3123 scout: chiral base-stalk guard and normal-fan/Cech route.

This is a synthesis scout, not an LRC14 proof.

It deliberately chooses two proof angles that are different from the broad
coordinate-resurrection calculus of HYP-3118:

1. A chiral/base-stalk guard: every mirror, converse, edge-sector, or
   rootless perspective quotient must either carry a Z/2 orientation sidecar
   or prove that the collapsed pair has the same terminal packet exit.
2. A normal-fan/Cech finite-ruler route: component control should be expressed
   in normalized slow/ruler coordinates with active normal-fan chambers,
   Cech/barcode data, owner currents, and a finite denominator-net threshold.

Tournament Analysis declaration:
  vertices: proof mechanisms and support carriers, not runners, arcs, roots,
            residues, or raw tournament classes;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: carrier A beats B if it wins more axes; ties use the fixed
                Hamiltonian path below;
  tie Hamiltonian path: chiral guard, normal-fan/Cech, first-obstruction,
                        observer payload, endpoint Phi, Lean frontier,
                        raw component count, raw chirality count.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Dict, Iterable, List, Sequence, Set, Tuple


BASE_STALK = frozenset(
    {
        "finite_address",
        "observer_gluing",
        "endpoint_owner",
        "uniformity",
    }
)

AXES = (
    "predicate_retention",
    "lost_coordinate_control",
    "frontier_adjacency",
    "finite_checkability",
    "cross_domain_transfer",
    "base_stalk_compatibility",
    "formalization_readiness",
    "false_positive_detection",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    angle: str
    anchors: Tuple[str, ...]
    tags: Set[str]
    coords: Set[str]
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_fields: Tuple[str, ...]
    theorem_shape: str


@dataclass(frozen=True)
class Obligation:
    name: str
    anchors: Tuple[str, ...]
    tags: Set[str]
    required: Set[str]
    repaired_by: str
    terminal_test: str


@dataclass(frozen=True)
class LadderStep:
    ladder: str
    level: int
    name: str
    retained: Tuple[str, ...]
    forgotten: str
    next_payload: str


def carrier(
    name: str,
    angle: str,
    anchors: Sequence[str],
    tags: Iterable[str],
    coords: Iterable[str],
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_fields: Sequence[str],
    theorem_shape: str,
) -> Carrier:
    assert len(scores) == len(AXES)
    return Carrier(
        name=name,
        angle=angle,
        anchors=tuple(anchors),
        tags=set(tags),
        coords=set(coords),
        axes=dict(zip(AXES, scores)),
        preserves=preserves,
        destroys=destroys,
        next_fields=tuple(next_fields),
        theorem_shape=theorem_shape,
    )


CARRIERS: List[Carrier] = [
    carrier(
        "chiral_base_stalk_guard",
        "angle_A_chiral_guard",
        ("HYP-3106", "HYP-3120", "S369", "S370", "HYP-3056", "HYP-3118"),
        {
            "chirality",
            "observer",
            "payload",
            "base_stalk",
            "f7",
            "endpoint_owner",
            "state_lift",
            "tournament_perspective",
        },
        BASE_STALK
        | {
            "chiral_guard_word",
            "cross_sector_orientation_word",
            "mirror_pair_id",
            "endpoint_owner_cocycle",
            "state_lift_sign",
            "rootless_fixed_term",
        },
        (5, 5, 5, 4, 5, 5, 4, 5),
        "terminal packet equivalence after mirror/converse identifications",
        "raw class count, raw A000568 equality, and mirror-pair labels unless the Z/2 sidecar is retained",
        (
            "chiral_guard_word",
            "mirror_pair_id",
            "cross_sector_orientation_word",
            "endpoint_owner_cocycle",
            "state_lift_sign",
            "terminal_exit_or_named_debt",
        ),
        "If a quotient identifies a converse/mirror pair, either both fibers have the same finite-address or observer-gluing exit, or the chiral guard word is nonzero and routes to endpoint-owner, first-obstruction, or F7/THM-572 debt.",
    ),
    carrier(
        "normal_fan_cech_finite_ruler",
        "angle_B_topological_ruler",
        ("HYP-3101", "HYP-3025", "HYP-3018", "HYP-3096", "THM-565", "THM-573"),
        {
            "normal_fan",
            "cech",
            "barcode",
            "finite_ruler",
            "component_bound",
            "owner_current",
            "topology",
            "lean",
        },
        BASE_STALK
        | {
            "component_bound",
            "normal_fan_chamber_id",
            "closed_arc_cech_beta",
            "open_tope_component_count",
            "barcode_persistence_word",
            "finite_ruler_denominator_threshold",
            "owner_current_word",
            "boundary_cocircuit_word",
        },
        (5, 4, 5, 5, 5, 4, 4, 4),
        "normalized component lower/upper control for the witness route",
        "exact scale, branch route, and endpoint owner if topological summaries are isolated",
        (
            "normal_fan_chamber_id",
            "closed_arc_cech_beta",
            "barcode_persistence_word",
            "finite_ruler_denominator_threshold",
            "owner_current_word",
            "component_bound_status",
            "terminal_exit_or_named_debt",
        ),
        "Every non-tight THM-573 residual row has a normalized finite-ruler component packet, is an AP/GW closed-boundary H1 atom, or emits a named F7/THM-572 good-cover quotient defect.",
    ),
    carrier(
        "first_obstruction_syndrome_bridge",
        "bridge_syndrome",
        ("HYP-3102", "HYP-3071", "HYP-3056", "HYP-3120", "THM-572"),
        {
            "obstruction",
            "cocycle",
            "syndrome",
            "observer",
            "state_lift",
            "f7",
            "coding_theory",
        },
        {
            "observer_gluing",
            "endpoint_owner",
            "cocycle_exactness",
            "first_obstruction_basis_vector",
            "certificate_cycle_image_status",
            "dual_annihilator_status",
            "F7_THM572_state_lift_status",
            "terminal_exit",
        },
        (5, 5, 4, 4, 5, 4, 4, 5),
        "the first hidden payload difference is exact, generated, killed, descended, or named as debt",
        "raw route identity after the payload is encoded as a syndrome",
        (
            "first_obstruction_basis_vector",
            "certificate_cycle_image_status",
            "dual_annihilator_status",
            "F7_THM572_state_lift_status",
            "syndrome_decoder_status",
            "terminal_exit_or_named_debt",
        ),
        "The chiral and Cech routes meet at a syndrome test: any leftover difference must be zero, generated by known cycles, dual-annihilated, or routed to state-lift debt.",
    ),
    carrier(
        "observer_cut_payload_orbit",
        "support_payload",
        ("HYP-3054", "HYP-3056", "HYP-3106", "HYP-3120"),
        {
            "observer",
            "payload",
            "orbit",
            "gluing",
            "perspective",
            "extension",
        },
        {
            "observer_gluing",
            "endpoint_owner",
            "observer_cut_payload_orbit",
            "cross_sector_orientation_word",
            "incident_sector_deck",
            "terminal_exit",
        },
        (5, 5, 5, 4, 4, 4, 4, 5),
        "the first forgotten observer-extension payload over a visible fiber",
        "node-only perspective type and raw class equality",
        (
            "observer_cut_payload_orbit",
            "incident_sector_deck",
            "cross_sector_orientation_word",
            "payload_orbit_terminal_status",
        ),
        "A quotient may forget root labels only after the observer-cut payload orbit and cross-sector word are attached.",
    ),
    carrier(
        "endpoint_phi_packet_receiver",
        "support_endpoint",
        ("HYP-2108", "HYP-2112", "HYP-3116", "HYP-3119"),
        {
            "endpoint",
            "phi",
            "activation",
            "finite_address",
            "proof_circuit",
            "lean",
        },
        {
            "finite_address",
            "endpoint_owner",
            "uniformity",
            "Phi_gap",
            "P_sign",
            "endpoint_cover_activation_vector",
            "proof_circuit_missing_input_vector",
        },
        (5, 4, 5, 5, 3, 5, 5, 4),
        "endpoint activation as the exact proof-circuit output wire",
        "observer payload and topological component data unless attached as sidecars",
        (
            "endpoint_cover_activation_vector",
            "Phi_gap",
            "P_sign",
            "proof_circuit_missing_input_vector",
            "endpoint_phi_terminal_status",
        ),
        "The exact gap/activation circuit is the receiver for chiral and Cech sidecars, not a replacement for them.",
    ),
    carrier(
        "lean_frontier_packet_bus",
        "support_formal",
        ("HYP-3107", "HYP-3095", "HYP-3096", "HYP-3098"),
        {
            "lean",
            "finite_address",
            "observer_gluing",
            "proof_frontier",
            "formal",
            "terminal_exit",
        },
        BASE_STALK
        | {
            "terminal_exit",
            "ObserverGluingCertificate",
            "FiniteAddressBranchPacket",
            "LRC14Statement",
        },
        (5, 4, 5, 5, 3, 5, 5, 5),
        "the formal packet interface that can consume either proof angle",
        "domain-specific sidecar payload if the Lean bus is used without producers",
        (
            "FiniteAddressBranchPacket",
            "ObserverGluingCertificate",
            "LRC14Statement_handoff",
            "terminal_exit_or_named_debt",
        ),
        "Both angles should emit rows into the existing Lean-facing packet bus rather than inventing a parallel theorem language.",
    ),
    carrier(
        "raw_component_count_shadow",
        "negative_control",
        ("HYP-3096", "HYP-3101"),
        {"component_count", "raw_time", "negative_control"},
        {"component_count", "raw_time_largest_arc"},
        (3, 2, 4, 5, 3, 2, 3, 2),
        "a useful alarm that a direct chart is complicated",
        "normalization, owner currents, Cech beta, and finite-ruler threshold",
        ("component_count",),
        "Raw component count is telemetry until normalized slow/ruler coordinates and owner/Cech sidecars are present.",
    ),
    carrier(
        "raw_chirality_count_shadow",
        "negative_control",
        ("S369", "S370", "HYP-3047", "HYP-3106"),
        {"chirality", "raw_count", "negative_control", "tournament"},
        {"self_converse_count", "chiral_count", "A000568_count"},
        (2, 1, 2, 5, 3, 1, 2, 1),
        "a count-level warning that mirror orientation was collapsed",
        "edge-sector coupling, endpoint owners, and terminal packet identity",
        ("self_converse_count", "chiral_count"),
        "The 12/44/56 numerology is a guardrail until cross-sector orientation and endpoint-owner payloads are attached.",
    ),
]


OBLIGATIONS: List[Obligation] = [
    Obligation(
        "HYP-3098_observer_gluing_rows",
        ("HYP-3098", "HYP-3095", "HYP-3096"),
        {"observer", "payload", "finite_address", "lean"},
        {
            "finite_address",
            "observer_gluing",
            "endpoint_owner",
            "uniformity",
            "observer_cut_payload_orbit",
            "terminal_exit",
        },
        "chiral guard supplies the missing observer payload; Cech route supplies normalized component packet",
        "ObserverGluingCertificate or finite-address exit",
    ),
    Obligation(
        "HYP-3101_component_bound",
        ("HYP-3101", "HYP-3025", "HYP-3018", "THM-565"),
        {"normal_fan", "cech", "component_bound", "finite_ruler", "topology"},
        {
            "finite_address",
            "component_bound",
            "normal_fan_chamber_id",
            "closed_arc_cech_beta",
            "finite_ruler_denominator_threshold",
            "owner_current_word",
            "terminal_exit",
        },
        "normal-fan/Cech route is the direct producer",
        "bounded components plus positive measure gives finite-ruler hit, or named AP/GW/F7 debt",
    ),
    Obligation(
        "HYP-3102_first_obstruction",
        ("HYP-3102", "HYP-3071", "THM-572"),
        {"obstruction", "cocycle", "observer", "f7", "state_lift"},
        {
            "observer_gluing",
            "endpoint_owner",
            "cocycle_exactness",
            "first_obstruction_basis_vector",
            "certificate_cycle_image_status",
            "dual_annihilator_status",
            "F7_THM572_state_lift_status",
        },
        "syndrome bridge tests residue from either main angle",
        "zero/exact/generated/killed/descended or F7/THM-572 debt",
    ),
    Obligation(
        "HYP-3112_lee_yang_ear_payload",
        ("HYP-3112", "HYP-3109", "HYP-3108"),
        {"lee_yang", "ear", "observer", "chirality", "payload"},
        {
            "root_ear_payload",
            "observer_gluing",
            "endpoint_owner",
            "chiral_guard_word",
            "cross_sector_orientation_word",
            "terminal_exit",
        },
        "chiral guard names the ear parity and mirror orientation that root motion forgets",
        "root motion reconstructed or routed to nonnested-ear / first-obstruction debt",
    ),
    Obligation(
        "HYP-3115_ising_domain_walls",
        ("HYP-3115", "HYP-3116", "HYP-3119"),
        {"ising", "relation_lattice", "domain_wall", "chirality", "proof_circuit"},
        {
            "relation_lattice_shape",
            "root_ear_payload",
            "chiral_guard_word",
            "endpoint_cover_activation_vector",
            "proof_circuit_missing_input_vector",
            "terminal_exit",
        },
        "chiral guard distinguishes legal wall orientation; endpoint Phi receives the gate stack",
        "domain wall becomes proof-circuit gate or named root/ear/chiral debt",
    ),
    Obligation(
        "THM-573_level7_residual",
        ("THM-573", "HYP-3083", "HYP-3096", "HYP-3101"),
        {"level7", "finite_ruler", "state_lift", "component_bound", "f7"},
        {
            "finite_address",
            "endpoint_owner",
            "component_bound",
            "normal_fan_chamber_id",
            "finite_ruler_denominator_threshold",
            "F7_THM572_state_lift_status",
        },
        "Cech route supplies component/finite-ruler packet; chiral guard stops anonymous F7",
        "finite grid witness, AP/GW atom, or named F7/THM-572 state-lift",
    ),
    Obligation(
        "Lean_PartA_finite_ruler",
        ("LRCWitnessPartA.lean", "LRCGoodSet.lean", "LRCGapReach.lean"),
        {"lean", "finite_ruler", "good_set", "approximation"},
        {
            "finite_ruler_denominator_threshold",
            "component_bound",
            "goodSet_margin",
            "arcCount",
            "rhoK_approximation",
            "terminal_exit",
        },
        "normal-fan/Cech route should emit the missing concrete arc bound and finite-ruler threshold",
        "Part-A finite approximation inequality for actual events",
    ),
]


LADDERS = [
    LadderStep(
        "chiral_guard",
        0,
        "raw count coincidence",
        ("A000568_count", "self_converse_count", "chiral_count"),
        "edge sectors, endpoint owners, and cross-sector orientation",
        "ordered-pair or directed-edge sector deck",
    ),
    LadderStep(
        "chiral_guard",
        1,
        "edge-sector deck",
        ("incident_sector_deck", "sector_size_deck", "sector_internal_deck"),
        "converse/mirror orientation between sectors",
        "cross_sector_orientation_word",
    ),
    LadderStep(
        "chiral_guard",
        2,
        "Z/2 orientation local system",
        ("cross_sector_orientation_word", "mirror_pair_id", "chiral_guard_word"),
        "endpoint-owner and state-lift consequences",
        "endpoint_owner_cocycle + state_lift_sign",
    ),
    LadderStep(
        "chiral_guard",
        3,
        "base-stalk terminal guard",
        ("finite_address", "observer_gluing", "endpoint_owner", "uniformity", "terminal_exit"),
        "nothing terminal; failure becomes named debt",
        "finite-address exit, observer-gluing exit, or F7/THM-572 debt",
    ),
    LadderStep(
        "normal_fan_cech",
        0,
        "raw direct component count",
        ("component_count", "largest_direct_arc", "raw_time_largest_arc"),
        "normalized slow/ruler scale and active owner geometry",
        "normalized chamber packet",
    ),
    LadderStep(
        "normal_fan_cech",
        1,
        "active normal-fan chamber",
        ("normal_fan_chamber_id", "active_support_word", "owner_current_word"),
        "closed/open topology and barcode persistence",
        "Cech/barcode packet",
    ),
    LadderStep(
        "normal_fan_cech",
        2,
        "Cech/barcode component packet",
        ("closed_arc_cech_beta", "open_tope_component_count", "barcode_persistence_word"),
        "finite denominator-net threshold",
        "finite_ruler_denominator_threshold",
    ),
    LadderStep(
        "normal_fan_cech",
        3,
        "finite-ruler terminal packet",
        ("component_bound", "finite_ruler_denominator_threshold", "terminal_exit"),
        "nothing terminal; failure becomes named AP/GW/F7 debt",
        "finite grid witness, AP/GW H1 atom, or F7/THM-572 debt",
    ),
]


TIE_PATH = [
    "chiral_base_stalk_guard",
    "normal_fan_cech_finite_ruler",
    "first_obstruction_syndrome_bridge",
    "observer_cut_payload_orbit",
    "endpoint_phi_packet_receiver",
    "lean_frontier_packet_bus",
    "raw_component_count_shadow",
    "raw_chirality_count_shadow",
]


def base_stalk_status(c: Carrier) -> Tuple[Tuple[str, ...], bool]:
    missing = tuple(sorted(BASE_STALK - c.coords))
    return missing, not missing


def bridge_score(c: Carrier, o: Obligation) -> int:
    coord_overlap = len(c.coords & o.required)
    tag_overlap = len(c.tags & o.tags)
    missing_base = len(BASE_STALK - c.coords)
    return sum(c.axes.values()) + 3 * coord_overlap + tag_overlap - 2 * missing_base


def bridge_matrix() -> List[Tuple[int, str, str, Tuple[str, ...], Tuple[str, ...]]]:
    rows = []
    for c in CARRIERS:
        for o in OBLIGATIONS:
            score = bridge_score(c, o)
            supplied = tuple(sorted(c.coords & o.required))
            missing = tuple(sorted(o.required - c.coords))
            rows.append((score, c.name, o.name, supplied, missing))
    return sorted(rows, reverse=True)


def preferred(a: Carrier, b: Carrier) -> Carrier:
    wins_a = sum(1 for axis in AXES if a.axes[axis] > b.axes[axis])
    wins_b = sum(1 for axis in AXES if b.axes[axis] > a.axes[axis])
    if wins_a != wins_b:
        return a if wins_a > wins_b else b
    return a if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b


def tournament_edges() -> Dict[str, Set[str]]:
    edges: Dict[str, Set[str]] = {c.name: set() for c in CARRIERS}
    by_name = {c.name: c for c in CARRIERS}
    for left, right in combinations(TIE_PATH, 2):
        winner = preferred(by_name[left], by_name[right]).name
        loser = right if winner == left else left
        edges[winner].add(loser)
    return edges


def directed_3cycles(edges: Dict[str, Set[str]]) -> List[Tuple[str, str, str]]:
    cycles = []
    names = list(edges)
    for a, b, c in combinations(names, 3):
        triples = ((a, b, c), (a, c, b))
        for x, y, z in triples:
            if y in edges[x] and z in edges[y] and x in edges[z]:
                cycles.append((x, y, z))
                break
    return cycles


def sccs(edges: Dict[str, Set[str]]) -> List[List[str]]:
    names = list(edges)
    reverse = {n: set() for n in names}
    for u, outs in edges.items():
        for v in outs:
            reverse[v].add(u)

    seen: Set[str] = set()
    order: List[str] = []

    def dfs(u: str) -> None:
        seen.add(u)
        for v in edges[u]:
            if v not in seen:
                dfs(v)
        order.append(u)

    for n in names:
        if n not in seen:
            dfs(n)

    comps: List[List[str]] = []
    seen.clear()

    def rdfs(u: str, comp: List[str]) -> None:
        seen.add(u)
        comp.append(u)
        for v in reverse[u]:
            if v not in seen:
                rdfs(v, comp)

    for n in reversed(order):
        if n not in seen:
            comp: List[str] = []
            rdfs(n, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_paths(edges: Dict[str, Set[str]]) -> List[Tuple[str, ...]]:
    names = list(edges)
    paths = []
    for perm in permutations(names):
        if all(perm[i + 1] in edges[perm[i]] for i in range(len(perm) - 1)):
            paths.append(perm)
    return paths


def print_header(title: str) -> None:
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def main() -> None:
    print("HYP-3123 chiral/Cech proof-angle scout -- codex 2026-06-27")
    print("This is route synthesis and payload bookkeeping, not a proof.")

    print_header("ANGLE SUMMARY")
    for c in CARRIERS:
        if c.angle.startswith("angle_"):
            missing, complete = base_stalk_status(c)
            print(f"- {c.name}")
            print(f"  anchors={', '.join(c.anchors)}")
            print(f"  base_stalk_complete={complete}; missing={','.join(missing) if missing else 'none'}")
            print(f"  preserves={c.preserves}")
            print(f"  destroys_or_hides={c.destroys}")
            print(f"  theorem_shape={c.theorem_shape}")
            print(f"  next_fields={', '.join(c.next_fields)}")

    print_header("CONTROLLED-FORGETTING LADDERS")
    for ladder in ("chiral_guard", "normal_fan_cech"):
        print(f"{ladder}:")
        for step in [s for s in LADDERS if s.ladder == ladder]:
            print(f"  L{step.level}: {step.name}")
            print(f"      retained={', '.join(step.retained)}")
            print(f"      forgotten={step.forgotten}")
            print(f"      next_payload={step.next_payload}")

    print_header("OBLIGATION BRIDGE MATRIX: TOP MATCHES")
    for score, carrier_name, obligation_name, supplied, missing in bridge_matrix()[:18]:
        print(f"- score={score:2d} {carrier_name} -> {obligation_name}")
        print(f"  supplied={', '.join(supplied) if supplied else 'none'}")
        print(f"  missing={', '.join(missing) if missing else 'none'}")

    print_header("FIELD PROMOTION CANDIDATES")
    field_counter: Counter[str] = Counter()
    for c in CARRIERS:
        for field in c.next_fields:
            field_counter[field] += 1
    for field, count in field_counter.most_common():
        print(f"{field}: {count}")

    print_header("TOURNAMENT ANALYSIS")
    edges = tournament_edges()
    scores = {n: len(outs) for n, outs in edges.items()}
    score_hist = Counter(scores.values())
    cycles = directed_3cycles(edges)
    comps = sccs(edges)
    hpaths = hamiltonian_paths(edges)
    print(f"vertices={len(CARRIERS)}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={len(cycles)}")
    if cycles:
        for cyc in cycles:
            print(f"  cycle={' -> '.join(cyc)} -> {cyc[0]}")
    print(f"scc_sizes={[len(comp) for comp in comps]}")
    print(f"hamiltonian_path_count={len(hpaths)}")
    if hpaths:
        best = max(hpaths, key=lambda path: sum(scores[name] for name in path))
        print("selected_priority_path=")
        for name in best:
            print(f"  {name} (outscore={scores[name]})")

    print_header("CROSS-DISCIPLINARY CONNECTIONS")
    connections = [
        (
            "Z/2 orientation local systems",
            "chiral guard word; mirror/converse pairs become a local coefficient system over observer fibers",
        ),
        (
            "oriented matroid topes and cocircuits",
            "normal-fan chambers, boundary cocircuit words, and endpoint-owner currents are the finite sign data of the route",
        ),
        (
            "persistent homology / nerve theorem",
            "closed arc-Cech beta and barcode persistence are legal only with owners and finite-ruler thresholds",
        ),
        (
            "syndrome decoding",
            "first obstruction is a finite syndrome: zero, generated, dual-annihilated, descended, or named as F7 debt",
        ),
        (
            "spin/domain-wall physics",
            "Ising walls need chirality and endpoint Phi gates before wall tension becomes proof currency",
        ),
        (
            "observability and Schur complements",
            "rootless perspective collisions should be tested by the first hidden observer-cut payload, not by deeper node memory",
        ),
    ]
    for domain, translation in connections:
        print(f"- {domain}: {translation}")

    print_header("ASSUMPTION CHALLENGE")
    print("Alternate tournament vertex sets considered:")
    print(
        "  runners, gaps, fixed circle sections, section boundaries, wall-crossing events, residues,"
    )
    print(
        "  cover arcs, Fourier modes, matroid circuits, normal-fan chambers, Cech cycles,"
    )
    print("  directed-edge sectors, mirror pairs, and proof obligations.")
    print("Chosen vertices: proof mechanisms and support carriers.")
    print(
        "Preserved predicate: progress toward a finite LRC14 certificate via direct witness, AP/GW boundary atom, finite-address packet, observer-gluing certificate, or named F7/THM-572 debt."
    )
    print(
        "Destroyed information: object-level runner labels and raw scalar counts; they are accepted only after the chiral guard or Cech finite-ruler payload repairs the first forgotten coordinate."
    )
    print(
        "Challenged assumption: the remaining proof route is not another extremal scalar. It is probably a base-stalk maintenance problem with a Z/2 chiral guard and a normalized topological finite-ruler packet."
    )


if __name__ == "__main__":
    main()
