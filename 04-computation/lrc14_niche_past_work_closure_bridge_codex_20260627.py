#!/usr/bin/env python3
"""Niche past-work closure bridge for the LRC14 proof frontier.

This script does not try to reprove LRC14.  It turns a repo scour into a
reproducible routing table: older "niche" proof devices are scored against the
current HYP-3098..HYP-3118 frontier obligations, then compared by Tournament
Analysis.

Tournament Analysis declaration:
  vertices: past-work proof carriers, not runners, speeds, roots, or residues;
  pairwise observable: retained proof-packet payload across eight axes;
  switch/gauge: majority of axis scores, with a fixed tie path;
  tie Hamiltonian path: finite-address carriers before raw analogy carriers.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Sequence, Set, Tuple


AXES = (
    "predicate_retention",
    "frontier_adjacency",
    "finite_checkability",
    "forgotten_coordinate_control",
    "loaded_denominator_control",
    "formalization_readiness",
    "scalar_false_positive_detection",
    "extension_payload_retention",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    anchors: Tuple[str, ...]
    tags: Set[str]
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_use: str


@dataclass(frozen=True)
class Obligation:
    name: str
    tags: Set[str]
    weights: Dict[str, int]
    theorem_shape: str


def carrier(
    name: str,
    anchors: Sequence[str],
    tags: Iterable[str],
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_use: str,
) -> Carrier:
    assert len(scores) == len(AXES)
    return Carrier(
        name=name,
        anchors=tuple(anchors),
        tags=set(tags),
        axes=dict(zip(AXES, scores)),
        preserves=preserves,
        destroys=destroys,
        next_use=next_use,
    )


CARRIERS: List[Carrier] = [
    carrier(
        "finite_address_phi_tuple",
        ("HYP-3083", "HYP-3091", "LTI-232", "LTI-234"),
        {
            "finite_address",
            "phi",
            "observer_gluing",
            "lean",
            "covering",
            "lmax",
            "formal",
            "loaded",
        },
        (5, 5, 5, 5, 5, 5, 5, 4),
        "the LRC packet tuple `(covering | D, 1/lmax, arc spectrum | meas)`",
        "raw time and scalar packet names unless terminal route is retained",
        "make this the receiving type for closure certificates",
    ),
    carrier(
        "observer_cut_payload_orbit",
        ("HYP-3054", "HYP-3056", "HYP-3102", "LTI-201", "LTI-203"),
        {
            "observer",
            "payload",
            "orbit",
            "gluing",
            "forgetting",
            "formal",
            "source",
            "extension",
        },
        (5, 5, 4, 5, 3, 5, 5, 5),
        "the hidden observer-extension payload over each visible fiber",
        "runner labels after replacing them by an automorphism-orbit payload",
        "attach to HYP-3098 observer gluing and HYP-3112 ear payload rows",
    ),
    carrier(
        "q27_q31_resource_descent",
        ("OPEN-Q-086", "OPEN-Q-087", "HYP-2470", "HYP-2471", "HYP-2480"),
        {
            "q27",
            "q31",
            "resource",
            "descent",
            "loaded",
            "unit_shell",
            "root_lattice",
            "state_lift",
        },
        (5, 5, 4, 4, 5, 4, 4, 4),
        "finite resource ledgers for 13-clock debt and deleted-core address",
        "direct geometry of individual arcs unless the resource row points back",
        "explain HYP-3115's finite-bank `apex7_error<=5` gate by old resource walls",
    ),
    carrier(
        "endpoint_credit_farkas",
        ("OPEN-Q-108 S156", "HYP-2970", "HYP-2995", "HYP-2997"),
        {
            "endpoint",
            "credit",
            "farkas",
            "dual",
            "observer_gluing",
            "cycle",
            "certificate",
        },
        (4, 4, 4, 5, 3, 4, 4, 3),
        "endpoint credit graph and directed-cycle winding obstruction",
        "root/payload curves if dual potentials are used alone",
        "turn proof-circuit gates into dual certificates or named directed cycles",
    ),
    carrier(
        "twist_ladder_dual",
        ("OPEN-Q-108 S155", "HYP-2972", "THM-575", "HYP-3088"),
        {
            "twist_ladder",
            "dual",
            "farkas",
            "denominator",
            "loaded",
            "finite_grid",
            "q42",
            "q27",
        },
        (4, 4, 4, 4, 5, 3, 4, 3),
        "dynamic finite denominator ladder and dual rescue certificate",
        "absolute direct-time scale after the slow/ruler coordinate changes",
        "normalize HYP-3114 intervals and THM-575 loaded-denominator rows",
    ),
    carrier(
        "source_perspective_worry_fiber",
        ("THM-381", "HYP-1977", "HYP-2486", "HYP-3106", "LTI-008"),
        {
            "source",
            "perspective",
            "a000568",
            "observer",
            "ear",
            "payload",
            "source_deck",
            "extension",
        },
        (4, 4, 3, 5, 2, 4, 5, 5),
        "source-perspective defect as an observer-coupling fiber",
        "unrooted tournament count and raw A000568 equality",
        "supply the missing observer coordinate for Lee-Yang one-runner ears",
    ),
    carrier(
        "endpoint_circuit_phi",
        ("HYP-2108", "HYP-2112", "THM-398", "HYP-3111", "HYP-3116"),
        {
            "endpoint",
            "circuit",
            "phi",
            "lean",
            "observer_gluing",
            "finite",
            "proof_circuit",
            "missing_input",
        },
        (4, 4, 4, 4, 2, 4, 3, 3),
        "exact endpoint-cover circuit positivity through the `Phi` functional",
        "global root/ear dynamics unless coupled to a payload ledger",
        "compile HYP-3115/HYP-3116/HYP-3117 proof-circuit sidecars into endpoint `Phi` gates",
    ),
    carrier(
        "proof_circuit_past_work_compiler",
        ("HYP-3117", "HYP-3116", "HYP-3102", "HYP-3107", "HYP-2108", "HYP-2112"),
        {
            "proof_circuit",
            "missing_input",
            "finite_address",
            "observer_gluing",
            "phi",
            "endpoint",
            "lean",
            "formal",
            "payload",
            "cocycle",
            "uniformity",
            "extension",
        },
        (5, 5, 4, 5, 2, 5, 5, 4),
        "fixed proof-state compiler with missing-input vectors and terminal exits",
        "raw row identity and scalar classifier labels after conversion to proof gates",
        "turn HYP-3117's gate map into required columns for the HYP-3120 packet receiver",
    ),
    carrier(
        "circuit_certificate_vector",
        ("HYP-3116", "HYP-2112", "HYP-2108", "HYP-3023", "HYP-3077", "HYP-3082"),
        {
            "proof_circuit",
            "missing_input",
            "phi",
            "endpoint",
            "finite_address",
            "observer_gluing",
            "route_purity",
            "bridge_safety",
            "horn_closure",
            "uniformity",
            "formal",
            "certificate",
            "extension",
        },
        (5, 5, 4, 5, 3, 5, 5, 5),
        "HYP-3116 circuit certificate vector: gate basis, sidecar closure, exact gap, route purity, bridge safety, uniform parameter, terminal exit",
        "endpoint owners, magnitude cocycle, root/ear payload, and uniform-family debt if the vector is scalarized",
        "add Phi/P/magnitude/Horn/protected-bridge columns before any HYP-3120 carrier can close a packet",
    ),
    carrier(
        "coordinate_resurrection_sheaf",
        ("HYP-3118", "HYP-3117", "HYP-3116", "HYP-3102", "HYP-3098", "HYP-3083"),
        {
            "coordinate_resurrection",
            "sheaf",
            "adjoint",
            "missing_input",
            "finite_address",
            "observer_gluing",
            "payload",
            "extension",
            "formal",
            "lean",
            "certificate",
        },
        (5, 5, 3, 5, 3, 4, 5, 5),
        "HYP-3118 base stalk plus live-section repair covers for destroyed proof coordinates",
        "the original quotient coordinate once replaced by a repair cover, adjoint section, or named debt",
        "attach coordinate_resurrection_cover, repair_cover_rank, live_section_type, and theorem-to-signal fields from S267",
    ),
    carrier(
        "ostrowski_beatty_pell_normal_form",
        ("HYP-1901", "HYP-1902", "HYP-2456", "HYP-3075"),
        {
            "ostrowski",
            "beatty",
            "pell",
            "continued_fraction",
            "approximation",
            "carry",
            "denominator",
            "loaded",
        },
        (3, 3, 4, 4, 5, 3, 4, 3),
        "continued-fraction carry word with Pell remainder and denominator spikes",
        "endpoint owner and status route if the carry word is scalarized",
        "add carry/Pell walls to HYP-3114 normalized witness intervals",
    ),
    carrier(
        "normal_fan_cech_barcode",
        ("HYP-3101", "LTI-240", "HYP-3025", "HYP-3018"),
        {
            "normal_fan",
            "cech",
            "barcode",
            "components",
            "topology",
            "normalized_interval",
            "finite",
        },
        (4, 4, 4, 4, 2, 3, 4, 3),
        "component-bound topology with active chamber and owner support",
        "arithmetic denominator debt unless exact scale is reattached",
        "bound normalized interval components before invoking finite grids",
    ),
    carrier(
        "signed_polymer_dirichlet_network",
        ("HYP-3073", "LTI-220", "HYP-3111", "HYP-3115"),
        {
            "dirichlet",
            "schur",
            "energy",
            "polymer",
            "f7",
            "state_lift",
            "ising",
            "network",
        },
        (4, 4, 3, 4, 3, 3, 4, 4),
        "renormalized signed-packet network and Schur-complement conductance",
        "the legal/phantom boundary distinction if F7 exits are anonymous",
        "compare Ising domain walls to Dirichlet sidecar conductance",
    ),
    carrier(
        "vitali_antipoisson_width",
        ("HYP-2104", "HYP-2105", "HYP-2154", "S599"),
        {
            "vitali",
            "anti_poisson",
            "width",
            "measure",
            "cap",
            "coverage",
            "normalized_interval",
        },
        (4, 3, 2, 4, 2, 2, 5, 3),
        "correlated-tail width obstruction instead of bare measure",
        "finite route certificates unless width is packet-keyed",
        "test coverage extremality against width/correlation, not only mass",
    ),
    carrier(
        "unit_endpoint_sieve",
        ("THM-358", "THM-359", "THM-360", "THM-369"),
        {
            "unit_endpoint",
            "sieve",
            "divisor",
            "finite_grid",
            "endpoint",
            "loaded",
            "lean",
        },
        (4, 3, 5, 3, 4, 4, 3, 2),
        "initial-segment safe set exactly at unit endpoints",
        "interior witness geometry if the unit endpoint is used alone",
        "add endpoint/unit endpoint guardrails to loaded denominator rows",
    ),
    carrier(
        "sexy_prime_residue_sieve_echo",
        ("HYP-3100", "OPEN-Q-108 residue-sieve note", "sexy-prime local channel"),
        {
            "sexy_prime",
            "residue_sieve",
            "local_obstruction",
            "admissibility",
            "prime_gap",
        },
        (2, 2, 3, 2, 2, 2, 3, 1),
        "local admissibility bookkeeping for residue channels",
        "analytic prime-distribution input and LRC packet route",
        "keep as analogy only; do not treat it as an LRC proof carrier",
    ),
]


OBLIGATIONS: List[Obligation] = [
    Obligation(
        "lean_frontier_packet",
        {"lean", "formal", "finite_address", "phi", "observer_gluing", "proof_circuit"},
        {
            "predicate_retention": 3,
            "frontier_adjacency": 3,
            "finite_checkability": 3,
            "forgotten_coordinate_control": 3,
            "formalization_readiness": 4,
        },
        "compile old certificates into `FiniteAddressBranchPacket` or `ObserverGluingCertificate`",
    ),
    Obligation(
        "observer_gluing_packet",
        {"observer", "gluing", "payload", "orbit", "endpoint", "certificate"},
        {
            "predicate_retention": 3,
            "frontier_adjacency": 4,
            "forgotten_coordinate_control": 4,
            "formalization_readiness": 3,
            "extension_payload_retention": 4,
        },
        "retain the first forgotten observer payload before quotienting a row",
    ),
    Obligation(
        "q27_resource_normalizer",
        {"q27", "q31", "resource", "descent", "loaded", "unit_shell"},
        {
            "predicate_retention": 3,
            "frontier_adjacency": 4,
            "finite_checkability": 3,
            "loaded_denominator_control": 4,
            "scalar_false_positive_detection": 2,
        },
        "turn loaded-denominator and apex-7 finite gates into resource ledgers",
    ),
    Obligation(
        "normalized_witness_interval",
        {"normalized_interval", "finite_grid", "denominator", "approximation", "loaded"},
        {
            "predicate_retention": 4,
            "finite_checkability": 4,
            "loaded_denominator_control": 4,
            "scalar_false_positive_detection": 3,
        },
        "replace direct-time witness intervals by THM-565 slow/ruler intervals",
    ),
    Obligation(
        "lee_yang_ear_payload",
        {"ear", "payload", "source", "observer", "root", "extension"},
        {
            "frontier_adjacency": 4,
            "forgotten_coordinate_control": 4,
            "scalar_false_positive_detection": 3,
            "extension_payload_retention": 4,
        },
        "explain root motion by retained one-runner payload, not final roots alone",
    ),
    Obligation(
        "root_lattice_relation_wall",
        {"root_lattice", "relation", "ising", "dirichlet", "minkowski", "network"},
        {
            "frontier_adjacency": 4,
            "finite_checkability": 2,
            "forgotten_coordinate_control": 3,
            "scalar_false_positive_detection": 4,
        },
        "separate relation pressure/root strata from fitted scalar classifiers",
    ),
    Obligation(
        "coverage_extremality_cap_debt",
        {"anti_poisson", "width", "coverage", "cap", "measure"},
        {
            "predicate_retention": 3,
            "frontier_adjacency": 3,
            "forgotten_coordinate_control": 3,
            "scalar_false_positive_detection": 4,
        },
        "measure width/correlation debt before accepting a cap extremal scalar",
    ),
    Obligation(
        "f7_state_lift_exit",
        {"f7", "state_lift", "dirichlet", "network", "q27", "resource"},
        {
            "predicate_retention": 3,
            "frontier_adjacency": 3,
            "forgotten_coordinate_control": 3,
            "formalization_readiness": 3,
            "extension_payload_retention": 3,
        },
        "name legal F7/THM-572 exits instead of leaving phantom boundary atoms anonymous",
    ),
    Obligation(
        "sexy_prime_local_channel",
        {"sexy_prime", "residue_sieve", "local_obstruction", "admissibility"},
        {
            "finite_checkability": 2,
            "forgotten_coordinate_control": 2,
            "scalar_false_positive_detection": 2,
        },
        "use LRC tournament machinery only for local admissibility bookkeeping",
    ),
]


def bridge_score(carrier_: Carrier, obligation: Obligation) -> Tuple[int, int, int]:
    axis = sum(min(carrier_.axes.get(axis_name, 0), weight) for axis_name, weight in obligation.weights.items())
    overlap = len(carrier_.tags & obligation.tags)
    anchor_bonus = sum(1 for a in carrier_.anchors if a.startswith("HYP-31") or a.startswith("OPEN-Q-108"))
    return (axis + 2 * overlap + anchor_bonus, axis, overlap)


def compare(a: Carrier, b: Carrier) -> int:
    awins = 0
    bwins = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            awins += 1
        elif b.axes[axis] > a.axes[axis]:
            bwins += 1
    if awins > bwins:
        return 1
    if bwins > awins:
        return -1
    return 0


def build_tournament(carriers: Sequence[Carrier]) -> Dict[Tuple[int, int], int]:
    priority = {c.name: i for i, c in enumerate(carriers)}
    edges: Dict[Tuple[int, int], int] = {}
    for i, j in combinations(range(len(carriers)), 2):
        c = compare(carriers[i], carriers[j])
        if c == 0:
            winner = i if priority[carriers[i].name] < priority[carriers[j].name] else j
        elif c > 0:
            winner = i
        else:
            winner = j
        edges[(i, j)] = winner
    return edges


def beats(edges: Dict[Tuple[int, int], int], i: int, j: int) -> bool:
    a, b = sorted((i, j))
    return edges[(a, b)] == i


def directed_3cycles(edges: Dict[Tuple[int, int], int], n: int) -> List[Tuple[int, int, int]]:
    cycles: List[Tuple[int, int, int]] = []
    for a, b, c in combinations(range(n), 3):
        if beats(edges, a, b) and beats(edges, b, c) and beats(edges, c, a):
            cycles.append((a, b, c))
        elif beats(edges, a, c) and beats(edges, c, b) and beats(edges, b, a):
            cycles.append((a, c, b))
    return cycles


def sccs(edges: Dict[Tuple[int, int], int], n: int) -> List[List[int]]:
    graph = {i: [j for j in range(n) if i != j and beats(edges, i, j)] for i in range(n)}
    index = 0
    stack: List[int] = []
    on_stack: Set[int] = set()
    indices: Dict[int, int] = {}
    lowlink: Dict[int, int] = {}
    comps: List[List[int]] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)

        for w in graph[v]:
            if w not in indices:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in on_stack:
                lowlink[v] = min(lowlink[v], indices[w])

        if lowlink[v] == indices[v]:
            comp: List[int] = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp))

    for v in range(n):
        if v not in indices:
            strongconnect(v)
    return sorted(comps, key=lambda comp: (-len(comp), comp[0]))


def count_hamiltonian_paths(edges: Dict[Tuple[int, int], int], n: int) -> int:
    dp: Dict[Tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if count == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(edges, last, nxt):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def main() -> None:
    print("LRC14 niche past-work closure bridge")
    print("session: codex-2026-06-27")
    print()
    print("Assumption challenge:")
    print("  vertices are past-work carriers and proof obligations, not runners.")
    print("  preserved predicate: LRC14 packet validity through a witness, observer-gluing,")
    print("    finite-address, root/ear, resource-descent, or named-residual exit.")
    print("  destroyed coordinate under audit: raw time, raw A000568 count, raw root")
    print("    locus, raw H value, raw residue-sieve analogy, and anonymous F7 debt.")
    print()

    bridge_rows = []
    for c in CARRIERS:
        for o in OBLIGATIONS:
            total, axis, overlap = bridge_score(c, o)
            bridge_rows.append((total, axis, overlap, c.name, o.name, c, o))
    bridge_rows.sort(reverse=True, key=lambda row: (row[0], row[1], row[2], row[3], row[4]))

    print("Top carrier -> obligation bridges:")
    for total, axis, overlap, cname, oname, c, o in bridge_rows[:14]:
        print(f"  {total:02d}  {cname} -> {oname}  axis={axis} tag_overlap={overlap}")
        print(f"      next: {c.next_use}")
        print(f"      target: {o.theorem_shape}")
    print()

    print("Best carrier per frontier obligation:")
    for o in OBLIGATIONS:
        ranked = sorted(
            ((bridge_score(c, o), c) for c in CARRIERS),
            reverse=True,
            key=lambda row: (row[0][0], row[0][1], row[0][2], row[1].name),
        )
        score, c = ranked[0]
        print(f"  {o.name}: {c.name} score={score[0]} anchors={','.join(c.anchors[:3])}")
    print()

    edges = build_tournament(CARRIERS)
    n = len(CARRIERS)
    outdegrees = {
        i: sum(1 for j in range(n) if i != j and beats(edges, i, j))
        for i in range(n)
    }
    hist: Dict[int, int] = {}
    for degree in outdegrees.values():
        hist[degree] = hist.get(degree, 0) + 1
    cycles = directed_3cycles(edges, n)
    comps = sccs(edges, n)
    hpaths = count_hamiltonian_paths(edges, n)
    path = sorted(range(n), key=lambda i: (-outdegrees[i], CARRIERS.index(CARRIERS[i])))
    flips = sum(
        1
        for i, j in combinations(range(n), 2)
        if CARRIERS.index(CARRIERS[i]) < CARRIERS.index(CARRIERS[j]) and beats(edges, j, i)
    )

    print("Tournament Analysis over past-work carriers:")
    print(f"  carrier_count={n}")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={len(cycles)}")
    print(f"  scc_sizes={[len(comp) for comp in comps]}")
    print(f"  hamiltonian_path_count={hpaths}")
    print(f"  edge_flips_against_tie_path={flips}")
    print("  priority_path:")
    for i in path:
        c = CARRIERS[i]
        print(f"    {outdegrees[i]:02d}  {c.name}")
    print()

    print("Closure synthesis:")
    synthesis = [
        "Use finite_address_phi_tuple as the receiver type; it scores first for Lean/frontier closure.",
        "Use observer_cut_payload_orbit plus source_perspective_worry_fiber to repair Lee-Yang ear forgetting.",
        "Use circuit_certificate_vector to require Phi/P, magnitude, Horn closure, protected-bridge, uniform-family, and terminal-exit columns.",
        "Use coordinate_resurrection_sheaf to require the HYP-3118 base stalk, a live section, repair-cover rank, and terminal debt status for every forgotten coordinate.",
        "Use q27_q31_resource_descent and twist_ladder_dual to explain loaded denominators and the apex-7 finite gate.",
        "Use proof_circuit_past_work_compiler with endpoint_credit_farkas and endpoint_circuit_phi to convert missing-input vectors into dual/endpoint certificates.",
        "Use Ostrowski/Beatty/Pell and unit_endpoint_sieve only after normalized witness intervals are in hand.",
        "Use signed_polymer_dirichlet_network for Ising/F7 boundary-energy exits, not as an anonymous residual.",
        "Keep sexy_prime_residue_sieve_echo as local admissibility bookkeeping; it lacks the analytic prime-distribution input and is not an LRC closure carrier.",
    ]
    for line in synthesis:
        print(f"  - {line}")


if __name__ == "__main__":
    main()
