#!/usr/bin/env python3
"""Circuit-complexity synthesis over older LRC14 proof carriers.

This is a small executable literature lookback inside the repository.  It does
not recompute the HYP-2963 packet bank.  Instead it records audited facts from
older hypotheses, gives them common proof-circuit fields, and runs Tournament
Analysis on carriers rather than runner rows.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Tuple


@dataclass(frozen=True)
class Carrier:
    name: str
    sources: Tuple[str, ...]
    role: str
    predicate: int
    exact_gap: int
    uniformity: int
    route_purity: int
    sidecar_closure: int
    bridge_safety: int
    proof_ready: int
    destroyed_discharge: int
    finite_bank_warning: int
    size_proxy: int
    depth_proxy: int
    retained_coordinate: str
    destroyed_coordinate: str
    next_test: str

    def proof_vector(self) -> Tuple[int, ...]:
        return (
            self.predicate,
            self.exact_gap,
            self.uniformity,
            self.route_purity,
            self.sidecar_closure,
            self.bridge_safety,
            self.proof_ready,
            self.destroyed_discharge,
            -self.finite_bank_warning,
        )

    def weighted_payload(self) -> int:
        weights = (4, 4, 3, 3, 3, 3, 4, 3, 3)
        return sum(w * v for w, v in zip(weights, self.proof_vector()))

    def small_circuit_key(self) -> Tuple[int, int, int]:
        return (self.size_proxy, self.depth_proxy, self.finite_bank_warning)


TIE_PATH = [
    "endpoint_phi_relu_gap",
    "endpoint_cover_circuit_positivity",
    "automatic_magnitude_zipper",
    "route_state_horn_median_hull",
    "protected_branch_graph_no_naked_bridge",
    "pde_weak_form_endpoint_compiler",
    "mciq_monotone_proof_frontier",
    "three_state_middle_automaton",
    "finite_bank_apex_threshold_warning",
    "raw_scalar_p0_shadow",
]

TIE_RANK = {name: rank for rank, name in enumerate(TIE_PATH)}


CARRIERS = [
    Carrier(
        name="endpoint_phi_relu_gap",
        sources=("HYP-2112", "HYP-2108"),
        role="exact numeric output gate",
        predicate=5,
        exact_gap=5,
        uniformity=4,
        route_purity=3,
        sidecar_closure=3,
        bridge_safety=2,
        proof_ready=5,
        destroyed_discharge=4,
        finite_bank_warning=0,
        size_proxy=4,
        depth_proxy=2,
        retained_coordinate="Phi(C)=mu(safe set) as a ReLU/ramp gap sum",
        destroyed_coordinate="which endpoint owners caused each ramp unless paired with HYP-2108",
        next_test="compile Phi, endpoint owners, and multiple-of-14 checks into HYP-2963 rows",
    ),
    Carrier(
        name="endpoint_cover_circuit_positivity",
        sources=("HYP-2108",),
        role="closed endpoint-cover sign gate",
        predicate=5,
        exact_gap=4,
        uniformity=4,
        route_purity=3,
        sidecar_closure=2,
        bridge_safety=2,
        proof_ready=4,
        destroyed_discharge=4,
        finite_bank_warning=0,
        size_proxy=3,
        depth_proxy=2,
        retained_coordinate="P(S)>0 iff loose; midpoint circuit resonance",
        destroyed_coordinate="gap magnitude; recovered by HYP-2112 Phi",
        next_test="pair P-sign with Phi-value so tight/worry kernels become explicit proof obligations",
    ),
    Carrier(
        name="three_state_middle_automaton",
        sources=("HYP-2109",),
        role="stateful wall-crossing gate",
        predicate=4,
        exact_gap=2,
        uniformity=3,
        route_purity=3,
        sidecar_closure=3,
        bridge_safety=2,
        proof_ready=3,
        destroyed_discharge=4,
        finite_bank_warning=0,
        size_proxy=3,
        depth_proxy=2,
        retained_coordinate="L/M/R terminal state, with M as wall or tie corridor",
        destroyed_coordinate="metric gap size and magnitude cocycle",
        next_test="prove or refute no-closed-middle around endpoint-cover circuits",
    ),
    Carrier(
        name="automatic_magnitude_zipper",
        sources=("HYP-3016", "HYP-3023"),
        role="finite-state quotient repaired by magnitude cocycle",
        predicate=4,
        exact_gap=3,
        uniformity=4,
        route_purity=5,
        sidecar_closure=4,
        bridge_safety=3,
        proof_ready=4,
        destroyed_discharge=5,
        finite_bank_warning=1,
        size_proxy=6,
        depth_proxy=3,
        retained_coordinate="(M, q_threshold, farey_excess, lacunary_tail_ratio)",
        destroyed_coordinate="route purity inside automatic/residue-terminal fibers",
        next_test="turn the magnitude cocycle purity into a family lemma for automatic fibers",
    ),
    Carrier(
        name="route_state_horn_median_hull",
        sources=("HYP-3074", "HYP-3077"),
        role="proof-circuit legality and closure gate",
        predicate=5,
        exact_gap=2,
        uniformity=4,
        route_purity=4,
        sidecar_closure=5,
        bridge_safety=3,
        proof_ready=4,
        destroyed_discharge=5,
        finite_bank_warning=0,
        size_proxy=34,
        depth_proxy=2,
        retained_coordinate="packet/route/certificate/sidecar/discharge Horn closure",
        destroyed_coordinate="specific terminal atom when only scheduler center survives",
        next_test="compile multi-premise sidecar laws or prove the current unary closure is enough",
    ),
    Carrier(
        name="protected_branch_graph_no_naked_bridge",
        sources=("HYP-3082",),
        role="terminal no-naked-bridge gate",
        predicate=5,
        exact_gap=2,
        uniformity=4,
        route_purity=5,
        sidecar_closure=4,
        bridge_safety=5,
        proof_ready=4,
        destroyed_discharge=5,
        finite_bank_warning=0,
        size_proxy=83,
        depth_proxy=4,
        retained_coordinate="protected branch graph after route/section/grid/no-lift/q-cusp exits",
        destroyed_coordinate="raw scalar-star bridges; five naked bridges before protection",
        next_test="prove every primitive residual emits the protected sidecar set before contraction",
    ),
    Carrier(
        name="mciq_monotone_proof_frontier",
        sources=("HYP-3111", "HYP-3115"),
        role="bounded monotone proof-circuit frontier",
        predicate=4,
        exact_gap=3,
        uniformity=3,
        route_purity=3,
        sidecar_closure=4,
        bridge_safety=3,
        proof_ready=4,
        destroyed_discharge=4,
        finite_bank_warning=1,
        size_proxy=8,
        depth_proxy=4,
        retained_coordinate="10-input monotone proof-frontier circuit, size 8, depth 4",
        destroyed_coordinate="uniform family parameter and non-fitted thresholds",
        next_test="replace finite-bank literals with a uniform circuit over named proof sidecars",
    ),
    Carrier(
        name="pde_weak_form_endpoint_compiler",
        sources=("HYP-3111", "HYP-2112", "HYP-2108"),
        role="weak-form compiler into endpoint Phi",
        predicate=4,
        exact_gap=4,
        uniformity=3,
        route_purity=3,
        sidecar_closure=4,
        bridge_safety=3,
        proof_ready=4,
        destroyed_discharge=4,
        finite_bank_warning=0,
        size_proxy=8,
        depth_proxy=5,
        retained_coordinate="mass/stiffness/boundary/zero-mode data routed to Phi",
        destroyed_coordinate="root/free-energy geometry before the weak-form packet is typed",
        next_test="annotate residual rows with mass, stiffness, zero modes, Phi, and terminal exit",
    ),
    Carrier(
        name="finite_bank_apex_threshold_warning",
        sources=("HYP-3115",),
        role="small fitted classifier warning",
        predicate=1,
        exact_gap=1,
        uniformity=1,
        route_purity=5,
        sidecar_closure=1,
        bridge_safety=0,
        proof_ready=1,
        destroyed_discharge=1,
        finite_bank_warning=5,
        size_proxy=1,
        depth_proxy=1,
        retained_coordinate="single bounded-bank literal apex7_error <= 5",
        destroyed_coordinate="uniformity, proof basis, and reason the threshold is not fitted",
        next_test="use only as an alarm that a genuine uniform circuit field is missing",
    ),
    Carrier(
        name="raw_scalar_p0_shadow",
        sources=("HYP-3103", "HYP-3115"),
        role="negative-control scalar output",
        predicate=1,
        exact_gap=1,
        uniformity=1,
        route_purity=1,
        sidecar_closure=0,
        bridge_safety=0,
        proof_ready=1,
        destroyed_discharge=0,
        finite_bank_warning=4,
        size_proxy=1,
        depth_proxy=1,
        retained_coordinate="single value p0=G(0)",
        destroyed_coordinate="root locus, route label, endpoint owners, and all proof sidecars",
        next_test="never use without a retained root/ear/gap/route certificate",
    ),
]


PAST_WORK_FACTS = [
    ("HYP-2112", "Phi(C)=mu(safe set) verified exactly 900/900 for n=6..14; Phi>0 iff loose."),
    ("HYP-2108", "P(S)>0 iff loose for endpoint-cover circuits; sign gate precedes Phi value."),
    ("HYP-2109", "L/M/R automaton proposes M as the unavoidable wall/tie corridor."),
    ("HYP-3016", "Automatic/residue languages mix AP/GW boundary atoms with open rows unless magnitude is retained."),
    ("HYP-3023", "Full HYP-2963 bank: magnitude_cocycle has 100.0% route purity after automatic fibers mix."),
    ("HYP-3077", "Horn proof-state hull: 41 features, 34 rules, 29791 triples, 0 illegal centers."),
    ("HYP-3082", "Raw scalar star has 5 naked bridges; protected branch graph has 0 naked bridges."),
    ("HYP-3111", "Monotone proof-frontier circuit has 10 inputs, size 8, depth 4, all inputs essential."),
    ("HYP-3115", "Finite-bank literal apex7_error<=5 isolates max p0, but is explicitly a uniformity warning."),
]


def majority_compare(a: Carrier, b: Carrier) -> int:
    votes = 0
    for av, bv in zip(a.proof_vector(), b.proof_vector()):
        if av > bv:
            votes += 1
        elif av < bv:
            votes -= 1
    if votes:
        return 1 if votes > 0 else -1
    if a.weighted_payload() != b.weighted_payload():
        return 1 if a.weighted_payload() > b.weighted_payload() else -1
    return 1 if TIE_RANK[a.name] < TIE_RANK[b.name] else -1


def small_circuit_compare(a: Carrier, b: Carrier) -> int:
    if a.small_circuit_key() != b.small_circuit_key():
        return 1 if a.small_circuit_key() < b.small_circuit_key() else -1
    return 1 if TIE_RANK[a.name] < TIE_RANK[b.name] else -1


def adjacency(
    carriers: List[Carrier], compare
) -> Dict[Tuple[str, str], bool]:
    edges: Dict[Tuple[str, str], bool] = {}
    for a, b in combinations(carriers, 2):
        a_wins = compare(a, b) > 0
        edges[(a.name, b.name)] = a_wins
        edges[(b.name, a.name)] = not a_wins
    return edges


def outdegree_hist(carriers: List[Carrier], edges: Dict[Tuple[str, str], bool]) -> Counter:
    hist: Counter = Counter()
    for a in carriers:
        score = sum(1 for b in carriers if a != b and edges[(a.name, b.name)])
        hist[score] += 1
    return hist


def directed_three_cycles(carriers: List[Carrier], edges: Dict[Tuple[str, str], bool]) -> int:
    count = 0
    for a, b, c in combinations([x.name for x in carriers], 3):
        if (
            edges[(a, b)]
            and edges[(b, c)]
            and edges[(c, a)]
            or edges[(a, c)]
            and edges[(c, b)]
            and edges[(b, a)]
        ):
            count += 1
    return count


def scc_sizes(carriers: List[Carrier], edges: Dict[Tuple[str, str], bool]) -> List[int]:
    names = [c.name for c in carriers]
    graph = {a: [b for b in names if a != b and edges[(a, b)]] for a in names}
    index = 0
    stack: List[str] = []
    on_stack = set()
    indices: Dict[str, int] = {}
    low: Dict[str, int] = {}
    sizes: List[int] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in graph[v]:
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for name in names:
        if name not in indices:
            strongconnect(name)
    return sizes


def hamiltonian_path_count(carriers: List[Carrier], edges: Dict[Tuple[str, str], bool]) -> int:
    names = [c.name for c in carriers]
    n = len(names)
    dp = [defaultdict(int) for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last, count in list(dp[mask].items()):
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edges[(names[last], names[nxt])]:
                    dp[mask | (1 << nxt)][nxt] += count
    full = (1 << n) - 1
    return sum(dp[full].values())


def proof_order(carriers: List[Carrier], edges: Dict[Tuple[str, str], bool]) -> List[str]:
    return [
        c.name
        for c in sorted(
            carriers,
            key=lambda x: (
                -sum(1 for y in carriers if x != y and edges[(x.name, y.name)]),
                TIE_RANK[x.name],
            ),
        )
    ]


def edge_flips(
    carriers: List[Carrier],
    left: Dict[Tuple[str, str], bool],
    right: Dict[Tuple[str, str], bool],
) -> List[Tuple[str, str]]:
    flips = []
    for a, b in combinations([c.name for c in carriers], 2):
        if left[(a, b)] != right[(a, b)]:
            flips.append((a, b))
    return flips


def format_hist(hist: Counter) -> str:
    return "{" + ", ".join(f"{k}:{hist[k]}" for k in sorted(hist)) + "}"


def main() -> None:
    proof_edges = adjacency(CARRIERS, majority_compare)
    circuit_edges = adjacency(CARRIERS, small_circuit_compare)
    flips = edge_flips(CARRIERS, proof_edges, circuit_edges)

    print("LRC14 CIRCUIT-COMPLEXITY PAST-WORK SYNTHESIS")
    print("status: evidence / executable synthesis; not a proof")
    print("vertex_set: proof carriers and sidecar gates, not runners or arcs")
    print(
        "pairwise_observable: retained LRC predicate, exact gap readout, uniform family, "
        "route purity, sidecar closure, bridge safety, proof readiness, destroyed-coordinate "
        "discharge, and finite-bank warning"
    )
    print("switch: majority comparison of the observable vector; ties use weighted payload then the declared tie path")
    print("tie_hamiltonian_path:")
    print("  " + " -> ".join(TIE_PATH))
    print()

    print("PAST-WORK FACTS")
    for source, fact in PAST_WORK_FACTS:
        print(f"- {source}: {fact}")
    print()

    print("CARRIER LEDGER")
    for carrier in CARRIERS:
        print(f"- {carrier.name}")
        print(f"  sources: {', '.join(carrier.sources)}")
        print(f"  role: {carrier.role}")
        print(f"  payload_score: {carrier.weighted_payload()}")
        print(f"  circuit_proxy: size={carrier.size_proxy}, depth={carrier.depth_proxy}")
        print(f"  retains: {carrier.retained_coordinate}")
        print(f"  destroys_or_warns: {carrier.destroyed_coordinate}")
        print(f"  next_test: {carrier.next_test}")
    print()

    print("TOURNAMENT FINGERPRINT")
    print(f"score_hist={format_hist(outdegree_hist(CARRIERS, proof_edges))}")
    print(f"directed_3cycles={directed_three_cycles(CARRIERS, proof_edges)}")
    print(f"scc_sizes={scc_sizes(CARRIERS, proof_edges)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(CARRIERS, proof_edges)}")
    print("proof_payload_order:")
    print("  " + " > ".join(proof_order(CARRIERS, proof_edges)))
    print()

    print("SMALLEST-CIRCUIT GAUGE CHECK")
    print(f"edge_flips_against_smallest_circuit_first={len(flips)}")
    for a, b in flips:
        winner = a if proof_edges[(a, b)] else b
        small = a if circuit_edges[(a, b)] else b
        print(f"- proof_payload prefers {winner}; smallest_circuit_first prefers {small}")
    print()

    print("SYNTHESIZED INVARIANT")
    print(
        "circuit_certificate_vector = (input_packet_schema, gate_basis, sidecar_closure, "
        "exact_gap_functional, route_purity, bridge_safety, uniform_family_parameter, terminal_exit)"
    )
    print("gate_assignment:")
    print("- output_numeric_gate: HYP-2112 Phi")
    print("- sign_resonance_gate: HYP-2108 endpoint-cover P")
    print("- stateful_wall_gate: HYP-2109 L/M/R middle automaton")
    print("- route_purity_gate: HYP-3023 magnitude cocycle")
    print("- legality_gate: HYP-3077 Horn sidecar closure")
    print("- terminal_bridge_gate: HYP-3082 protected branch graph")
    print("- frontier_uniformity_gate: HYP-3111/HYP-3115 proof circuit")
    print()

    print("ASSUMPTION CHALLENGE")
    print("- challenged: circuit vertices should be runner rows or Boolean literals.")
    print("- replacement: vertices should be proof obligations/gates with declared retained coordinates.")
    print("- preserves: LRC14 packet validity up to exact gap, route, bridge, and terminal-exit sidecars.")
    print("- destroys if omitted: endpoint owners, magnitude cocycle, root/ear payload, and uniform-family debt.")
    print()

    print("NEXT PROOF TESTS")
    print("1. Build a HYP-2963 row ledger with columns Phi, P-sign, L/M/R terminal state, magnitude cocycle, Horn closure, protected-bridge status, and proof-depth stage.")
    print("2. Prove the no-closed-middle lemma or exhibit the first endpoint-cover circuit whose terminal M cycle survives all Phi and magnitude sidecars.")
    print("3. Replace the finite apex7_error<=5 literal with a uniform circuit family over named sidecar gates; reject it if size/depth grows by fitting the bank.")
    print("4. Use protected-branch status as the terminal gate: a row closes only after no naked bridge remains or a named THM-572/F7 residual is emitted.")


if __name__ == "__main__":
    main()
