#!/usr/bin/env python3
"""S580: formal certificate calculus for the current n=14 LRC frontier.

This script is intentionally small.  It does not run another speed census; it
turns the newest proof fragments into explicit gates, residual classes, and a
Tournament Analysis fingerprint whose vertices are proof obligations.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations


@dataclass(frozen=True)
class Gate:
    name: str
    status: str
    output: str
    residual: str
    evidence: str
    inputs: tuple[str, ...]
    open_rank: int
    residual_weight: int
    dependency_count: int
    evidence_strength: int
    tie_order: int
    extension: str

    def hardness_key(self) -> tuple[int, int, int, int, int]:
        """Higher key means harder/more proof-critical in the current ledger."""
        return (
            self.open_rank,
            self.residual_weight,
            self.dependency_count,
            -self.evidence_strength,
            self.tie_order,
        )


@dataclass(frozen=True)
class Extension:
    name: str
    target: str
    route: str
    first_test: str
    payoff: str
    risk: int
    leverage: int
    novelty: int

    def priority_key(self) -> tuple[int, int, int]:
        return (self.leverage, self.novelty, -self.risk)


GATES = [
    Gate(
        name="G0_no_multiple_n_clock",
        status="PROVED",
        output="closed witness t=1/n",
        residual="none",
        evidence="THM-369; reused by THM-398",
        inputs=("speed row S", "no v divisible by n"),
        open_rank=0,
        residual_weight=0,
        dependency_count=1,
        evidence_strength=5,
        tie_order=0,
        extension="Keep as the first clearance gate; it removes all no-multiple rows.",
    ),
    Gate(
        name="G1_Cprime_reduction",
        status="PROVED reduction",
        output="LRC(n) follows from Cprime",
        residual="Cprime itself: n|v implies M(S)>1/n",
        evidence="THM-398 formalizes HYP-2102",
        inputs=("G0_no_multiple_n_clock", "Cprime multiple branch"),
        open_rank=0,
        residual_weight=1,
        dependency_count=2,
        evidence_strength=5,
        tie_order=1,
        extension="Use Cprime as the global spine; every other n=14 route should feed it or bypass it with a cheap witness.",
    ),
    Gate(
        name="G2_dominance_or_long_interval_dodge",
        status="PROVED partial Cprime",
        output="positive measure",
        residual="all-short Vitali alignment for a small multiple",
        evidence="THM-398 Lemma B/Bprime; HYP-2104 quantifies 96.8% proved at n=14 samples",
        inputs=("LRC(n-1)", "component of G(S\\{v}) longer than 2/(n v)"),
        open_rank=0,
        residual_weight=2,
        dependency_count=2,
        evidence_strength=4,
        tie_order=2,
        extension="Replace crude measure subtraction by component-length accounting.",
    ),
    Gate(
        name="G3_all_short_Cprime_residual",
        status="OPEN",
        output="positive measure for small multiples",
        residual="one AP of thin arcs is arc-aligned enough to cover a fixed open safe set",
        evidence="HYP-2103/HYP-2104; all-short tight rows remain 0 in S573 sample",
        inputs=("small multiple v=nw", "all components of G(S\\{v}) <= 2/(n^2 w)"),
        open_rank=3,
        residual_weight=7,
        dependency_count=3,
        evidence_strength=2,
        tie_order=9,
        extension="Translate the AP-cover residual into endpoint-owner cover circuits.",
    ),
    Gate(
        name="G4_paired_or_anchored_cheap_pair",
        status="REDUCTION plus verified split",
        output="unblocked small-pair pinch witness",
        residual="block-all rows must be proved positive measure",
        evidence="HYP-2095; THM-396/397 blocker dichotomy",
        inputs=("measure-zero row", "small reduced-sum pair ledger"),
        open_rank=2,
        residual_weight=5,
        dependency_count=3,
        evidence_strength=3,
        tie_order=5,
        extension="Make block-all equivalent to an endpoint cover circuit.",
    ),
    Gate(
        name="G5_fixed_boundary_owner_functor",
        status="BOUNDED scaffold",
        output="64 fixed fibres with speed-owner labels",
        residual="unbounded realization independence",
        evidence="HYP-2099; AP and Vstar floor rows both cheap in canonical fibre",
        inputs=("64 fixed boundary", "unit spine", "endpoint owners"),
        open_rank=2,
        residual_weight=6,
        dependency_count=4,
        evidence_strength=3,
        tie_order=7,
        extension="Attach HYP-2101 sections to each fixed-owner fibre.",
    ),
    Gate(
        name="G6_unit_lift_cheap_sieve",
        status="BOUNDED zero residual",
        output="cheap pair fires before exchange",
        residual="simultaneous multi-unit no-cheap rows",
        evidence="HYP-2100: 13169/13169 full one-unit lifts cheap",
        inputs=("unit-shell lift", "full D/U/N cover"),
        open_rank=1,
        residual_weight=3,
        dependency_count=2,
        evidence_strength=4,
        tie_order=3,
        extension="Prove full cover => cheap pair before attempting unit lowering.",
    ),
    Gate(
        name="G7_apex_lift_certificate_sheaf",
        status="BOUNDED zero restriction residual",
        output="global cheap section or positive-measure restriction",
        residual="unbounded apex monodromy over fixed fibres",
        evidence="HYP-2101: restriction_residual=0 on S579 site",
        inputs=("mod-7 tie wall", "unit-shell lift site", "local sections"),
        open_rank=1,
        residual_weight=4,
        dependency_count=4,
        evidence_strength=4,
        tie_order=4,
        extension="Search for residual acyclicity of section transports.",
    ),
    Gate(
        name="G8_endpoint_cover_circuit_positivity",
        status="OPEN bridge",
        output="positive measure from failed gluing",
        residual="cover circuits with no private pivot/owner",
        evidence="HYP-2095 + THM-397 + HYP-2101 synthesis",
        inputs=("block-all cheap pairs", "shield/anchor ledger", "D/U/N pivots"),
        open_rank=3,
        residual_weight=8,
        dependency_count=4,
        evidence_strength=2,
        tie_order=10,
        extension="This is the natural common target of all-short Cprime and failed sheaf gluing.",
    ),
    Gate(
        name="G9_tie_wall_limit_functor",
        status="OPEN formalization",
        output="boundary perturbations become a labelled tie-wall object",
        residual="AP/Vstar perturbation-direction ambiguity",
        evidence="HYP-2098 and HYP-2101",
        inputs=("mod-7 lanes", "apex seam", "left/right perturbation directions"),
        open_rank=2,
        residual_weight=5,
        dependency_count=3,
        evidence_strength=2,
        tie_order=6,
        extension="Treat time-wall limits as sheaf stalks, not as one naked tournament class.",
    ),
    Gate(
        name="G10_certificate_calculus_closure",
        status="PROPOSED",
        output="every row has a section or a named positive-measure residual",
        residual="terminal residual acyclicity",
        evidence="S580 synthesis of THM-398, HYP-2095, HYP-2099-HYP-2103",
        inputs=("G0-G9",),
        open_rank=3,
        residual_weight=9,
        dependency_count=5,
        evidence_strength=1,
        tie_order=11,
        extension="Prove no terminal residual can be closed under restrictions.",
    ),
]


EXTENSIONS = [
    Extension(
        name="E1_all_short_to_endpoint_cover",
        target="HYP-2103/HYP-2104 all-short Vitali alignment residual",
        route="turn thin AP cover of G(S') into an endpoint-owner cover circuit",
        first_test="enumerate all-short v=14 rows and record endpoint owners of every covered component",
        payoff="would prove the remaining Cprime branch for n=14 if every cover circuit peels",
        risk=2,
        leverage=5,
        novelty=4,
    ),
    Extension(
        name="E2_fixed_fibre_sheaf_table",
        target="HYP-2099/HYP-2101 bridge",
        route="attach local section data to AP, Vstar, transversal flips, and minimal gcd-3/gcd-9 lifts",
        first_test="build the promised S580 table: cheap pair, reduced sum, mod-7 owner, mod-27 shell, blockers, pivots",
        payoff="turns the 64-class overcount into labelled proof obligations",
        risk=2,
        leverage=4,
        novelty=3,
    ),
    Extension(
        name="E3_residual_acyclicity",
        target="apex-lift monodromy",
        route="make a directed graph of certificate transports and prove terminal SCCs are positive measure",
        first_test="use S579 restrictions; add union-only rows as two-step composites",
        payoff="formal sheaf gluing proof without class-by-class guessing",
        risk=3,
        leverage=4,
        novelty=4,
    ),
    Extension(
        name="E4_component_discrepancy_bound",
        target="HYP-2104 Criterion Bprime all-short complement",
        route="Erdos-Turan or three-distance bound for {k/(nw)} against G(S') components",
        first_test="measure component count and endpoint denominator spectrum in bounded small-multiple rows",
        payoff="analytic proof of Cprime independent of the 64-class machinery",
        risk=4,
        leverage=5,
        novelty=3,
    ),
    Extension(
        name="E5_level3_cascade_product",
        target="certificate product closure",
        route="model each gate as a conditional clearance factor in a cascade product",
        first_test="assign each S579 route a factor type: witness, cheap, positive_measure, ledger_failure, residual",
        payoff="connects the user's cascade language to a checkable proof ledger",
        risk=2,
        leverage=3,
        novelty=3,
    ),
]


def beats(a: Gate, b: Gate) -> bool:
    return a.hardness_key() > b.hardness_key()


def tournament_edges(gates: list[Gate]) -> dict[tuple[int, int], int]:
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(gates)), 2):
        if beats(gates[i], gates[j]):
            edges[(i, j)] = 1
            edges[(j, i)] = 0
        else:
            edges[(i, j)] = 0
            edges[(j, i)] = 1
    return edges


def sccs(gates: list[Gate], edges: dict[tuple[int, int], int]) -> list[list[str]]:
    n = len(gates)
    graph = [[] for _ in range(n)]
    rev = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and edges[(i, j)]:
                graph[i].append(j)
                rev[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    comps: list[list[str]] = []
    seen.clear()

    def rdfs(v: int, comp: list[int]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            comps.append([gates[i].name for i in sorted(comp)])
    return comps


def directed_3_cycles(gates: list[Gate], edges: dict[tuple[int, int], int]) -> int:
    count = 0
    for i, j, k in combinations(range(len(gates)), 3):
        eij = edges[(i, j)]
        ejk = edges[(j, k)]
        eki = edges[(k, i)]
        if eij and ejk and eki:
            count += 1
        if (not eij) and (not ejk) and (not eki):
            count += 1
    return count


def hamiltonian_path_count(gates: list[Gate], edges: dict[tuple[int, int], int]) -> int:
    n = len(gates)

    @lru_cache(None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if prev_mask & (1 << prev) and edges[(prev, last)]:
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def print_gate_table() -> None:
    print("GATES")
    print("name | status | output | residual")
    for gate in GATES:
        print(f"{gate.name} | {gate.status} | {gate.output} | {gate.residual}")
    print()


def print_formal_calculus() -> None:
    print("S580 LRC n=14 certificate calculus")
    print()
    print("FORMAL OBJECTS")
    print("SpeedRow := (n, S, time-site labels, owner labels, quotient labels).")
    print("CertificateSection := witness_1_over_n | cheap_pair | positive_measure | ledger_failure | residual.")
    print("Gate := a partial rule SpeedRow -> CertificateSection or NamedResidual, with restriction maps between rows.")
    print("Vitali handoff := the divisibility wall n|v, where THM-398 moves from construction to measure.")
    print("Cascade product := ordered product of conditional clearances supplied by gates G0..G10.")
    print("Target theorem form := every row exits through a witness/positive_measure section, or every terminal residual SCC is impossible.")
    print()


def print_tournament() -> None:
    edges = tournament_edges(GATES)
    scores = []
    for i, gate in enumerate(GATES):
        score = sum(edges[(i, j)] for j in range(len(GATES)) if i != j)
        scores.append((score, gate.name))
    scores_sorted = sorted(scores, reverse=True)
    score_hist = Counter(score for score, _ in scores)

    print("TOURNAMENT ANALYSIS")
    print("vertices: proof obligations/gates, not runners and not naked round classes")
    print("pair observable: (open_rank, residual_weight, dependency_count, -evidence_strength, tie_order)")
    print("switch: harder unresolved proof obligation beats easier certified gate")
    print(f"score_hist: {dict(sorted(score_hist.items()))}")
    print(f"directed_3_cycles: {directed_3_cycles(GATES, edges)}")
    print(f"sccs: {sccs(GATES, edges)}")
    print(f"hamiltonian_path_count: {hamiltonian_path_count(GATES, edges)}")
    print("hardness_path:")
    for score, name in scores_sorted:
        print(f"  {score:2d}  {name}")
    print()


def print_residuals() -> None:
    residual_groups: dict[str, list[str]] = defaultdict(list)
    for gate in GATES:
        residual_groups[gate.residual].append(gate.name)

    print("NAMED RESIDUALS")
    for residual, names in residual_groups.items():
        if residual != "none":
            print(f"- {residual}: {', '.join(names)}")
    print()


def print_extensions() -> None:
    print("EXTENSION QUEUE")
    for ext in sorted(EXTENSIONS, key=lambda e: e.priority_key(), reverse=True):
        print(f"{ext.name}: target={ext.target}")
        print(f"  route: {ext.route}")
        print(f"  first_test: {ext.first_test}")
        print(f"  payoff: {ext.payoff}")
        print(f"  priority_key: {ext.priority_key()}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("Considered vertex sets: runners, gaps, fixed sections, wall events, residues, cover arcs, Fourier modes, certificate germs, and proof obligations.")
    print("Chosen vertices here: proof obligations/gates.")
    print("Preserved predicate: whether every row receives a witness, cheap section, positive-measure section, or named residual.")
    print("Destroyed data: exact speed values, full round-class identity, and most Hamiltonian-path complexity.")
    print("Challenged assumption: the n=14 proof does not have to be a 64-class table first; the table is useful only after owner labels and section restrictions are attached.")
    print()


def main() -> None:
    print_formal_calculus()
    print_gate_table()
    print_residuals()
    print_tournament()
    print_extensions()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
