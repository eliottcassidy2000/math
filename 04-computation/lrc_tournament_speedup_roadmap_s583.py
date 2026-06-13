#!/usr/bin/env python3
"""Roadmap of tournament computations that can speed up LRC proof search.

This is a synthesis script, not a new verifier.  It turns repo archaeology into
an ordered list of reusable engines, with a Tournament Analysis fingerprint on
the engines themselves.
"""

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations


@dataclass(frozen=True)
class Engine:
    name: str
    certificate_safety: int
    algorithmic_leverage: int
    proof_payoff: int
    implementation_cost: int
    maturity: int
    source: str
    preserves: str
    destroys: str
    next_move: str

    def key(self):
        return (
            self.certificate_safety,
            self.proof_payoff,
            self.algorithmic_leverage,
            -self.implementation_cost,
            self.maturity,
        )

    def priority(self):
        return (
            3 * self.certificate_safety
            + 2 * self.proof_payoff
            + 2 * self.algorithmic_leverage
            + self.maturity
            - self.implementation_cost
        )


ENGINES = [
    Engine(
        name="observer_source_marked_reachability",
        certificate_safety=5,
        algorithmic_leverage=4,
        proof_payoff=5,
        implementation_cost=2,
        maturity=5,
        source="THM-381, lrc_observer_source_tournament_s511.py",
        preserves="exact LRC predicate: observer lonely iff marked observer is a source",
        destroys="continuous time geometry after cell/class projection",
        next_move="Build a marked-class transition automaton with source targets precomputed.",
    ),
    Engine(
        name="Cprime_gate_tournament_calculus",
        certificate_safety=5,
        algorithmic_leverage=5,
        proof_payoff=5,
        implementation_cost=3,
        maturity=5,
        source="THM-398, HYP-2105/2108/2112, lrc_n14_certificate_calculus_s580.py, lrc_circuit_to_gap_functional_s576.py",
        preserves="proof obligations, residual gates, Phi/P gap functionals, and named exits",
        destroys="fine speed realization labels unless endpoint owners are carried",
        next_move="Compile gates into a decision DAG whose terminal middle circuits call Phi(C)/P(S).",
    ),
    Engine(
        name="A2_iso_fingerprint_lookup",
        certificate_safety=3,
        algorithmic_leverage=5,
        proof_payoff=3,
        implementation_cost=2,
        maturity=4,
        source="fractal-isomorphism-via-a2.md, a2_applications_s339.py",
        preserves="tournament isomorphism class through n<=9; strong hash beyond",
        destroys="proof certainty beyond verified range unless backed by canonical fallback",
        next_move="Use sorted A^2 rows as first-pass cache key for LRC quotient walks.",
    ),
    Engine(
        name="goodcut_SCC_protection_core",
        certificate_safety=4,
        algorithmic_leverage=4,
        proof_payoff=4,
        implementation_cost=2,
        maturity=5,
        source="THM-354, THM-395, good-cut-count.md, tournament_tda.py",
        preserves="endpoint-protection core shape via SCC defect and backward-wedge debt",
        destroys="arithmetic labels on protector intervals unless attached separately",
        next_move="Attach SCC/private-leaf peeling to endpoint-owner and pressure graphs.",
    ),
    Engine(
        name="three_state_middle_automata",
        certificate_safety=3,
        algorithmic_leverage=4,
        proof_payoff=5,
        implementation_cost=3,
        maturity=3,
        source="HYP-2109, tournament_three_state_automaton_s582.py",
        preserves="left/right owner exits and middle wall survival before projection",
        destroys="metric margin unless M cells keep length/residue labels",
        next_move="Make endpoint-owner and residue cells emit L/M/R event words.",
    ),
    Engine(
        name="threshold_decorated_quotient_stack",
        certificate_safety=4,
        algorithmic_leverage=4,
        proof_payoff=4,
        implementation_cost=3,
        maturity=4,
        source="lrc_restricted_tournament_mappings_s535.py, lrc_section_boundary_functors_s539.py",
        preserves="good-only class fibers after observer/source/gap/threshold labels",
        destroys="purity if labels are dropped; raw phase classes mix good and bad",
        next_move="Use good-only colored fibers as fast certificate targets before exact search.",
    ),
    Engine(
        name="cheap_TDA_prefilter",
        certificate_safety=2,
        algorithmic_leverage=4,
        proof_payoff=3,
        implementation_cost=1,
        maturity=5,
        source="tournament_tda.py, tournament_tda_residue_features_s355.out",
        preserves="score, c3, SCC defect, residue rank, Omega packet summaries",
        destroys="the LRC threshold predicate; cannot certify alone",
        next_move="Emit LRC JSON fingerprints for each quotient cell and rank residuals.",
    ),
    Engine(
        name="near_transitive_perturbation_DP",
        certificate_safety=2,
        algorithmic_leverage=5,
        proof_payoff=3,
        implementation_cost=3,
        maturity=3,
        source="applications-probe-2026-05-30.md, algorithmic_breakthroughs_s90ai.out",
        preserves="small upset set relative to a base order or tie-wall perturbation",
        destroys="global class exactness when the upset set is not small",
        next_move="Parameterize AP/Vstar/tie-wall neighborhoods by upset set k and update locally.",
    ),
    Engine(
        name="Burnside_restricted_orbit_counts",
        certificate_safety=3,
        algorithmic_leverage=5,
        proof_payoff=4,
        implementation_cost=4,
        maturity=4,
        source="HYP-2097, HYP-2099, A000568/Burnside quotient scripts",
        preserves="symmetry quotient counts and orbit representatives",
        destroys="endpoint/pair certificates until labels are reattached",
        next_move="Refine A000568(13)->round316->merged190->fixed64 with owner labels.",
    ),
]


def build_tournament():
    names = [e.name for e in ENGINES]
    matrix = {u: {v: False for v in names if v != u} for u in names}
    for a, b in combinations(ENGINES, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        matrix[winner.name][loser.name] = True
        matrix[loser.name][winner.name] = False
    return names, matrix


def directed_3cycles(names, matrix):
    cycles = []
    for triple in combinations(names, 3):
        for a, b, c in permutations(triple):
            if matrix[a][b] and matrix[b][c] and matrix[c][a]:
                cycles.append((a, b, c))
                break
    return cycles


def sccs(names, matrix):
    index = 0
    stack = []
    on_stack = set()
    indices = {}
    low = {}
    out = []

    def visit(v):
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in names:
            if w == v or not matrix[v][w]:
                continue
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(tuple(sorted(comp)))

    for name in names:
        if name not in indices:
            visit(name)
    return out


def hamiltonian_paths(names, matrix):
    count = 0
    examples = []
    for path in permutations(names):
        if all(matrix[path[i]][path[i + 1]] for i in range(len(path) - 1)):
            count += 1
            if len(examples) < 2:
                examples.append(path)
    return count, examples


def main():
    print("LRC tournament-computation speedup roadmap (S583)")
    print("Scores: certificate_safety, proof_payoff, algorithmic_leverage, -cost, maturity.")
    print()
    ordered = sorted(ENGINES, key=lambda e: (e.priority(), e.key()), reverse=True)
    print("Priority order")
    for rank, engine in enumerate(ordered, 1):
        print(
            f"{rank:2d}. {engine.name:<36} priority={engine.priority():2d} "
            f"key={engine.key()} source={engine.source}"
        )
        print(f"    preserves: {engine.preserves}")
        print(f"    destroys:  {engine.destroys}")
        print(f"    next:      {engine.next_move}")

    names, matrix = build_tournament()
    scores = {name: sum(matrix[name].values()) for name in names}
    cycles = directed_3cycles(names, matrix)
    hp_count, hp_examples = hamiltonian_paths(names, matrix)
    print()
    print("Tournament Analysis fingerprint on speedup engines")
    print("  vertices: speedup engines")
    print("  observable: (certificate_safety, proof_payoff, leverage, -cost, maturity)")
    print("  switch: lexicographically larger observable beats smaller")
    print("  tie Hamiltonian path: declaration order")
    print(f"  score histogram: {dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed 3-cycles: {len(cycles)}")
    print(f"  SCCs: {sccs(names, matrix)}")
    print(f"  Hamiltonian path count: {hp_count}")
    print(f"  Hamiltonian path examples: {hp_examples}")

    print()
    print("Recommended implementation stack")
    print("  1. exact source-marked reachability target (THM-381)")
    print("  2. certificate-gate DAG with Cprime/Phi/P(S) terminal calls (THM-398/HYP-2112/HYP-2108)")
    print("  3. A^2 cache for small quotient-class lookup, canonical fallback outside range")
    print("  4. SCC/good-cut/private-leaf peeling on endpoint-owner proof graphs")
    print("  5. L/M/R middle automata so only live middle cells receive CRT/metric detail")
    print("  6. threshold-colored section/fiber pre-sieves before exact interval search")


if __name__ == "__main__":
    main()
