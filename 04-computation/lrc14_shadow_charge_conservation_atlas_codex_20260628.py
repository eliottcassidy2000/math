#!/usr/bin/env python3
"""HYP-3400 scout: LRC14 shadow-charge conservation atlas.

This is a synthesis computation, not a proof.  It turns the current web of
LRC14 proof routes into a finite "proof charge" ledger.  A shadow/quotient is
legal only when it preserves a witness-bearing payload, transfers that payload
to a dual shadow, or emits an explicit named debt.

The model deliberately treats proof reservoirs as Tournament Analysis vertices,
not runners, arcs, lags, graph nodes, or roots.  Those lower-level objects are
allowed to appear only as payload fields inside a reservoir.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


OBLIGATION_WEIGHT = {
    "open_hole_witness": 8,
    "core_phi14_witness": 8,
    "bulk_density_floor": 7,
    "finite_trap_discharge": 7,
    "tight_locus_finiteness": 9,
    "odd_negative_sign": 7,
    "quotient_descent_legality": 8,
    "formalizable_exit": 6,
    "scale_lift_legality": 6,
    "analytic_zero_control": 5,
    "endpoint_owner_memory": 6,
    "green_current_certificate": 5,
    "contact_holonomy_memory": 6,
}


@dataclass(frozen=True)
class Reservoir:
    name: str
    metaphor: str
    shadows: tuple[str, ...]
    preserves: frozenset[str]
    destroys: frozenset[str]
    transfers_to: tuple[str, ...]
    proof_exit: str
    compression_failure: str
    priority: int


RESERVOIRS = [
    Reservoir(
        name="index_theorem_degree_charge",
        metaphor="analytic Euler/Cech index equals topological Borsuk-Ulam degree",
        shadows=("analytic chi", "topological degree", "Gauss-sum index"),
        preserves=frozenset(
            {
                "open_hole_witness",
                "core_phi14_witness",
                "odd_negative_sign",
                "tight_locus_finiteness",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor"}),
        transfers_to=("cech_euler_hole_charge", "d7_borsuk_ulam_sign_charge"),
        proof_exit="nonzero index plus S-dependent floor, or named forcing-gap debt",
        compression_failure="ambient index describes the saddle but does not supply the S-dependent floor",
        priority=0,
    ),
    Reservoir(
        name="uniform_margin_floor_charge",
        metaphor="finite tight locus, contact-graph case split, and uniform margin away from AP/GW/dilations",
        shadows=("tight-locus census", "unit contact graph", "covering margin", "decorrelation floor"),
        preserves=frozenset(
            {
                "bulk_density_floor",
                "tight_locus_finiteness",
                "finite_trap_discharge",
                "scale_lift_legality",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"core_phi14_witness", "endpoint_owner_memory"}),
        transfers_to=("bulk_discrepancy_density_charge", "state_lift_obstruction_charge"),
        proof_exit="finite tight locus plus uniform positive margin for all other rows",
        compression_failure="random/census margin evidence still needs a rigorous uniform theorem",
        priority=1,
    ),
    Reservoir(
        name="contact_holonomy_curvature_charge",
        metaphor="shell-lag quotient curvature repaired by zeta7 contact holonomy",
        shadows=("shell-lag commutator", "zeta7 contact holonomy", "connection coordinate"),
        preserves=frozenset(
            {
                "core_phi14_witness",
                "finite_trap_discharge",
                "endpoint_owner_memory",
                "contact_holonomy_memory",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor", "analytic_zero_control"}),
        transfers_to=("cyclotomic_witness_address_charge", "tiling_lift_descent_charge"),
        proof_exit="zero curvature, holonomy-killed curvature, endpoint-cell lift, or named curvature debt",
        compression_failure="first Fourier holonomy is a connection coordinate, not a global terminal quotient",
        priority=2,
    ),
    Reservoir(
        name="cech_euler_hole_charge",
        metaphor="cap as measured Euler characteristic; lonely witness as cover hole",
        shadows=("danger-cover nerve", "Cech/Betti barcode", "inclusion-exclusion"),
        preserves=frozenset(
            {
                "open_hole_witness",
                "bulk_density_floor",
                "endpoint_owner_memory",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"odd_negative_sign"}),
        transfers_to=("d7_borsuk_ulam_sign_charge", "green_dirichlet_current_charge"),
        proof_exit="open safe component, cover-hole witness, or named closed-boundary debt",
        compression_failure="Euler characteristic can cancel holes unless Betti/owner sidecars survive",
        priority=3,
    ),
    Reservoir(
        name="d7_borsuk_ulam_sign_charge",
        metaphor="odd antipodal index on the heptagon reflection quotient",
        shadows=("D7 sign irrep", "i*sqrt(7) Gauss sign", "Borsuk-Ulam index"),
        preserves=frozenset(
            {
                "core_phi14_witness",
                "odd_negative_sign",
                "endpoint_owner_memory",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor"}),
        transfers_to=("cyclotomic_witness_address_charge", "cech_euler_hole_charge"),
        proof_exit="free-Z2 antipodal witness pair or odd-sign contradiction",
        compression_failure="positive/even quotients erase the sign-isotypic payload",
        priority=4,
    ),
    Reservoir(
        name="cyclotomic_witness_address_charge",
        metaphor="Phi_14 core witnesses and Phi_{14d} dilation witness grid",
        shadows=("Phi_14 units", "Phi_{14d} promoted grid", "de Moivre cubic"),
        preserves=frozenset(
            {
                "core_phi14_witness",
                "scale_lift_legality",
                "odd_negative_sign",
                "endpoint_owner_memory",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor", "tight_locus_finiteness"}),
        transfers_to=("lee_yang_root_motion_charge", "d7_borsuk_ulam_sign_charge"),
        proof_exit="explicit primitive witness or promoted dilation witness",
        compression_failure="single-field norm scalars miss the binding-row dip",
        priority=5,
    ),
    Reservoir(
        name="tiling_lift_descent_charge",
        metaphor="fixed-path tiling cover descends through half-tiling quotient with sidecars",
        shadows=(
            "path-presentation fiber",
            "parent automorphism word orbit",
            "rectangle/hourglass residue",
        ),
        preserves=frozenset(
            {
                "quotient_descent_legality",
                "endpoint_owner_memory",
                "finite_trap_discharge",
                "tight_locus_finiteness",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor"}),
        transfers_to=("normal_fan_toeplitz_slack_charge", "state_lift_obstruction_charge"),
        proof_exit="lifted witness, certified descent, or named rectangle/hourglass debt",
        compression_failure="orbit counts forget canary/filler fibers unless sidecars are attached",
        priority=6,
    ),
    Reservoir(
        name="autocorrelation_transport_charge",
        metaphor="AP Fejer autocorrelation; non-AP traps transport mass outward",
        shadows=("Fejer F7", "Gram/Bochner lag profile", "out-correlation residual"),
        preserves=frozenset(
            {
                "finite_trap_discharge",
                "green_current_certificate",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"core_phi14_witness", "endpoint_owner_memory"}),
        transfers_to=("green_dirichlet_current_charge", "normal_fan_toeplitz_slack_charge"),
        proof_exit="lag-transport inequality routed to trap class sidecar",
        compression_failure="lag profiles do not remember which endpoint owner caused the ripple",
        priority=7,
    ),
    Reservoir(
        name="green_dirichlet_current_charge",
        metaphor="positive covariance conductance, Fiedler cut, and Thomson current",
        shadows=("lambda2", "effective resistance", "Dirichlet energy"),
        preserves=frozenset(
            {
                "finite_trap_discharge",
                "green_current_certificate",
                "bulk_density_floor",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"odd_negative_sign", "core_phi14_witness"}),
        transfers_to=("autocorrelation_transport_charge", "cech_euler_hole_charge"),
        proof_exit="strict resistance excess or Fiedler bottleneck discharge",
        compression_failure="positive-part conductance clips negative covariance and odd payload",
        priority=8,
    ),
    Reservoir(
        name="normal_fan_toeplitz_slack_charge",
        metaphor="moment-cone face, Toeplitz lambda slack, shell dual L_y",
        shadows=("normal fan", "Toeplitz PSD margin", "k=8 shell magic dual"),
        preserves=frozenset(
            {
                "finite_trap_discharge",
                "tight_locus_finiteness",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"core_phi14_witness", "bulk_density_floor"}),
        transfers_to=("autocorrelation_transport_charge", "lee_yang_root_motion_charge"),
        proof_exit="separating hyperplane, shell contact, or finite chamber inequality",
        compression_failure="PSD and shell faces can overprice q3 and forget root witnesses",
        priority=9,
    ),
    Reservoir(
        name="lee_yang_root_motion_charge",
        metaphor="PGF root confinement, Joukowski/Hermite-Biehler, and off-circle phi4 dip",
        shadows=("Lee-Yang roots", "Joukowski fold", "Hermite-Biehler interlacing"),
        preserves=frozenset(
            {
                "analytic_zero_control",
                "finite_trap_discharge",
                "core_phi14_witness",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"endpoint_owner_memory"}),
        transfers_to=("normal_fan_toeplitz_slack_charge", "cyclotomic_witness_address_charge"),
        proof_exit="zero-free region, HB interlacing, or named root-collision debt",
        compression_failure="circle-root ideals fail on overcrowded tail blocks without SPEC debt",
        priority=10,
    ),
    Reservoir(
        name="bulk_discrepancy_density_charge",
        metaphor="Roth/Halasz/Erdos-Turan bulk floor plus Hensel/Krasner local lifts",
        shadows=("discrepancy height", "ET clocks", "p-adic unit-root lift"),
        preserves=frozenset(
            {
                "bulk_density_floor",
                "scale_lift_legality",
                "analytic_zero_control",
                "quotient_descent_legality",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"core_phi14_witness", "endpoint_owner_memory"}),
        transfers_to=("cyclotomic_witness_address_charge", "cech_euler_hole_charge"),
        proof_exit="positive measure floor, ET/Hensel split, or local-lift debt",
        compression_failure="density is blind at the Vitali core and to exact owners",
        priority=11,
    ),
    Reservoir(
        name="state_lift_obstruction_charge",
        metaphor="finite state lift to forbidden H=7 / residual named debt",
        shadows=("H=7 obstruction", "observer gluing", "finite-address packet"),
        preserves=frozenset(
            {
                "tight_locus_finiteness",
                "odd_negative_sign",
                "quotient_descent_legality",
                "endpoint_owner_memory",
                "formalizable_exit",
            }
        ),
        destroys=frozenset({"bulk_density_floor", "analytic_zero_control"}),
        transfers_to=("tiling_lift_descent_charge", "d7_borsuk_ulam_sign_charge"),
        proof_exit="residual atom lifts to forbidden state or becomes named debt",
        compression_failure="state labels prove nothing unless packet-to-state functor is faithful",
        priority=12,
    ),
    Reservoir(
        name="law_defect_entropy_charge",
        metaphor="failed laws emit conditional-entropy sidecars instead of being quotiented away",
        shadows=("commutativity defect", "associativity defect", "distributivity defect"),
        preserves=frozenset({"quotient_descent_legality", "formalizable_exit"}),
        destroys=frozenset(
            {
                "open_hole_witness",
                "core_phi14_witness",
                "bulk_density_floor",
                "finite_trap_discharge",
                "endpoint_owner_memory",
            }
        ),
        transfers_to=("tiling_lift_descent_charge", "state_lift_obstruction_charge"),
        proof_exit="typed sidecar ledger for every illegal compression",
        compression_failure="this is the failure detector, not the proof-carrying charge itself",
        priority=13,
    ),
    Reservoir(
        name="raw_scalar_shadow",
        metaphor="single number such as p0, lambda2, entropy, or root radius",
        shadows=("diagnostic scalar",),
        preserves=frozenset(),
        destroys=frozenset(OBLIGATION_WEIGHT),
        transfers_to=(),
        proof_exit="none without a legal carrier",
        compression_failure="collapses the proof packet below the LRC predicate",
        priority=14,
    ),
]


CONSERVATION_LAWS = [
    (
        "index_equality_law",
        "Analytic Cech/Euler charge and topological Borsuk-Ulam degree are "
        "two proposed shadows of the same index, but the S-dependent floor "
        "and forcing gap must be supplied separately.",
    ),
    (
        "no_naked_quotient",
        "A shadow is legal only if every destroyed payload is reconstructed, "
        "dual-annihilated, or named as debt.",
    ),
    (
        "curvature_holonomy_law",
        "A quotient with nonzero shell-lag commutator carries curvature; a legal "
        "compression either kills it by zeta7 contact holonomy, lifts to endpoint "
        "cells, or names curvature debt.",
    ),
    (
        "contact_graph_case_law",
        "The six unit contacts form a proof graph: surviving global contacts route "
        "to the AP/GW core, surviving non-global contacts give a strict open row, "
        "and killed contacts require a promoted covering/floor witness.",
    ),
    (
        "observability_morse_law",
        "A finite chamber scalarization is legal only after observability columns "
        "separate hidden coordinates and the vector Morse energy either descends "
        "or names the critical debt.",
    ),
    (
        "bulk_core_transfer",
        "Positive-measure discrepancy charge can vanish at the Vitali core only "
        "when Phi_14/Phi_{14d} witness charge turns on.",
    ),
    (
        "even_odd_positive_negative_balance",
        "Even/positive Fejer, covariance, and conductance shadows must carry or "
        "discharge the odd/negative D7 sign charge.",
    ),
    (
        "boundary_current_law",
        "Green current, Cech boundary, normal-fan slack, and autocorrelation "
        "transport are four shadows of the same finite-boundary discharge.",
    ),
    (
        "lift_descent_law",
        "A proof built in a witness-rich cover may descend through a quotient "
        "only with fiber, automorphism, coboundary, and deletion sidecars.",
    ),
    (
        "root_radius_law",
        "Circular Lee-Yang/Joukowski zeros explain low-defect structure, while "
        "off-circle phi4/root-collision debt must be priced explicitly.",
    ),
]


PROOF_EXITS = [
    (
        "index_nonvanishing",
        "analytic chi / topological degree / Gauss-sum index is nonzero and "
        "the forcing gap from degree to lonely point is filled by an S-dependent floor",
    ),
    (
        "open_tope_or_hole",
        "the danger-cover complement has a certified open component/hole carrying a lonely witness",
    ),
    (
        "core_phi14_or_phi14d",
        "the packet is AP/GW/core-tight or a covering-tight dilation with explicit promoted witness",
    ),
    (
        "bulk_density_floor",
        "ET/Roth/Halasz/SPEC/Hensel data gives positive mass away from the Vitali core",
    ),
    (
        "finite_trap_discharge",
        "Toeplitz, Green, autocorrelation, normal-fan, or root-motion sidecars strictly discharge the trap",
    ),
    (
        "contact_graph_case_split",
        "unit-contact graph survives globally, survives with strict off-unit peak, or is killed with promoted covering witness",
    ),
    (
        "observable_morse_descent",
        "observability matrix separates hidden packet coordinates and vector Morse energy descends or names critical debt",
    ),
    (
        "curvature_holonomy_discharge",
        "lag/residue residuals have zero curvature, holonomy-killed curvature, or endpoint-cell lift",
    ),
    (
        "legal_lift_descent",
        "tiling/observer/finite-address lift descends through all quotient sidecars",
    ),
    (
        "state_lift_contradiction",
        "the surviving residual atom maps faithfully to the forbidden H=7 state",
    ),
    (
        "named_residual_debt",
        "the route fails honestly and records the exact missing payload",
    ),
]


def reservoir_score(reservoir: Reservoir) -> int:
    keep = sum(OBLIGATION_WEIGHT[k] for k in reservoir.preserves)
    lose = sum(OBLIGATION_WEIGHT[k] for k in reservoir.destroys)
    transfer_bonus = len(reservoir.transfers_to)
    shadow_bonus = min(len(reservoir.shadows), 3)
    return 2 * keep + transfer_bonus + shadow_bonus - lose


def beats(left: Reservoir, right: Reservoir) -> bool:
    left_key = (
        reservoir_score(left),
        len(left.preserves),
        -len(left.destroys),
        -left.priority,
    )
    right_key = (
        reservoir_score(right),
        len(right.preserves),
        -len(right.destroys),
        -right.priority,
    )
    return left_key > right_key


def adjacency() -> dict[str, set[str]]:
    graph = {r.name: set() for r in RESERVOIRS}
    for a, b in combinations(RESERVOIRS, 2):
        if beats(a, b):
            graph[a.name].add(b.name)
        else:
            graph[b.name].add(a.name)
    return graph


def directed_3cycles(graph: dict[str, set[str]]) -> int:
    count = 0
    names = list(graph)
    for a, b, c in combinations(names, 3):
        if b in graph[a] and c in graph[b] and a in graph[c]:
            count += 1
        if c in graph[a] and b in graph[c] and a in graph[b]:
            count += 1
    return count


def strongly_connected_components(graph: dict[str, set[str]]) -> list[list[str]]:
    names = list(graph)
    reverse = {name: set() for name in names}
    for src, dsts in graph.items():
        for dst in dsts:
            reverse[dst].add(src)

    seen: set[str] = set()
    order: list[str] = []

    def dfs1(node: str) -> None:
        seen.add(node)
        for nxt in graph[node]:
            if nxt not in seen:
                dfs1(nxt)
        order.append(node)

    def dfs2(node: str, comp: list[str]) -> None:
        seen.add(node)
        comp.append(node)
        for nxt in reverse[node]:
            if nxt not in seen:
                dfs2(nxt, comp)

    for name in names:
        if name not in seen:
            dfs1(name)
    seen.clear()
    comps = []
    for name in reversed(order):
        if name not in seen:
            comp: list[str] = []
            dfs2(name, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(graph: dict[str, set[str]]) -> int:
    names = list(graph)
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            src = names[last]
            for dst in graph[src]:
                nxt = index[dst]
                if not mask & (1 << nxt):
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[-1])


def priority_path(graph: dict[str, set[str]]) -> list[str]:
    unused = set(graph)
    path: list[str] = []
    while unused:
        candidates = [
            r
            for r in RESERVOIRS
            if r.name in unused
            and all(other in graph[r.name] for other in unused if other != r.name)
        ]
        if not candidates:
            candidates = [r for r in RESERVOIRS if r.name in unused]
        chosen = min(candidates, key=lambda r: r.priority)
        path.append(chosen.name)
        unused.remove(chosen.name)
    return path


def coverage_table() -> list[tuple[str, list[str]]]:
    rows = []
    for obligation in OBLIGATION_WEIGHT:
        carriers = [
            r.name
            for r in sorted(RESERVOIRS, key=lambda item: (-reservoir_score(item), item.priority))
            if obligation in r.preserves
        ]
        rows.append((obligation, carriers[:4]))
    return rows


def main() -> None:
    graph = adjacency()
    scores = {r.name: reservoir_score(r) for r in RESERVOIRS}
    hist = Counter(scores.values())
    sccs = strongly_connected_components(graph)
    path = priority_path(graph)

    print("HYP-3400 shadow-charge conservation atlas")
    print("=" * 78)
    print("claim=proof routes are shadows of conserved witness-bearing charge")
    print("status=synthesis scout / proof-route router; not an LRC14 proof")
    print(
        "principle=each quotient must preserve charge, transfer it to a dual "
        "shadow, or emit named debt"
    )
    print()

    print("CONSERVATION LAWS")
    for name, law in CONSERVATION_LAWS:
        print(f"  {name}: {law}")
    print()

    print("PROOF EXITS")
    for name, description in PROOF_EXITS:
        print(f"  {name}: {description}")
    print()

    print("RESERVOIR LEDGER")
    for r in sorted(RESERVOIRS, key=lambda item: (-scores[item.name], item.priority)):
        print(f"  {r.name} score={scores[r.name]}")
        print(f"    metaphor={r.metaphor}")
        print(f"    shadows={'; '.join(r.shadows)}")
        print(f"    preserves={','.join(sorted(r.preserves)) or 'none'}")
        print(f"    destroys={','.join(sorted(r.destroys)) or 'none'}")
        print(f"    transfers_to={','.join(r.transfers_to) or 'none'}")
        print(f"    exit={r.proof_exit}")
        print(f"    compression_failure={r.compression_failure}")
    print()

    print("OBLIGATION COVERAGE")
    for obligation, carriers in coverage_table():
        print(f"  {obligation}: {' | '.join(carriers) if carriers else 'UNCOVERED'}")
    print()

    print("PROPOSED FINITE THEOREM SCHEMA")
    print("  primitive LRC14 packet")
    print("    -> index_nonvanishing")
    print("    -> open_tope_or_hole")
    print("    or core_phi14_or_phi14d")
    print("    or bulk_density_floor")
    print("    or contact_graph_case_split")
    print("    or observable_morse_descent")
    print("    or finite_trap_discharge")
    print("    or curvature_holonomy_discharge")
    print("    or legal_lift_descent")
    print("    or state_lift_contradiction")
    print("    or named_residual_debt")
    print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof-charge reservoirs, not runners/arcs/lags/scalars")
    print(
        "  pairwise_observable=which reservoir preserves more weighted LRC "
        "obligations with fewer destroyed payloads"
    )
    print("  switch/gauge=A->B iff A carries more conserved charge after debt pricing")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(graph)}")
    print(f"  scc_sizes={[len(comp) for comp in sccs]}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(graph)}")
    print(f"  priority_path={' -> '.join(path)}")
    print()

    print("ASSUMPTION-CHALLENGE")
    print(
        "  alternate vertices considered: runners, arcs, endpoint walls, "
        "Fourier modes, autocorrelation lags, graph cuts, Cech holes, "
        "unit-contact graphs, observability columns, Morse chambers, state-lift "
        "atoms, tiling fibers, Hensel lift classes, contact-support holonomy "
        "classes, and proof obligations."
    )
    print(
        "  chosen vertices: reservoirs of proof charge; lower-level objects are "
        "payload fields inside reservoirs."
    )
    print(
        "  preserved predicate: LRC14 witness existence or a legal route to one "
        "of the finite proof exits."
    )
    print(
        "  destroyed information: endpoint owner, odd sign, bulk/core wall, "
        "root radius sidecar, quotient curvature, quotient fiber, and state-lift functor unless "
        "explicitly retained."
    )
    print(
        "  challenged assumption: a scalar extremality certificate is the proof. "
        "The atlas treats scalars as shadows of a higher-charge packet."
    )
    print()

    print("NEXT PROOF-FACING TESTS")
    print("  1. Turn every HYP-3202 non-AP trap into a charge-discharge row.")
    print("  2. Add a Cech/Green/autocorrelation boundary-current identity for traps.")
    print("  3. Attach ET/Hensel bulk charge to the Vitali-wall bulk/core splitter.")
    print("  4. Test index equality as descriptor plus the S-dependent floor as proof.")
    print("  5. Thread contact-holonomy curvature through every lag/residue quotient.")
    print("  6. Formalize no-naked-quotient as a reusable lemma over packet maps.")


if __name__ == "__main__":
    main()
