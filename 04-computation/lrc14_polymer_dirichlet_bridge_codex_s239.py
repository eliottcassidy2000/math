#!/usr/bin/env python3
"""S239 scout: renormalized signed polymers and Dirichlet sidecar energy.

This is a proof-interface scout, not a proof of LRC14.  It cross-reads older
Riesz-product / signed-polymer work, Poisson/theta finite-cell work, residual
capacitor cuts, and the S237 cycle-class observability matrix.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import sqrt


@dataclass(frozen=True)
class RelationPacket:
    name: str
    r6_box2: int
    corr: Fraction
    source: str
    route: str

    @property
    def corr_per_relation(self) -> Fraction:
        return self.corr / self.r6_box2

    @property
    def density(self) -> Fraction:
        return Fraction(self.r6_box2, 982)

    @property
    def abs_corr_float(self) -> float:
        return abs(float(self.corr))

    @property
    def coherence(self) -> float:
        return abs(float(self.corr)) / sqrt(self.r6_box2)


def f_decimal(x: Fraction, places: int = 6) -> str:
    return f"{float(x):.{places}f}"


def build_relation_packets() -> list[RelationPacket]:
    # HYP-2645 table.  Fractions use the displayed decimal precision only; the
    # scout uses them as typed evidence, not as exact theorem constants.
    return [
        RelationPacket(
            "AP_core",
            982,
            Fraction(302731, 1_000_000),
            "HYP-2645",
            "finite AP / bounded-core certificate",
        ),
        RelationPacket(
            "near_AP",
            924,
            Fraction(183854, 1_000_000),
            "HYP-2645",
            "signed wall transport / Freiman pocket",
        ),
        RelationPacket(
            "odd_AP",
            546,
            Fraction(213476, 1_000_000),
            "HYP-2645",
            "odd-sector owner/current sidecar",
        ),
        RelationPacket(
            "primes_ish",
            414,
            Fraction(9649, 1_000_000),
            "HYP-2645",
            "wide finite-cell / Poisson-theta route",
        ),
        RelationPacket(
            "squares_sidon",
            156,
            Fraction(536, 1_000_000),
            "HYP-2645",
            "wide/Sidon fast signed convergence",
        ),
    ]


@dataclass(frozen=True)
class Edge:
    u: str
    v: str
    conductance: Fraction


def laplacian(nodes: list[str], edges: list[Edge]) -> list[list[Fraction]]:
    idx = {n: i for i, n in enumerate(nodes)}
    L = [[Fraction(0) for _ in nodes] for _ in nodes]
    for e in edges:
        i = idx[e.u]
        j = idx[e.v]
        c = e.conductance
        L[i][i] += c
        L[j][j] += c
        L[i][j] -= c
        L[j][i] -= c
    return L


def solve_linear(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    n = len(b)
    M = [row[:] + [rhs] for row, rhs in zip(A, b)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("singular linear system")
        M[col], M[pivot] = M[pivot], M[col]
        div = M[col][col]
        M[col] = [x / div for x in M[col]]
        for row in range(n):
            if row == col:
                continue
            factor = M[row][col]
            if factor:
                M[row] = [x - factor * y for x, y in zip(M[row], M[col])]
    return [M[i][-1] for i in range(n)]


def effective_conductance(
    nodes: list[str], edges: list[Edge], source: str, sink: str
) -> Fraction:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        adjacency[e.u].add(e.v)
        adjacency[e.v].add(e.u)
    stack = [source]
    component = {source}
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if v not in component:
                component.add(v)
                stack.append(v)
    if sink not in component:
        return Fraction(0)

    internal = [n for n in nodes if n in component and n not in {source, sink}]
    idx = {n: i for i, n in enumerate(internal)}
    A = [[Fraction(0) for _ in internal] for _ in internal]
    b = [Fraction(0) for _ in internal]

    for e in edges:
        u, v, c = e.u, e.v, e.conductance
        for a, z in [(u, v), (v, u)]:
            if a in idx:
                ia = idx[a]
                A[ia][ia] += c
                if z in idx:
                    A[ia][idx[z]] -= c
                elif z == source:
                    b[ia] += c
                elif z == sink:
                    b[ia] += 0

    potentials = {source: Fraction(1), sink: Fraction(0)}
    potentials.update(dict(zip(internal, solve_linear(A, b))))
    current = Fraction(0)
    for e in edges:
        if e.u == source:
            current += e.conductance * (potentials[source] - potentials[e.v])
        elif e.v == source:
            current += e.conductance * (potentials[source] - potentials[e.u])
    return current


def min_cut_capacity(
    nodes: list[str], directed_edges: list[Edge], source: str, sink: str
) -> tuple[Fraction, tuple[str, ...]]:
    best = None
    middle = [n for n in nodes if n not in {source, sink}]
    for r in range(len(middle) + 1):
        for subset in combinations(middle, r):
            S = set(subset) | {source}
            if sink in S:
                continue
            cap = sum(e.conductance for e in directed_edges if e.u in S and e.v not in S)
            if best is None or cap < best[0]:
                best = (cap, tuple(sorted(S)))
    assert best is not None
    return best


def build_energy_networks() -> tuple[list[str], list[Edge], list[Edge]]:
    nodes = [
        "coarse_residual_charge",
        "raw_route_label",
        "arc_topology_tooth",
        "coarse_safe_stalk",
        "owner_strip_exit",
        "cycle_class_image",
        "known_cycle_span",
        "discharged",
        "phantom_f7_class",
    ]
    raw = [
        Edge("coarse_residual_charge", "raw_route_label", Fraction(1)),
        Edge("raw_route_label", "discharged", Fraction(1)),
    ]
    legal = [
        # S237 first-tooth split: 13 topology fibers, 2 coarse-stalk fibers.
        Edge("coarse_residual_charge", "arc_topology_tooth", Fraction(13)),
        Edge("arc_topology_tooth", "owner_strip_exit", Fraction(13)),
        Edge("coarse_residual_charge", "coarse_safe_stalk", Fraction(2)),
        Edge("coarse_safe_stalk", "owner_strip_exit", Fraction(2)),
        Edge("owner_strip_exit", "discharged", Fraction(15)),
        # HYP-3071 cycle-class template: 12 known generators, one named F7 atom.
        Edge("coarse_residual_charge", "cycle_class_image", Fraction(12)),
        Edge("cycle_class_image", "known_cycle_span", Fraction(12)),
        Edge("known_cycle_span", "discharged", Fraction(12)),
        Edge("cycle_class_image", "phantom_f7_class", Fraction(1)),
    ]
    return nodes, raw, legal


@dataclass(frozen=True)
class CarrierVertex:
    name: str
    vector: tuple[int, ...]


def build_tournament_vertices() -> list[CarrierVertex]:
    # Vector coordinates:
    # retained_payload, positive_test_measure, signed_cancellation,
    # boundary_energy, finite_cell_exactness, discharge_specificity,
    # cross_carrier_reuse, raw_scalar_penalty.
    return [
        CarrierVertex("renormalized_signed_polymer", (9, 8, 10, 5, 7, 8, 9, 0)),
        CarrierVertex("dirichlet_schur_certificate", (9, 7, 6, 10, 8, 9, 8, 0)),
        CarrierVertex("cross_carrier_pullback_portfolio", (8, 7, 7, 8, 8, 8, 10, 0)),
        CarrierVertex("cycle_class_observability", (8, 6, 6, 8, 9, 10, 7, 0)),
        CarrierVertex("riesz_positive_test_measure", (7, 10, 8, 4, 5, 6, 7, 0)),
        CarrierVertex("residual_capacitor_min_cut", (7, 5, 5, 9, 6, 8, 7, 0)),
        CarrierVertex("poisson_finite_cell", (6, 5, 8, 3, 10, 5, 6, 0)),
        CarrierVertex("repeated_residue_character", (6, 4, 9, 3, 7, 6, 6, 0)),
        CarrierVertex("absolute_mayer_shadow", (3, 2, 1, 1, 2, 2, 3, -8)),
        CarrierVertex("raw_route_scalar", (1, 1, 0, 0, 1, 1, 1, -10)),
    ]


def tournament_fingerprint(vertices: list[CarrierVertex]) -> dict[str, object]:
    n = len(vertices)
    wins = Counter()
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(n), 2):
        vi = vertices[i].vector
        vj = vertices[j].vector
        winner = i if vi > vj else j
        loser = j if winner == i else i
        wins[winner] += 1
        edges[(winner, loser)] = 1

    directed_3cycles = 0
    for a, b, c in combinations(range(n), 3):
        if (
            ((a, b) in edges and (b, c) in edges and (c, a) in edges)
            or ((a, c) in edges and (c, b) in edges and (b, a) in edges)
        ):
            directed_3cycles += 1

    order = sorted(range(n), key=lambda i: vertices[i].vector, reverse=True)
    hp_count = 1 if directed_3cycles == 0 else None
    return {
        "score_hist": dict(sorted(Counter(wins[i] for i in range(n)).items())),
        "directed_3cycles": directed_3cycles,
        "scc_sizes": [1] * n if directed_3cycles == 0 else "not computed",
        "hamiltonian_path_count": hp_count,
        "tie_path": [vertices[i].name for i in order],
    }


def print_relation_packet_section() -> None:
    packets = build_relation_packets()
    print("=== Signed relation-packet budget scout ===")
    print("source: HYP-2645 finite x-cell / Poisson-theta table")
    print(
        "packet                 R6_box2  corr       corr/R6     density_vs_AP  coherence   route"
    )
    for p in packets:
        print(
            f"{p.name:22s} {p.r6_box2:7d}  {f_decimal(p.corr):>8s}  "
            f"{f_decimal(p.corr_per_relation, 8):>10s}  "
            f"{f_decimal(p.density):>12s}  {p.coherence:9.6f}  {p.route}"
        )
    print()
    print("Interpretation:")
    print("- Raw relation density is not a theorem currency: odd_AP has fewer R6")
    print("  relations than near_AP but larger signed correction.")
    print("- Wide/Sidon packets have tiny signed activities and should be routed")
    print("  through finite-cell / Poisson-theta control, not an absolute Mayer sum.")
    print("- AP-like packets are the slow signed-convergence sector and need bounded")
    print("  core, Riesz-product, or cross-carrier pullback certificates.")
    print("- Therefore the polymer grammar should be typed by packet class and")
    print("  signed activity, not by raw support size or relation count.")
    print()


def print_dirichlet_section() -> None:
    nodes, raw, legal = build_energy_networks()
    raw_g = effective_conductance(nodes, raw, "coarse_residual_charge", "discharged")
    legal_g = effective_conductance(nodes, legal, "coarse_residual_charge", "discharged")
    raw_cut = min_cut_capacity(nodes, raw, "coarse_residual_charge", "discharged")
    legal_cut = min_cut_capacity(nodes, legal, "coarse_residual_charge", "discharged")
    f7_cut = min_cut_capacity(nodes, legal, "coarse_residual_charge", "phantom_f7_class")

    print("=== Dirichlet / Schur sidecar network ===")
    print("Toy conductances use S237 first-tooth counts and cycle-template rank.")
    print(f"raw_route_effective_conductance = {raw_g} = {float(raw_g):.6f}")
    print(f"legal_sidecar_effective_conductance = {legal_g} = {float(legal_g):.6f}")
    print(f"conductance_ratio_legal/raw = {legal_g / raw_g} = {float(legal_g / raw_g):.6f}")
    print(f"raw_min_cut = {raw_cut[0]} via source side {raw_cut[1]}")
    print(f"legal_min_cut_to_discharge = {legal_cut[0]} via source side {legal_cut[1]}")
    print(f"legal_min_cut_to_phantom_f7 = {f7_cut[0]} via source side {f7_cut[1]}")
    print()
    print("Interpretation:")
    print("- Collapsing to raw route labels leaves only a unit-capacity fake cut.")
    print("- Keeping topology/stalk/cycle sidecars Schur-complements to a positive")
    print(f"  energy certificate: the exact toy effective conductance is {legal_g}.")
    print("- The named phantom F7 coordinate remains a one-unit side exit, matching")
    print("  HYP-3071's 'outside the known cycle image' reading.")
    print("- Proof target: legal quotients must preserve this boundary energy or")
    print("  explicitly route the lost current to AP/GW, owner-strip, cycle span,")
    print("  or named THM-572/F7 debt.")
    print()


def print_cross_disciplinary_section() -> None:
    rows = [
        ("stat_mech", "signed polymer gas", "HYP-2540", "renormalize typed activities; reject absolute Mayer"),
        ("harmonic", "Poisson/theta finite cell", "HYP-2645", "use signed finite x-cell, not divergent dual envelope"),
        ("additive", "repeated-residue character kernel", "HYP-2632", "route saturated packets by mod-7 character channels"),
        ("electrical", "residual capacitor min-cut", "HYP-3037", "sidecars are first cuts, not scalar rankings"),
        ("Hodge", "cycle-class image", "HYP-3066/HYP-3071", "generated rows discharge; phantom F7 is a basis atom"),
        ("portfolio", "cross-carrier pullback", "HYP-3072", "score carrier rows by preserved predicate and destroyed coordinate"),
    ]
    print("=== Cross-disciplinary pullback table ===")
    for field, object_name, source, use in rows:
        print(f"{field:11s} | {object_name:34s} | {source:18s} | {use}")
    print()


def print_tournament_section() -> None:
    vertices = build_tournament_vertices()
    fp = tournament_fingerprint(vertices)
    print("=== Tournament Analysis ===")
    print("Vertex candidates considered:")
    print("  runners, route labels, packet rows, relation-lattice supports, polymers,")
    print("  residual currents, sidecar boundary conditions, cycle-class atoms,")
    print("  cross-carrier rows, and proof obligations.")
    print("Chosen vertices: proof carriers and renormalization/energy obligations.")
    print("Pairwise observable:")
    print("  retained payload, positive test measure, signed cancellation, boundary")
    print("  energy, finite-cell exactness, discharge specificity, cross-carrier reuse.")
    print("Gauge: orient toward the carrier with lexicographically stronger payload.")
    print(f"score_hist = {fp['score_hist']}")
    print(f"directed_3cycles = {fp['directed_3cycles']}")
    print(f"scc_sizes = {fp['scc_sizes']}")
    print(f"hamiltonian_path_count = {fp['hamiltonian_path_count']}")
    print("tie_hamiltonian_path:")
    for item in fp["tie_path"]:
        print(f"  {item}")
    print()


def main() -> None:
    print("S239 LRC14 RENORMALIZED POLYMER / DIRICHLET BRIDGE SCOUT")
    print("status: proof-interface synthesis; not a proof of LRC14")
    print("namespace: HYP-3073 / T1155 / LTI-220 / LTT-118")
    print("upstream rebase signal: HYP-3072 claimed cross-carrier pullback portfolio")
    print()
    print_cross_disciplinary_section()
    print_relation_packet_section()
    print_dirichlet_section()
    print_tournament_section()
    print("=== Assumption Challenge ===")
    print("Do not assume vertices are runners, routes, or median centers.")
    print("For this scout, the useful vertices are proof carriers, typed polymers,")
    print("and sidecar-energy obligations.  The preserved LRC predicate is: the")
    print("sidecar-completed residual cannot cover/discharge illegally.  The")
    print("destroyed data are individual runner labels, raw relation counts, and")
    print("raw route labels.  The quotient is legal only if signed activity is")
    print("renormalized by packet type or boundary current is preserved by a")
    print("Dirichlet/Schur certificate.")
    print()
    print("=== Next concrete theorem targets ===")
    print("1. Polymer target: define packet activities by (packet_type, signed_activity,")
    print("   finite_cell_route), then prove wide/Sidon and repeated-residue activities")
    print("   are summable after AP-like cores are removed or certified.")
    print("2. Dirichlet target: build the actual HYP-2963 residual sidecar graph and")
    print("   show every admissible Schur complement preserves positive conductance")
    print("   to named exits or isolates phantom F7 as a concrete boundary atom.")


if __name__ == "__main__":
    main()
