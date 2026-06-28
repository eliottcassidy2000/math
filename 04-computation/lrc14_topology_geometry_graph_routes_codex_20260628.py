#!/usr/bin/env python3
"""HYP-3243 scout: topology, geometry, and graph proof-route atlas.

This is a structural scout, not a numerical proof of LRC(14).  It packages
the current topology/geometry/graph ideas as proof carriers and audits which
LRC payload each carrier preserves or destroys.

Tournament Analysis uses proof carriers as vertices:

* circle endpoint arrangements and wall crossings;
* oriented-matroid topes/cocircuits;
* Cech nerves of safe components;
* D7/Borsuk-Ulam index packets;
* cyclotomic witness strata Phi_14 / Phi_{14d};
* Green conductance graphs and normal-fan/Toeplitz faces;
* Lee-Yang roots, ear payload graphs, finite chamber atlases, state lifts.

The point is to keep visual proof routes honest: a pretty picture is useful
only when its quotient preserves the LRC predicate or names the lost sidecar.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


OBLIGATION_WEIGHT = {
    "base_phi14_witness": 6,
    "dilation_phi14d_witness": 5,
    "bulk_equidistribution": 6,
    "tight_locus_finiteness": 7,
    "boundary_cocircuit_classification": 6,
    "odd_negative_sign": 5,
    "green_trap_discharge": 4,
    "finite_formalizable_atlas": 5,
    "compression_legality": 6,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    visual: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    route: str
    guardrail: str
    priority: int


CARRIERS = [
    Carrier(
        "cyclotomic_witness_strata",
        "Phi_14 / Phi_{14d} witness packet",
        frozenset(
            {
                "base_phi14_witness",
                "dilation_phi14d_witness",
                "odd_negative_sign",
                "compression_legality",
            }
        ),
        frozenset({"bulk_equidistribution"}),
        "constructs the measure-zero core witnesses directly",
        "does not prove tight-locus finiteness or bulk slack by itself",
        0,
    ),
    Carrier(
        "oriented_matroid_topes_cocircuits",
        "endpoint arrangement topes, walls, and boundary cocircuits",
        frozenset(
            {
                "boundary_cocircuit_classification",
                "finite_formalizable_atlas",
                "compression_legality",
                "tight_locus_finiteness",
            }
        ),
        frozenset({"bulk_equidistribution"}),
        "turns open lonely intervals and equality atoms into sign cells",
        "must retain endpoint owners; raw cell counts are too coarse",
        1,
    ),
    Carrier(
        "circle_endpoint_arrangement",
        "danger-arc endpoints on R/Z",
        frozenset(
            {
                "base_phi14_witness",
                "boundary_cocircuit_classification",
                "finite_formalizable_atlas",
                "compression_legality",
            }
        ),
        frozenset({"odd_negative_sign"}),
        "draws the wall-crossing picture where safe topes are visible",
        "forgets D7 parity unless antipodal owners are recorded",
        2,
    ),
    Carrier(
        "d7_borsuk_ulam_index_packet",
        "free Z2 antipodal witness pairs on the heptagon",
        frozenset(
            {
                "base_phi14_witness",
                "odd_negative_sign",
                "boundary_cocircuit_classification",
                "compression_legality",
            }
        ),
        frozenset({"bulk_equidistribution"}),
        "turns HYP-3241 index 3 into a topological sign certificate",
        "needs a finite atlas to say which row carries the packet",
        3,
    ),
    Carrier(
        "cech_nerve_safe_components",
        "safe-component nerve and barcode",
        frozenset(
            {
                "bulk_equidistribution",
                "boundary_cocircuit_classification",
                "finite_formalizable_atlas",
                "compression_legality",
            }
        ),
        frozenset({"odd_negative_sign"}),
        "separates positive open mass from endpoint-only equality",
        "closed boundary atoms must not be counted as open components",
        4,
    ),
    Carrier(
        "normal_fan_toeplitz_fejer",
        "moment polytope normal fan and Toeplitz/Fejer faces",
        frozenset(
            {
                "green_trap_discharge",
                "finite_formalizable_atlas",
                "tight_locus_finiteness",
                "compression_legality",
            }
        ),
        frozenset({"base_phi14_witness"}),
        "builds separating hyperplanes and finite chamber inequalities",
        "positive semidefinite faces can overcompress the core witness",
        5,
    ),
    Carrier(
        "green_conductance_laplacian",
        "positive covariance conductance graph and Fiedler currents",
        frozenset(
            {
                "green_trap_discharge",
                "finite_formalizable_atlas",
                "compression_legality",
            }
        ),
        frozenset({"odd_negative_sign", "base_phi14_witness"}),
        "visualizes AP as the high-connectivity positive face",
        "clips negative covariance and odd Worpitzky payload",
        6,
    ),
    Carrier(
        "state_lift_forbidden_H7",
        "forbidden state-lift endpoint I(Omega,2)=7",
        frozenset(
            {
                "tight_locus_finiteness",
                "finite_formalizable_atlas",
                "odd_negative_sign",
                "compression_legality",
            }
        ),
        frozenset({"bulk_equidistribution", "base_phi14_witness"}),
        "turns a remaining bad atom into the named tournament contradiction",
        "requires a sharp packet-to-state-lift construction",
        7,
    ),
    Carrier(
        "lee_yang_pgf_root_motion",
        "danger-count PGF roots and zero-free regions",
        frozenset(
            {
                "bulk_equidistribution",
                "finite_formalizable_atlas",
                "compression_legality",
                "green_trap_discharge",
            }
        ),
        frozenset({"base_phi14_witness"}),
        "detects when roots hit the obstruction interval or stay confined",
        "unit-disk zero-freeness fails on overcrowded R-blocks",
        8,
    ),
    Carrier(
        "ear_payload_graph",
        "directed/odd/nested ear payload graph",
        frozenset(
            {
                "odd_negative_sign",
                "finite_formalizable_atlas",
                "compression_legality",
                "green_trap_discharge",
            }
        ),
        frozenset({"bulk_equidistribution"}),
        "tracks one-runner additions, parity ears, and root motion payloads",
        "raw ear counts forget endpoint-owner geometry",
        9,
    ),
    Carrier(
        "bulk_spec_equidistribution",
        "fixed-lattice SPEC / equidistribution floor",
        frozenset({"bulk_equidistribution", "compression_legality"}),
        frozenset({"base_phi14_witness", "boundary_cocircuit_classification"}),
        "proves positive measure away from the Vitali core",
        "measure information is blind exactly at the Phi core",
        10,
    ),
    Carrier(
        "raw_scalar_shadow",
        "lambda2, p0, q0+q6, or measure as a single number",
        frozenset(),
        frozenset(OBLIGATION_WEIGHT),
        "diagnostic only",
        "destroys the topology, owners, parity, and core witness sidecars",
        11,
    ),
]


def carrier_score(carrier: Carrier) -> int:
    keep = sum(OBLIGATION_WEIGHT[k] for k in carrier.preserves)
    lose = sum(OBLIGATION_WEIGHT[k] for k in carrier.destroys)
    return 2 * keep - lose


def beats(left: Carrier, right: Carrier) -> bool:
    left_key = (carrier_score(left), -len(left.destroys), -left.priority)
    right_key = (carrier_score(right), -len(right.destroys), -right.priority)
    return left_key > right_key


def adjacency() -> dict[str, set[str]]:
    graph = {c.name: set() for c in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
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


def sccs(graph: dict[str, set[str]]) -> list[list[str]]:
    names = list(graph)
    reverse = {name: set() for name in names}
    for src, dsts in graph.items():
        for dst in dsts:
            reverse[dst].add(src)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
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


def hamiltonian_path() -> list[Carrier]:
    return sorted(
        CARRIERS,
        key=lambda c: (carrier_score(c), -len(c.destroys), -c.priority),
        reverse=True,
    )


def score_hist() -> dict[int, int]:
    hist: dict[int, int] = {}
    for c in CARRIERS:
        score = carrier_score(c)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def main() -> None:
    graph = adjacency()
    path = hamiltonian_path()

    print("HYP-3243 LRC14 topology/geometry/graph proof-route atlas")
    print("status: structural scout, not a proof")
    print()

    print("open frontier after incoming HYP-3240/HYP-3241:")
    print("  witness construction: base Phi_14 + dilation Phi_{14d} mostly organized")
    print("  remaining hard work: tight-locus finiteness + bulk equidistribution")
    print()

    print("obligation weights:")
    for key, value in OBLIGATION_WEIGHT.items():
        print(f"  {key}: {value}")
    print()

    print("carrier scores:")
    for carrier in path:
        print(
            f"  {carrier.name}: score={carrier_score(carrier)} "
            f"preserves={sorted(carrier.preserves)} destroys={sorted(carrier.destroys)}"
        )
    print()

    print("Tournament Analysis:")
    print(f"  vertices={len(CARRIERS)}")
    print(f"  score_hist={score_hist()}")
    print(f"  directed_3cycles={directed_3cycles(graph)}")
    print(f"  scc_sizes={[len(c) for c in sccs(graph)]}")
    print("  hamiltonian_path=" + " -> ".join(c.name for c in path))
    print()

    print("route packets:")
    print(
        "  visual_topology: circle_endpoint_arrangement -> "
        "oriented_matroid_topes_cocircuits -> cech_nerve_safe_components -> "
        "finite chamber proof or state_lift_forbidden_H7"
    )
    print(
        "  witness_core: cyclotomic_witness_strata -> "
        "d7_borsuk_ulam_index_packet -> oriented cocircuit classification"
    )
    print(
        "  graph_energy: green_conductance_laplacian -> "
        "normal_fan_toeplitz_fejer -> trap discharge, with odd-negative sidecar"
    )
    print(
        "  root_motion: lee_yang_pgf_root_motion -> ear_payload_graph -> "
        "normal fan / finite atlas"
    )
    print(
        "  bulk_route: bulk_spec_equidistribution -> cech_nerve_safe_components -> "
        "boundary cocircuit exclusion"
    )
    print()

    print("assumption challenge:")
    print(
        "  vertices considered: proof carriers, topes, Cech components, endpoint owners, "
        "witness strata, conductance cuts, normal fan chambers, PGF roots, ears, "
        "and state-lift obligations."
    )
    print(
        "  challenged assumption: runners or raw scalar metrics are not enough; "
        "they destroy owner, parity, cocircuit, and witness-construction payload."
    )
    print()

    print("mermaid:")
    print("flowchart LR")
    print('  Bulk["bulk SPEC / equidistribution"] --> Cech["safe-component Cech nerve"]')
    print('  Phi["Phi_14 / Phi_{14d} witnesses"] --> D7["D7 Borsuk-Ulam index"]')
    print('  D7 --> Cocircuit["oriented topes / boundary cocircuits"]')
    print('  Cech --> Cocircuit')
    print('  Cocircuit --> Atlas["finite chamber atlas"]')
    print('  Green["Green conductance graph"] --> Fan["Toeplitz / Fejer normal fan"]')
    print('  Fan --> Atlas')
    print('  LeeYang["Lee-Yang PGF roots"] --> Ears["ear payload graph"]')
    print('  Ears --> Atlas')
    print('  Atlas --> StateLift["H=7 state-lift contradiction or named debt"]')


if __name__ == "__main__":
    main()
