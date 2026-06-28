#!/usr/bin/env python3
"""
HYP-3423 scout: q-uniform topology versus q-specific magnitude arithmetic.

This script does not attempt a new LRC14 proof.  It makes one proof guardrail
executable:

    C2/Borsuk-Ulam topology is q-uniform, so it can certify only the uniform
    residue half.  The Goddyn-Wong magnitude break is q-specific and must be
    routed through arithmetic, owner packets, or the covering-floor theorem.

The q-table uses the HYP-3413 empirical criterion as an input:

    canonical GW doubling (n-2)->2(n-2) is ON iff q == 1 mod 3.

Tournament Analysis vertices are proof routes/obligations, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def c6_orbits_mod_14() -> list[list[int]]:
    units = [a for a in range(1, 14) if gcd(a, 14) == 1]
    seen = set()
    orbits = []
    for x in range(1, 14):
        if x in seen:
            continue
        orbit = sorted({(u * x) % 14 for u in units})
        seen.update(orbit)
        orbits.append(orbit)
    return sorted(orbits, key=lambda o: (-len(o), o[0]))


@dataclass(frozen=True)
class QRow:
    q: int
    n: int
    bu_residue_charge: str
    q_mod_3: int
    gw_switch: str
    prime_gate_note: str
    requested_note: str


def q_rows(q_min: int = 3, q_max: int = 22) -> list[QRow]:
    rows = []
    for q in range(q_min, q_max + 1):
        gw_on = q % 3 == 1
        prime = is_prime(q)
        if prime and gw_on:
            prime_note = "prime: unit group has C3 / Eisenstein split"
        elif prime:
            prime_note = "prime: no C3 gate"
        elif gw_on:
            prime_note = "composite: same q mod 3 switch, not prime-field evidence"
        else:
            prime_note = "composite/off-row"

        if q in (4, 7):
            requested = "requested ON example"
        elif q in (5, 6):
            requested = "requested OFF example"
        elif q == 3:
            requested = "canonical doubling off; small-n exotic is separate"
        else:
            requested = ""

        rows.append(
            QRow(
                q=q,
                n=2 * q,
                bu_residue_charge="present",
                q_mod_3=q % 3,
                gw_switch="ON" if gw_on else "off",
                prime_gate_note=prime_note,
                requested_note=requested,
            )
        )
    return rows


@dataclass(frozen=True)
class Route:
    name: str
    critical_path: int
    closes_residue: int
    sees_q_specific_magnitude: int
    supplies_floor: int
    supplies_owner_certificate: int
    quotient_guardrail: int
    preserves_labels: int
    false_route: int
    note: str

    def score(self) -> int:
        return (
            18 * self.critical_path
            + 14 * self.sees_q_specific_magnitude
            + 13 * self.supplies_floor
            + 11 * self.supplies_owner_certificate
            + 9 * self.quotient_guardrail
            + 8 * self.closes_residue
            + 5 * self.preserves_labels
            - 35 * self.false_route
        )


def route_vertices() -> list[Route]:
    return [
        Route(
            "HYP3415_decorrelation_floor_Rprime",
            1,
            0,
            1,
            1,
            0,
            0,
            1,
            0,
            "critical-path covering theorem: arithmetic/floor, not topology",
        ),
        Route(
            "HYP3413_q_mod_3_GW_arithmetic_switch",
            0,
            0,
            1,
            0,
            0,
            0,
            1,
            0,
            "q-specific magnitude/census switch; ON iff q==1 mod 3",
        ),
        Route(
            "HYP3417_labelled_owner_current_packet",
            0,
            0,
            1,
            0,
            1,
            0,
            1,
            0,
            "local finite certificate: one even-cover label plus binding labels",
        ),
        Route(
            "HYP3416_recursive_quotient_guardrail",
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            0,
            "legal forgetting policy: mixed fibers must emit sidecars/debt",
        ),
        Route(
            "HYP3411_C6_two_orbits_fixed_apex_packet",
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            "structural split: units, evens, fixed 7",
        ),
        Route(
            "HYP3312_C2_BU_topological_charge",
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            "q-uniform residue charge; cannot select q=4,7 over q=5,6",
        ),
        Route(
            "real_cubic_C3_trace_equioscillation",
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            "cap/equioscillation half in Q(cos 2pi/7)",
        ),
        Route(
            "raw_topology_closes_magnitude_false_route",
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            1,
            "forbidden inference: uniform topological charge proves q-specific break",
        ),
    ]


def tournament_fingerprint(routes: list[Route]) -> dict[str, object]:
    priority = {r.name: i for i, r in enumerate(routes)}
    scores = {r.name: r.score() for r in routes}
    edges: dict[tuple[str, str], str] = {}

    def beats(a: Route, b: Route) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return priority[a.name] < priority[b.name]

    for a, b in combinations(routes, 2):
        winner, loser = (a, b) if beats(a, b) else (b, a)
        edges[(winner.name, loser.name)] = winner.name

    names = [r.name for r in routes]

    def has_edge(a: str, b: str) -> bool:
        return (a, b) in edges

    directed_3cycles = 0
    for triple in combinations(names, 3):
        out = []
        for x in triple:
            out.append(sum(1 for y in triple if x != y and has_edge(x, y)))
        if sorted(out) == [1, 1, 1]:
            directed_3cycles += 1

    adjacency = defaultdict(list)
    reverse = defaultdict(list)
    for a, b in edges:
        adjacency[a].append(b)
        reverse[b].append(a)

    seen = set()
    order = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adjacency[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    scc_sizes = []

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
            scc_sizes.append(len(comp))

    ham_count = 0
    first_path: tuple[str, ...] | None = None
    for perm in permutations(names):
        if all(has_edge(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            ham_count += 1
            if first_path is None:
                first_path = perm

    return {
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3cycles": directed_3cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_path_count": ham_count,
        "selected_path": first_path,
        "scores": scores,
    }


def main() -> None:
    print("HYP-3423: q-uniform topology versus q-specific magnitude arithmetic")
    print("=" * 76)
    print()
    print("C6 orbit decomposition of the 13 nonzero residues mod 14:")
    for orbit in c6_orbits_mod_14():
        if orbit == [7]:
            role = "fixed apex"
        elif all(gcd(x, 14) == 1 for x in orbit):
            role = "unit/binding orbit"
        else:
            role = "even/covering orbit"
        print(f"  {orbit}: {role}")
    print()
    print("Subfield/proof-half ledger:")
    print("  Q(cos 2pi/7): real cubic / C3 trace / cap-equioscillation residue half")
    print("  Qsqrt(-7): imaginary quadratic / C2-BU sign / Gauss-floor orientation")
    print("  Qsqrt(-3): Eisenstein gate / q mod 3 switch for GW magnitude doubling")
    print()

    rows = q_rows()
    print("q-table using the HYP-3413 criterion for canonical GW doubling:")
    print("  q  n   C2/BU residue charge  q%3  GW     arithmetic note")
    for row in rows:
        suffix = f"  [{row.requested_note}]" if row.requested_note else ""
        print(
            f"  {row.q:2d} {row.n:2d}  {row.bu_residue_charge:>7s}"
            f"{'':14s} {row.q_mod_3:>2d}   {row.gw_switch:>3s}    "
            f"{row.prime_gate_note}{suffix}"
        )
    print()

    total = len(rows)
    uniform_count = sum(1 for row in rows if row.bu_residue_charge == "present")
    gw_on_count = sum(1 for row in rows if row.gw_switch == "ON")
    print("Immediate guardrail readout:")
    print(f"  C2/BU residue charge present on {uniform_count}/{total} q-rows.")
    print(f"  GW magnitude switch ON on {gw_on_count}/{total} q-rows.")
    print("  Requested contrast: q=4,7 are ON; q=5,6 are off.")
    print("  Therefore a q-uniform topological flag cannot be the magnitude switch.")
    print()

    print("Quotient guardrail theorem target:")
    print("  A quotient may use C2/BU or C6 topology to close a residue/equioscillation")
    print("  obligation.  If it claims a q-specific magnitude conclusion, it must retain")
    print("  or restore at least one of:")
    print("    (1) the HYP-3413 q mod 3 / Eisenstein arithmetic switch;")
    print("    (2) the HYP-3417 labelled owner-current packet for local mixed fibers;")
    print("    (3) the S259/HYP-3418 two-adic covering-floor descent signal;")
    print("    (4) the HYP-3415 decorrelation floor / Rprime inequality.")
    print()

    routes = route_vertices()
    fp = tournament_fingerprint(routes)
    print("Tournament Analysis over proof-route vertices")
    print("  Pairwise observable: preserved proof coordinate plus forbidden forgetting debt.")
    print("  Switch/gauge: higher route score; ties use the declared priority path.")
    print("  Vertices and scores:")
    for route in routes:
        print(f"    {route.name}: score={route.score():3d}; {route.note}")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("  selected_path=")
    assert fp["selected_path"] is not None
    for i, name in enumerate(fp["selected_path"], 1):
        print(f"    {i}. {name}")
    print()

    print("Assumption challenge:")
    print("  Considered vertices: runners, residues, C6 orbits, subfields, q-rows,")
    print("  owner labels, floor packets, proof obligations, and quotient policies.")
    print("  Chosen vertices: proof-route obligations.  This preserves the LRC predicate")
    print("  each route can actually certify and records what it destroys before any")
    print("  scalar/topological compression is allowed.")


if __name__ == "__main__":
    main()
