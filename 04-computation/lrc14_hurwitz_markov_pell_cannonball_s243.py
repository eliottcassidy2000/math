#!/usr/bin/env python3
"""Hurwitz/Markov/Pell/cannonball carrier scout for LRC14.

This is not an LRC proof.  It tests whether the user's arithmetic cues should
enter the current HYP-3072/HYP-3074 proof-carrier ledger as direct quotients or
as sidecar coordinates.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from math import isqrt, sqrt
from pathlib import Path


MAX_MARKOV = 10000
MAX_CANNONBALL_N = 100000
OUT = Path("05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out")


def markov_triples(max_coord: int) -> set[tuple[int, int, int]]:
    """Generate positive Markov triples x^2+y^2+z^2=3xyz by Vieta mutation."""
    seen: set[tuple[int, int, int]] = set()
    q: deque[tuple[int, int, int]] = deque([(1, 1, 1)])
    while q:
        triple = tuple(sorted(q.popleft()))
        if triple in seen or triple[-1] > max_coord:
            continue
        seen.add(triple)
        x, y, z = triple
        for nxt in (
            (3 * y * z - x, y, z),
            (x, 3 * x * z - y, z),
            (x, y, 3 * x * y - z),
        ):
            nxt = tuple(sorted(nxt))
            if min(nxt) > 0 and max(nxt) <= max_coord and nxt not in seen:
                q.append(nxt)
    return seen


def pell_numbers(limit: int) -> list[int]:
    vals = [0, 1]
    while vals[-1] <= limit:
        vals.append(2 * vals[-1] + vals[-2])
    return vals


def cannonball_solutions(max_n: int) -> list[tuple[int, int]]:
    """Solve 1^2+...+n^2=m^2 by search up to max_n."""
    hits = []
    for n in range(1, max_n + 1):
        s2 = n * (n + 1) * (2 * n + 1) // 6
        m = isqrt(s2)
        if m * m == s2:
            hits.append((n, m))
    return hits


def pell_wall_hits(limit: int) -> dict[str, list[tuple[int, int]]]:
    """Endpoint wall families from the triangular-tower archive."""
    hits: dict[str, list[tuple[int, int]]] = {
        "same_end_tail_prefix_m_m1_eq_2n_n1": [],
        "right_block_begins_left_m2_eq_2n2_2n_1": [],
        "left_block_begins_left_m2_eq_n_2n1": [],
    }
    for n in range(1, limit + 1):
        for m in range(1, 2 * limit + 1):
            if m * (m + 1) == 2 * n * (n + 1):
                hits["same_end_tail_prefix_m_m1_eq_2n_n1"].append((n, m))
            if m * m == 2 * n * n + 2 * n + 1:
                hits["right_block_begins_left_m2_eq_2n2_2n_1"].append((n, m))
            if m * m == n * (2 * n + 1):
                hits["left_block_begins_left_m2_eq_n_2n1"].append((n, m))
    return {k: v[:6] for k, v in hits.items()}


def lagrange_value(markov_number: int) -> float:
    """Markov/Lagrange spectrum value below 3: sqrt(9 - 4/m^2)."""
    return sqrt(9.0 - 4.0 / (markov_number * markov_number))


@dataclass(frozen=True)
class Carrier:
    name: str
    retained: frozenset[str]
    destroyed: frozenset[str]
    note: str


def carrier_bank() -> list[Carrier]:
    return [
        Carrier(
            "labelled_lrc_packet_ledger",
            frozenset(
                {
                    "anti_bohr_predicate",
                    "endpoint_owner",
                    "exact_scale",
                    "route",
                    "legal_exit",
                    "finite_certificate",
                    "controlled_forgetting",
                }
            ),
            frozenset(),
            "HYP-2963-style packet state; theorem currency.",
        ),
        Carrier(
            "route_state_closure_median",
            frozenset(
                {
                    "anti_bohr_predicate",
                    "route",
                    "legal_exit",
                    "sidecar_closure",
                    "finite_certificate",
                    "controlled_forgetting",
                }
            ),
            frozenset({"endpoint_owner_location"}),
            "S240/HYP-3074 interface after legal sidecar closure.",
        ),
        Carrier(
            "cross_carrier_portfolio",
            frozenset(
                {
                    "anti_bohr_predicate",
                    "exact_scale",
                    "route",
                    "sidecar_closure",
                    "formal_exit",
                    "controlled_forgetting",
                }
            ),
            frozenset({"packet_identity"}),
            "S238/HYP-3072 local portfolios glued by exits.",
        ),
        Carrier(
            "markov_three_leg_resonance",
            frozenset(
                {
                    "lagrange_markov_depth",
                    "continued_fraction_period",
                    "quadratic_unit",
                    "three_leg_residue",
                    "aggregate_warning",
                }
            ),
            frozenset({"endpoint_owner", "route", "anti_bohr_predicate"}),
            "HYP-2745/HYP-2753 three-leg language, not the LRC predicate.",
        ),
        Carrier(
            "hurwitz_threshold",
            frozenset(
                {
                    "approximation_exponent",
                    "badly_approximable_wall",
                    "continued_fraction_period",
                    "height_fence",
                }
            ),
            frozenset({"endpoint_owner", "route", "anti_bohr_predicate"}),
            "Best-approximant fence; LRC needs anti-approximation.",
        ),
        Carrier(
            "markov_pell_fixed_two_branch",
            frozenset(
                {
                    "lagrange_markov_depth",
                    "quadratic_unit",
                    "pell_cassini_gap",
                    "three_leg_residue",
                    "continued_fraction_period",
                }
            ),
            frozenset({"endpoint_owner", "route", "finite_certificate"}),
            "Fixed coordinate 2 turns Markov mutation into a Pell branch.",
        ),
        Carrier(
            "cannonball_square_pyramid_gate",
            frozenset(
                {
                    "square_pyramid_scalar",
                    "pell_cassini_gap",
                    "quadratic_unit",
                    "endpoint_wall",
                }
            ),
            frozenset({"route", "anti_bohr_predicate", "owner_support"}),
            "70^2 is visible residue; the Pell/elliptic address is the sidecar.",
        ),
        Carrier(
            "beatty_pell_endpoint_wall",
            frozenset(
                {
                    "endpoint_wall",
                    "quadratic_unit",
                    "carry_residue",
                    "shell_address",
                    "visible_token_warning",
                }
            ),
            frozenset({"route", "anti_bohr_predicate"}),
            "HYP-2456/T800: visible word is carry-decorated, not raw Sturmian.",
        ),
    ]


CRITICAL_AXES = {
    "anti_bohr_predicate",
    "endpoint_owner",
    "exact_scale",
    "route",
    "legal_exit",
    "finite_certificate",
    "sidecar_closure",
    "controlled_forgetting",
}


def carrier_score(c: Carrier) -> tuple[int, int, int, str]:
    critical = len(c.retained & CRITICAL_AXES)
    retained = len(c.retained)
    destroyed = len(c.destroyed)
    return (critical, retained, -destroyed, c.name)


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    scores = {c.name: carrier_score(c) for c in carriers}
    edges: dict[tuple[str, str], str] = {}
    for i, a in enumerate(carriers):
        for b in carriers[i + 1 :]:
            winner = a.name if scores[a.name] > scores[b.name] else b.name
            edges[(a.name, b.name)] = winner
    out_degree = Counter()
    for (a, b), winner in edges.items():
        out_degree[winner] += 1
        out_degree[a if winner == b else b] += 0
    ordered = sorted(carriers, key=lambda c: scores[c.name], reverse=True)
    directed_3cycles = 0
    for i, a in enumerate(carriers):
        for j, b in enumerate(carriers):
            for k, c in enumerate(carriers):
                if not (i < j < k):
                    continue
                names = [a.name, b.name, c.name]
                wins = set()
                for x, y in ((names[0], names[1]), (names[1], names[2]), (names[2], names[0])):
                    key = (x, y) if (x, y) in edges else (y, x)
                    wins.add(edges[key])
                if len(wins) == 3:
                    directed_3cycles += 1
    return {
        "score_order": [c.name for c in ordered],
        "score_hist": dict(Counter(out_degree.values())),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path_count": 1 if directed_3cycles == 0 else "not computed",
    }


def main() -> None:
    triples = markov_triples(MAX_MARKOV)
    markov_nums = sorted({x for tri in triples for x in tri})
    pell = pell_numbers(MAX_MARKOV)
    pell_set = set(pell)
    markov_pell = sorted(set(markov_nums) & pell_set)
    fixed_two = []
    for j in range(1, len(pell) // 2):
        a = pell[2 * j - 1]
        b = pell[2 * j + 1] if 2 * j + 1 < len(pell) else None
        if b is None or b > MAX_MARKOV:
            break
        tri = tuple(sorted((2, a, b)))
        if tri in triples:
            fixed_two.append(tri)

    cannon = cannonball_solutions(MAX_CANNONBALL_N)
    walls = pell_wall_hits(250)

    lines: list[str] = []
    emit = lines.append
    emit("LRC14 HURWITZ / MARKOV / PELL / CANNONBALL CARRIER SCOUT")
    emit("=" * 72)
    emit("")
    emit("Scope")
    emit(f"  Markov triples generated with max coordinate <= {MAX_MARKOV}: {len(triples)}")
    emit(f"  distinct Markov numbers: {len(markov_nums)}")
    emit(f"  cannonball scan: n <= {MAX_CANNONBALL_N}")
    emit("")
    emit("Hurwitz/Markov spectrum payload")
    emit("  Hurwitz constant: 1/sqrt(5), equality obstruction: golden ratio class")
    emit("  Markov refinement below 3: L_m = sqrt(9 - 4/m^2)")
    emit("  first Markov numbers:")
    emit("    " + ", ".join(map(str, markov_nums[:24])))
    emit("  first Markov-Pell intersection:")
    emit("    " + ", ".join(map(str, markov_pell[:12])))
    emit("  sample Lagrange values for Markov-Pell branch:")
    for m in markov_pell[:8]:
        if m == 0:
            continue
        emit(f"    m={m:4d}  L_m={lagrange_value(m):.12f}  approx_coeff={1/lagrange_value(m):.12f}")
    emit("")
    emit("Fixed-2 Markov-Pell branch")
    for tri in fixed_two[:8]:
        x, y, z = tri
        lhs = x * x + y * y + z * z
        rhs = 3 * x * y * z
        emit(f"  {tri}: lhs=rhs={lhs}, recurrence coordinate={z}")
    emit("  Readout: fixing one Markov coordinate at 2 turns Vieta mutation into")
    emit("  the odd-index Pell branch 1,5,29,169,985,... with recurrence a(n+1)=6a(n)-a(n-1).")
    emit("")
    emit("Cannonball gate")
    emit("  square-pyramidal solutions found:")
    for n, m in cannon:
        emit(f"    sum_{{i<= {n}}} i^2 = {m}^2")
    p = pell
    if 7 < len(p):
        emit(f"  Pell context: P5={p[5]}, P6={p[6]}, P7={p[7]}")
        emit(f"  Cassini gap: P5*P7 - P6^2 = {p[5] * p[7] - p[6] * p[6]}")
        emit(f"  Cannonball square: 70^2 = {70 * 70}, and 29*169 = {29 * 169}")
        emit(f"  n=24 side observations: 24 = P5-P3 = {p[5] - p[3]} = 2*P4 = {2 * p[4]}")
    emit("  Readout: the nontrivial cannonball scalar is an elliptic/global square")
    emit("  event whose visible square root is a Pell carry between two Markov-Pell")
    emit("  wall numbers.  Use it as a wall/carry sidecar, not as a direct quotient.")
    emit("")
    emit("Triangular-tower Pell endpoint walls")
    for key, vals in walls.items():
        emit(f"  {key}: {vals}")
    emit("  Readout: S712/T800 already give the LRC transfer grammar:")
    emit("  visible token = projection(shell address, quadratic carry, endpoint atom).")
    emit("")
    emit("Carrier ledger")
    for c in carrier_bank():
        emit(f"  {c.name}")
        emit(f"    retained={sorted(c.retained)}")
        emit(f"    destroyed={sorted(c.destroyed)}")
        emit(f"    note={c.note}")
    emit("")
    emit("Tournament Analysis")
    emit("  Candidate vertex sets considered: runners, speed gaps, fixed circle sections,")
    emit("  section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,")
    emit("  matroid circuits, proof obligations, Markov triples, Pell wall atoms,")
    emit("  cannonball square-pyramid events, and proof-carrier ledgers.")
    emit("  Chosen vertices: proof carriers and arithmetic sidecar types, not runners.")
    emit("  Pairwise observable: retained critical LRC axes, total retained payload,")
    emit("  and destroyed owner/route/anti-Bohr coordinates.")
    emit("  Switch/gauge: orient toward the carrier preserving more critical LRC axes;")
    emit("  tie by total payload and fewer destroyed coordinates.")
    fp = tournament_fingerprint(carrier_bank())
    for key, value in fp.items():
        emit(f"  {key}: {value}")
    emit("")
    emit("Proof-interface conclusion")
    emit("  Hurwitz/Markov tells us how best-approximation exceptions are stratified.")
    emit("  LRC14 asks for anti-approximation endpoint survival, so the direct scalar")
    emit("  theorem is facing the wrong way.  The transferable object is the sidecar:")
    emit("  continued-fraction period, Markov depth, quadratic unit, Pell carry, endpoint")
    emit("  shell address, and named destroyed coordinates.  In HYP-3072 language this")
    emit("  belongs to the blindness-pair/resonance-portfolio ledger.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
