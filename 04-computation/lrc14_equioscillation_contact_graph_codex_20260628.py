#!/usr/bin/env python3
"""HYP-3265 scout: LRC14 equioscillation contact graph.

This is a structural check, not a proof of LRC(14).  It verifies the exact
six-contact picture for the AP and Goddyn-Wong rows, then packages the
topology/geometry/graph proof routes suggested by that picture.

Tournament Analysis deliberately uses proof carriers as vertices, not runners.
The pairwise observable is which carrier preserves more of the LRC14
equioscillation packet while destroying fewer sidecar coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


MODULUS = 14
UNITS = tuple(a for a in range(1, MODULUS) if gcd(a, MODULUS) == 1)
UNIT_PAIRS = tuple((a, MODULUS - a) for a in UNITS if a < MODULUS - a)


def circle_distance_mod(residue: int, modulus: int = MODULUS) -> Fraction:
    r = residue % modulus
    return Fraction(min(r, modulus - r), modulus)


def safety_at_grid(row: tuple[int, ...], numerator: int, modulus: int = MODULUS) -> Fraction:
    return min(circle_distance_mod(s * numerator, modulus) for s in row)


def inverse_mod_unit(a: int, modulus: int = MODULUS) -> int:
    for x in range(1, modulus):
        if (a * x) % modulus == 1:
            return x
    raise ValueError(f"{a} is not a unit mod {modulus}")


def complement_pair(x: int, modulus: int = MODULUS) -> tuple[int, int]:
    y = (-x) % modulus
    return tuple(sorted((x % modulus, y)))


def binders_at_unit(row: tuple[int, ...], unit: int) -> tuple[int, ...]:
    return tuple(sorted(s for s in row if (s * unit) % MODULUS in (1, MODULUS - 1)))


def contact_rows(row: tuple[int, ...]) -> list[dict[str, object]]:
    rows = []
    for a in UNITS:
        inv = inverse_mod_unit(a)
        expected_pair = complement_pair(inv)
        rows.append(
            {
                "unit": a,
                "t": Fraction(a, MODULUS),
                "safety": safety_at_grid(row, a),
                "binders": binders_at_unit(row, a),
                "expected_pair": expected_pair,
                "antipodal_unit": MODULUS - a,
            }
        )
    return rows


def unit_pair_contact_map() -> dict[tuple[int, int], tuple[int, int]]:
    mapping = {}
    for a, b in UNIT_PAIRS:
        mapping[(a, b)] = complement_pair(inverse_mod_unit(a))
    return mapping


OBLIGATION_WEIGHT = {
    "unit_contact_exactness": 7,
    "antipodal_pair_index": 6,
    "case_split_power": 7,
    "danger_nerve_hole": 5,
    "kolmogorov_local_min": 5,
    "bulk_core_gluing": 6,
    "covering_floor_route": 5,
    "compression_legality": 6,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    note: str
    priority: int


CARRIERS = [
    Carrier(
        "unit_contact_matching",
        frozenset(
            {
                "unit_contact_exactness",
                "antipodal_pair_index",
                "case_split_power",
                "compression_legality",
            }
        ),
        frozenset({"bulk_core_gluing"}),
        "six unit contacts collapse to three complement-pair binders",
        0,
    ),
    Carrier(
        "antipodal_three_pair_quotient",
        frozenset(
            {
                "antipodal_pair_index",
                "case_split_power",
                "danger_nerve_hole",
                "compression_legality",
            }
        ),
        frozenset({"covering_floor_route"}),
        "quotient by t -> 1-t leaves three active pair slots",
        1,
    ),
    Carrier(
        "d7_index_degree_packet",
        frozenset(
            {
                "unit_contact_exactness",
                "antipodal_pair_index",
                "kolmogorov_local_min",
                "compression_legality",
            }
        ),
        frozenset({"bulk_core_gluing"}),
        "reads the three pairs as the odd D7/Borsuk-Ulam index",
        2,
    ),
    Carrier(
        "danger_nerve_boundary_holes",
        frozenset(
            {
                "danger_nerve_hole",
                "bulk_core_gluing",
                "case_split_power",
                "compression_legality",
            }
        ),
        frozenset({"antipodal_pair_index"}),
        "turns the six contacts into zero-measure holes of the danger cover",
        3,
    ),
    Carrier(
        "kolmogorov_convex_hull",
        frozenset(
            {
                "unit_contact_exactness",
                "kolmogorov_local_min",
                "compression_legality",
            }
        ),
        frozenset({"covering_floor_route", "bulk_core_gluing"}),
        "local Chebyshev test: active gradients should leave no descent",
        4,
    ),
    Carrier(
        "covering_kill_switch",
        frozenset(
            {
                "case_split_power",
                "covering_floor_route",
                "bulk_core_gluing",
                "compression_legality",
            }
        ),
        frozenset({"unit_contact_exactness"}),
        "a multiple of 14 kills all unit contacts and forces off-unit proof",
        5,
    ),
    Carrier(
        "bulk_equidistribution_glue",
        frozenset(
            {
                "bulk_core_gluing",
                "covering_floor_route",
                "compression_legality",
            }
        ),
        frozenset({"unit_contact_exactness", "antipodal_pair_index"}),
        "handles positive-measure bulk once core contacts are classified",
        6,
    ),
    Carrier(
        "raw_safety_scalar",
        frozenset(),
        frozenset(OBLIGATION_WEIGHT),
        "diagnostic only; forgets contacts, owners, topology, and sign",
        7,
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


def scc_sizes(graph: dict[str, set[str]]) -> list[int]:
    names = list(graph)
    rev = {n: set() for n in names}
    for u, outs in graph.items():
        for v in outs:
            rev[v].add(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(u: str) -> None:
        seen.add(u)
        for v in graph[u]:
            if v not in seen:
                dfs(v)
        order.append(u)

    for n in names:
        if n not in seen:
            dfs(n)

    seen.clear()
    sizes = []

    def rdfs(u: str) -> int:
        seen.add(u)
        total = 1
        for v in rev[u]:
            if v not in seen:
                total += rdfs(v)
        return total

    for n in reversed(order):
        if n not in seen:
            sizes.append(rdfs(n))
    return sorted(sizes, reverse=True)


def hamiltonian_path() -> list[str]:
    return [c.name for c in sorted(CARRIERS, key=lambda c: (-carrier_score(c), len(c.destroys), c.priority))]


def fmt_fraction(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def main() -> None:
    rows = {
        "AP": tuple(range(1, 14)),
        "GW": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24),
    }

    print("HYP-3265 LRC14 equioscillation contact graph")
    print("status: exact structural scout, not an LRC14 proof")
    print()
    print(f"units mod 14: {UNITS}")
    print(f"antipodal unit pairs: {UNIT_PAIRS}")
    print(f"unit-pair -> binding complement-pair map: {unit_pair_contact_map()}")
    print()

    for name, row in rows.items():
        print(f"{name} row: {row}")
        contacts = contact_rows(row)
        all_touch = all(c["safety"] == Fraction(1, 14) for c in contacts)
        print(f"  all six unit contacts touch 1/14: {all_touch}")
        for c in contacts:
            print(
                "  "
                f"t={fmt_fraction(c['t'])} "
                f"safety={fmt_fraction(c['safety'])} "
                f"binders={c['binders']} "
                f"expected_pair={c['expected_pair']} "
                f"antipode={c['antipodal_unit']}/14"
            )
        print()

    print("elementary case split at unit times:")
    print("  if no speed is 0 mod 14, every unit a has f_S(a/14) >= 1/14")
    print("  if some speed is 0 mod 14, every unit a is killed by that runner")
    print("  tight means the six unit contacts are global maxima with no higher off-unit peak")
    print()

    graph = adjacency()
    score_hist: dict[int, int] = {}
    for c in CARRIERS:
        score_hist[carrier_score(c)] = score_hist.get(carrier_score(c), 0) + 1

    print("Tournament Analysis over proof carriers")
    print("  pairwise observable: retained equioscillation proof payload minus destroyed sidecars")
    print("  switch/gauge: orient toward higher weighted retained payload; ties use fewer destroyed coordinates")
    print(f"  vertices={len(CARRIERS)}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(graph)}")
    print(f"  scc_sizes={scc_sizes(graph)}")
    print(f"  hamiltonian_path={' -> '.join(hamiltonian_path())}")
    print()
    for c in sorted(CARRIERS, key=lambda x: (-carrier_score(x), x.priority)):
        print(
            f"  {c.name}: score={carrier_score(c)} "
            f"preserves={sorted(c.preserves)} destroys={sorted(c.destroys)}"
        )


if __name__ == "__main__":
    main()
