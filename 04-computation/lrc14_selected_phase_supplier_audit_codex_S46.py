#!/usr/bin/env python3
"""Exact-cell audit of the three selected-phase LRC(14) suppliers (S46).

For ten harmonic quotient speeds ``a`` and three detuned speeds ``delta`` at
sheet modulus ``g``, a selected witness is a phase ``u`` in ``[0,1]`` such
that

* ``||a_i u|| >= 1/14`` for every harmonic coordinate; and
* for some branch ``c in {0,...,g-1}``, all three quantities
  ``||delta_j (u+c)/g||`` are at least ``1/14``.

Every predicate changes only at a rational endpoint.  The program therefore
tests every endpoint and one midpoint of every complementary cell, giving an
exact existence decision and exact positive-cell measure.  It searches the
corrected q-patterns ``(3,3,3)``, ``(2,4,4)``, and ``(2,2,q)``.

Tournament-analysis audit: the useful vertices are harmonic-safe cells and
branch obligations, with incidence "branch c clears this cell".  That is a
bipartite hypergraph, not an oriented pair relation.  Any arbitrary tournament
gauge destroys the existential intersection predicate, so no tournament
fingerprint is promoted to proof evidence here.  This explicitly challenges
the assumption that runners or residue rows should be tournament vertices.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


RADIUS = F(1, 14)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: F) -> int:
    return -floor_fraction(-x)


def distance_to_integer(x: F) -> F:
    r = x - floor_fraction(x)
    return min(r, 1 - r)


def harmonic_good(a: tuple[int, ...], u: F) -> bool:
    return all(distance_to_integer(F(abs(x)) * u) >= RADIUS for x in a)


def branch_good(delta: tuple[int, int, int], g: int, u: F) -> tuple[bool, int | None]:
    for c in range(g):
        if all(
            distance_to_integer(F(abs(d)) * (u + c) / g) >= RADIUS
            for d in delta
        ):
            return True, c
    return False, None


def bad_branches(delta: int, g: int, u: F) -> tuple[int, ...]:
    """Exact strict bad row used by the pair-tower equality-wall fixture."""
    return tuple(
        c
        for c in range(g)
        if distance_to_integer(F(abs(delta)) * (u + c) / g) < RADIUS
    )


def harmonic_endpoints(a: tuple[int, ...]) -> set[F]:
    endpoints = {F(0), F(1)}
    for value in a:
        d = abs(value)
        for n in range(-1, d + 2):
            for sign in (-1, 1):
                u = (F(n) + sign * RADIUS) / d
                if 0 <= u <= 1:
                    endpoints.add(u)
    return endpoints


def detuned_endpoints(delta: tuple[int, int, int], g: int) -> set[F]:
    endpoints: set[F] = set()
    for value in delta:
        d = abs(value)
        for c in range(g):
            low = F(d * c, g)
            high = F(d * (c + 1), g)
            for n in range(floor_fraction(low) - 1, ceil_fraction(high) + 2):
                for sign in (-1, 1):
                    u = F(g, d) * (F(n) + sign * RADIUS) - c
                    if 0 <= u <= 1:
                        endpoints.add(u)
    return endpoints


@dataclass(frozen=True)
class Decision:
    positive_measure: F
    passing_endpoints: int
    witness: F | None
    branch: int | None
    cells: int


def decide(a: tuple[int, ...], delta: tuple[int, int, int], g: int) -> Decision:
    endpoints = sorted(harmonic_endpoints(a) | detuned_endpoints(delta, g))
    positive_measure = F(0)
    witness = None
    branch = None
    passing_endpoints = 0

    for u in endpoints:
        good, c = branch_good(delta, g, u)
        if harmonic_good(a, u) and good:
            passing_endpoints += 1
            if witness is None:
                witness, branch = u, c

    for left, right in zip(endpoints, endpoints[1:]):
        if left == right:
            continue
        middle = (left + right) / 2
        good, c = branch_good(delta, g, middle)
        if harmonic_good(a, middle) and good:
            positive_measure += right - left
            if witness is None:
                witness, branch = middle, c

    return Decision(positive_measure, passing_endpoints, witness, branch, len(endpoints) - 1)


@dataclass(frozen=True)
class Pattern:
    name: str
    g: int
    delta: tuple[int, int, int]


def reduced_denominator(delta: int, g: int) -> int:
    return g // math.gcd(abs(delta), g)


def patterns(delta_max: int, g_max: int) -> tuple[Pattern, ...]:
    rows: list[Pattern] = []

    q3 = [d for d in range(1, delta_max + 1) if reduced_denominator(d, 3) == 3]
    rows.extend(Pattern("333", 3, ds) for ds in combinations(q3, 3))

    q2 = [d for d in range(1, delta_max + 1) if reduced_denominator(d, 4) == 2]
    q4 = [d for d in range(1, delta_max + 1) if reduced_denominator(d, 4) == 4]
    rows.extend(Pattern("244", 4, (d2, *pair)) for d2 in q2 for pair in combinations(q4, 2))

    for g in range(2, g_max + 1):
        q2g = [d for d in range(1, delta_max + 1) if reduced_denominator(d, g) == 2]
        nonmultiples = [
            d for d in range(1, delta_max + 1) if reduced_denominator(d, g) != 1
        ]
        for pair in combinations(q2g, 2):
            for dx in nonmultiples:
                rows.append(Pattern("22q", g, (pair[0], pair[1], dx)))

    # Exact predicate is invariant under permuting delta rows.  Deduplicate the
    # generated (2,2,q) overlaps and keep a stable order.
    unique = {(row.name, row.g, tuple(sorted(row.delta))): row for row in rows}
    return tuple(unique[key] for key in sorted(unique))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-harmonic", type=int, default=14)
    parser.add_argument("--delta-max", type=int, default=15)
    parser.add_argument("--g-max", type=int, default=10)
    parser.add_argument("--max-cases", type=int, default=25_000)
    parser.add_argument("--stop-on-counterexample", action="store_true")
    args = parser.parse_args()
    if args.max_harmonic < 10 or args.delta_max < 3 or args.g_max < 2:
        raise SystemExit("bounds too small")

    harmonic_bank = tuple(combinations(range(1, args.max_harmonic + 1), 10))
    pattern_bank = patterns(args.delta_max, args.g_max)
    counts = {"333": 0, "244": 0, "22q": 0}
    zeros: list[tuple[tuple[int, ...], Pattern, Decision]] = []
    minimum: tuple[F, tuple[int, ...], Pattern, Decision] | None = None
    digest = hashlib.sha256()
    total = 0

    for a in harmonic_bank:
        for pattern in pattern_bank:
            if total >= args.max_cases:
                break
            decision = decide(a, pattern.delta, pattern.g)
            total += 1
            counts[pattern.name] += 1
            payload = (
                f"{a};{pattern.name};{pattern.g};{pattern.delta};"
                f"{decision.positive_measure};{decision.passing_endpoints};"
                f"{decision.witness};{decision.branch};{decision.cells}\n"
            )
            digest.update(payload.encode())
            if decision.witness is None:
                zeros.append((a, pattern, decision))
                if args.stop_on_counterexample:
                    break
            if decision.positive_measure > 0:
                candidate = (decision.positive_measure, a, pattern, decision)
                if minimum is None or candidate[0] < minimum[0]:
                    minimum = candidate
        if total >= args.max_cases or (zeros and args.stop_on_counterexample):
            break

    print("LRC14 SELECTED-PHASE SUPPLIER EXACT-CELL AUDIT")
    print("=" * 76)
    print(
        f"harmonic_bank=C([1,{args.max_harmonic}],10)={len(harmonic_bank)} "
        f"pattern_bank={len(pattern_bank)} delta_max={args.delta_max} g_max={args.g_max}"
    )
    print(f"cases={total}/{len(harmonic_bank) * len(pattern_bank)} by_pattern={counts}")
    print(f"counterexamples={len(zeros)} manifest_sha256={digest.hexdigest()}")
    if minimum is not None:
        measure, a, pattern, decision = minimum
        print(
            "minimum_positive_cell_measure="
            f"{measure} at A={a}, pattern={pattern.name}, g={pattern.g}, "
            f"delta={pattern.delta}, witness={decision.witness}, branch={decision.branch}"
        )
    if zeros:
        print("FIRST COUNTEREXAMPLES:")
        for a, pattern, decision in zeros[:20]:
            print(
                f"  A={a}; pattern={pattern.name}; g={pattern.g}; delta={pattern.delta}; "
                f"cells={decision.cells}; passing_endpoints={decision.passing_endpoints}"
            )
        print("VERDICT: AT LEAST ONE CORRECTED SELECTED-WITNESS PROP IS FALSE ON THIS BANK")
    else:
        print("VERDICT: NO COUNTEREXAMPLE ON THE STATED FINITE BANK; NOT A UNIFORM PROOF")
    print("Tournament Analysis: not promoted; carrier is safe-cell/branch incidence hypergraph")
    print("kept: exact joint witness predicate; destroyed by tournament: cell-branch incidence")
    print("challenged vertices: proof obligations and activation cells, not runners/residue rows")
    gap_two_rows = tuple(
        (delta, bad_branches(delta, 8, F(1, 6)))
        for delta in (178, 142, 33, 39)
    )
    expected_gap_two_rows = (
        (178, (1, 5)),
        (142, (0, 4)),
        (33, (2, 3)),
        (39, (6, 7)),
    )
    if gap_two_rows != expected_gap_two_rows:
        raise AssertionError((gap_two_rows, expected_gap_two_rows))
    partition = sorted(c for _, row in gap_two_rows for c in row) == list(range(8))
    print(
        "gap2_equality_obstruction="
        f"G=8,u=1/6,rows={gap_two_rows},partition={partition}"
    )
    adversarial_a = (1, 2, 3, 5, 7, 8, 9, 11, 12, 13)
    adversarial_anchors = (
        ("244", (34, 13, 47), 4, F(13042, 423423), 28, F(1, 14), 0, 289),
        ("333", (4, 46, 50), 3, F(451681, 14504490), 28, None, None, None),
        ("22q", (39, 51, 2), 6, F(18103, 504504), 26, None, None, None),
    )
    for name, delta, g, measure, endpoints, witness, branch, cells in adversarial_anchors:
        decision = decide(adversarial_a, delta, g)
        if decision.positive_measure != measure or decision.passing_endpoints != endpoints:
            raise AssertionError((name, decision, measure, endpoints))
        if witness is not None and (decision.witness, decision.branch, decision.cells) != (
            witness,
            branch,
            cells,
        ):
            raise AssertionError((name, decision, witness, branch, cells))
        print(
            f"adversarial_anchor={name},A={adversarial_a},g={g},delta={delta},"
            f"measure={decision.positive_measure},passing_endpoints="
            f"{decision.passing_endpoints},witness={decision.witness},"
            f"branch={decision.branch},cells={decision.cells}"
        )
    print(f"source_sha256={hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
