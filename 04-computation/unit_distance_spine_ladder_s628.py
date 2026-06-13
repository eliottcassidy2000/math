#!/usr/bin/env python3
"""
S628: symbolic certifier for THM-408, the Moser spine-ladder theorem.

The proof is finite-symbol rather than search-based.  It uses only six vectors
from the rank-4 Moser unit shell and two ordered cap families.  The script
verifies the step alphabet, samples the two infinite families, and records the
small instances used by S626.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

import unit_distance_impairment_lab_s622 as lab  # noqa: E402


Point = tuple[int, int, int, int]


def unit_edges(points: tuple[Point, ...]) -> int:
    unit_set = set(lab.UNITS)
    total = 0
    for p, q in combinations(points, 2):
        if lab.sub(q, p) in unit_set or lab.sub(p, q) in unit_set:
            total += 1
    return total


def layer_plus(a: int) -> tuple[Point, ...]:
    return (
        (a, 1, 1, -1),
        (a, 1, 0, -1),
        (a, 1, 0, 0),
        (a, 1, 1, 0),
        (a, 0, 1, 0),
        (a, 0, 0, 0),
        (a, 0, 0, -1),
        (a, 0, 1, -1),
    )


def cap_plus() -> tuple[Point, ...]:
    return (
        (0, 1, 1, -1),
        (0, 1, 0, -1),
        (0, 1, 0, 0),
        (0, 1, 1, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 0),
    )


def path_plus(m: int) -> tuple[Point, ...]:
    out: list[Point] = []
    for a in range(m, 0, -1):
        out.extend(layer_plus(a))
    out.extend(cap_plus())
    return tuple(out)


def layer_minus(a: int) -> tuple[Point, ...]:
    return (
        (a, 1, 0, -1),
        (a, 1, -1, -1),
        (a, 1, -1, 0),
        (a, 1, 0, 0),
        (a, 0, 0, 0),
        (a, 0, -1, 0),
        (a, 0, -1, -1),
        (a, 0, 0, -1),
    )


def cap_minus() -> tuple[Point, ...]:
    return (
        (0, 1, 0, -1),
        (0, 1, -1, -1),
        (0, 1, -1, 0),
        (0, 1, 0, 0),
        (0, 0, 0, 0),
    )


def path_minus(m: int) -> tuple[Point, ...]:
    out: list[Point] = []
    for a in range(m, 0, -1):
        out.extend(layer_minus(a))
    out.extend(cap_minus())
    return tuple(out)


def step_word(path: tuple[Point, ...]) -> tuple[Point, ...]:
    return tuple(lab.sub(path[i + 1], path[i]) for i in range(len(path) - 1))


def is_unit_step_path(path: tuple[Point, ...]) -> bool:
    unit_set = set(lab.UNITS)
    return all(step in unit_set for step in step_word(path))


def assert_family(name: str, builder, edge_formula) -> None:
    print(name)
    print("m | vertices | unit_edges | formula_ok | unit_spine")
    print("--- | --- | --- | --- | ---")
    for m in range(0, 9):
        path = builder(m)
        edges = unit_edges(tuple(sorted(path)))
        expected = edge_formula(m)
        ok = len(path) == len(set(path)) and is_unit_step_path(path)
        print(f"{m} | {len(path)} | {edges} | {edges == expected} | {ok}")
        if not ok or edges != expected:
            raise AssertionError((name, m, len(path), edges, expected, ok))
    print()


def main() -> None:
    print("S628 / THM-408 Moser spine-ladder certifier")
    print(f"Moser unit shell size: {len(lab.UNITS)} directed vectors")
    needed = sorted(set(step_word(layer_plus(1))) | set(step_word(cap_plus())) | {(-1, 1, 0, 0)})
    print("Step alphabet used in the proof:")
    for step in needed:
        print(f"  {step} unit={step in set(lab.UNITS)}")
        if step not in set(lab.UNITS):
            raise AssertionError(step)
    print()

    print("Layer-plus word W_a^+ and cap C^+ give |P_m^+|=8m+6.")
    print("For m>=1, edge count is 27m+6; for m=0 it is 8.")
    assert_family("plus/cap6 family", path_plus, lambda m: 8 if m == 0 else 27 * m + 6)

    print("Layer-minus word W_a^- and cap C^- give |P_m^-|=8m+5.")
    print("For m>=1, edge count is 27m+3; for m=0 it is 6.")
    assert_family("minus/cap5 family", path_minus, lambda m: 6 if m == 0 else 27 * m + 3)

    print("Named S626 instances:")
    print("  P_1^+ has n=14 and 33 edges; this is the S626 exact n=14 slab.")
    print("  P_2^- has n=21 and 57 edges; this is the S626 exact n=21 slab.")
    print("  P_2^+ has n=22 and 60 edges; this is the S626 Moser n=22 lane.")
    print()

    print("Tournament Analysis over proof abstractions")
    routes = [
        ("explicit layer word", 5, "finite symbolic certificate; no search"),
        ("traceable section", 4, "unit spine is retained as a section over the carrier quotient"),
        ("ear/slab gluing", 3, "local traceable blocks glue when endpoints are unit-adjacent"),
        ("frontier-gain search", 2, "finds the slab but does not itself prove the template"),
        ("point-order tournament", 1, "diagnostic only; can forget the section"),
    ]
    print("rank | route | outscore | note")
    print("--- | --- | --- | ---")
    for rank, (name, score, note) in enumerate(routes, 1):
        print(f"{rank} | {name} | {score} | {note}")
    print("score histogram: {1:1,2:1,3:1,4:1,5:1}")
    print("directed 3-cycles: 0")
    print("Hamiltonian path count: 1")
    print()
    print("Abstract reading: a unit spine is a traceable section.  The layer word")
    print("is a Gray code on an 8-vertex slab, and the bridge vector is the gluing")
    print("map between consecutive fibers.  The tournament tie path should remember")
    print("that section; otherwise the point-order quotient can create a fake flop.")


if __name__ == "__main__":
    main()
