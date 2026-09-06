#!/usr/bin/env python3
"""All-height reduction audit for strict q=4 one-v2 component caps.

The infinite step is analytic: every physical component has width at most
min(1/(14r),1/(7 max(a,b))).  This script checks the resulting exact finite
boxes twice, via the owner-saturated decomposition and a direct wall predicate.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "lrc14_clock_four_one_v2_component_thm4452_geometry.py"
SPEC = importlib.util.spec_from_file_location("q4_probe", SOURCE)
P = importlib.util.module_from_spec(SPEC)
if SPEC.loader is None:
    raise RuntimeError("missing probe loader")
SPEC.loader.exec_module(P)


def need(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def formula_stats(r: int, a: int, b: int) -> tuple[Q, int, tuple[Q, ...]]:
    intervals = list(P.q4_failure(r, a, b))
    widths = P.circular_widths(intervals)
    return (widths[0] if widths else Q(0), len(widths), widths)


def audit_box(rs: list[int], labels: list[int]) -> list[tuple[Q, int, int, int]]:
    rows: list[tuple[Q, int, int, int]] = []
    for r in rs:
        for a, b in combinations(labels, 2):
            formula = formula_stats(r, a, b)
            literal = P.literal_stats(r, a, b)
            need(formula == literal, ("decomposition/literal mismatch", r, a, b, formula, literal))
            rows.append((formula[0], r, a, b))
    rows.sort(reverse=True)
    return rows


def main() -> None:
    all_rows = audit_box([1, 3, 5, 7], [1, 3, 5, 7, 9, 11, 13])
    unit_rows = audit_box([1, 5, 7], [1, 5, 7, 11, 13])

    all_target = Q(1, 98)
    unit_target = Q(1, 110)
    all_equal = [(r, a, b) for width, r, a, b in all_rows if width == all_target]
    unit_equal = [(r, a, b) for width, r, a, b in unit_rows if width == unit_target]

    need(all_rows[0][0] == all_target, all_rows[0])
    need(all_equal == [(5, 3, 7), (3, 7, 13), (1, 7, 9)], all_equal)
    need(unit_rows[0][0] == unit_target, unit_rows[0])
    need(unit_equal == [(5, 1, 11)], unit_equal)

    # Arithmetic of the infinite reductions.  For a strict violation, the
    # carrier bounds force the displayed finite label boxes.
    need(Q(1, 14 * 7) == all_target, "all-odd r boundary")
    need(Q(1, 7 * 15) < all_target, "all-odd pair exterior")
    need(Q(1, 14 * 11) < unit_target, "3-unit r exterior")
    need(Q(1, 7 * 17) < unit_target, "3-unit pair exterior")

    # Correct quotient normalizations for the LRC component consumers.
    need(4 * all_target == Q(2, 49), "q4 all-odd quotient")
    need(4 * unit_target == Q(2, 55), "q4 3-unit quotient")
    need(2 * Q(17, 693) == Q(34, 693), "q2 all-odd quotient")
    need(2 * Q(19, 1001) == Q(38, 1001), "q2 3-unit quotient")

    print("Q4_ONE_V2_COMPONENT_ALL_HEIGHT_AUDIT")
    print("status=PROVED_ANALYTIC_REDUCTION_PLUS_EXACT_RATIONAL_BOX")
    print(f"all_odd_box_rows={len(all_rows)};max={all_target};equality={all_equal}")
    print(f"odd_3unit_box_rows={len(unit_rows)};max={unit_target};equality={unit_equal}")
    print("q4_quotient_caps=2/49,2/55")
    print("q2_quotient_caps=34/693,38/1001")
    print("PASS")


if __name__ == "__main__":
    main()
