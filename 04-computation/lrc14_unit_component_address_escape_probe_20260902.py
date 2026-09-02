#!/usr/bin/env python3
"""Exact controls for the unit component-address escape lemma.

Let I=[L,R] be a positive component of the 1/14-safe set of C.  For
2 <= d <= 7, a d-unit tail w is active on the d lifts of x exactly when

    ||w*x|| < d/14.

Thus an unsafe row d*C union T with |T|=d would put I strictly inside one
open effective w-tooth of width d/(7w).  Reduced endpoint denominators
Q_L,Q_R force gaps at least 1/(w Q_L),1/(w Q_R).  This script audits the
arithmetic, the initial-segment family, a p=5 minority-anchor control, and
an exact equality boundary for the strict additive certificate.

All arithmetic is rational.  This is a certificate audit, not LRC(14).
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


DELTA = Q(1, 14)
CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(distance(speed * phase) for speed in speeds)


def walls(speeds: tuple[int, ...]) -> tuple[Q, ...]:
    result = {Q(0), Q(1)}
    for speed in speeds:
        for k in range(speed):
            result.add(Q(14 * k + 1, 14 * speed))
            result.add(Q(14 * k + 13, 14 * speed))
    return tuple(sorted(result))


def global_wall_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    components: list[tuple[Q, Q]] = []
    cuts = walls(speeds)
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left + right) / 2
        if clearance(speeds, midpoint) < DELTA:
            continue
        require(clearance(speeds, left) >= DELTA, ("left-safe", speeds, left))
        require(clearance(speeds, right) >= DELTA, ("right-safe", speeds, right))
        if components and components[-1][1] == left:
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return tuple(components)


def bad_sheets(d: int, speed: int, quotient_phase: Q) -> tuple[int, ...]:
    return tuple(
        sheet
        for sheet in range(d)
        if distance(Q(speed) * (quotient_phase + sheet) / d) < DELTA
    )


def is_active(d: int, speed: int, quotient_phase: Q) -> bool:
    return distance(speed * quotient_phase) < Q(d, 14)


def least_unit_at_least(value: int, d: int) -> int:
    while gcd(value, d) != 1:
        value += 1
    return value


def ceil_fraction(value: Q) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def strict_integer_cutoff(value: Q) -> int:
    return value.numerator // value.denominator + 1


def certificate_ratios(d: int, component: tuple[Q, Q]) -> tuple[Q, Q, Q]:
    left, right = component
    width = right - left
    raw = Q(d, 7) / width
    simple = (Q(d, 7) - Q(1, min(left.denominator, right.denominator))) / width
    additive = (
        Q(d, 7) - Q(1, left.denominator) - Q(1, right.denominator)
    ) / width
    return raw, simple, additive


def audit_literal_witness(
    d: int,
    core: tuple[int, ...],
    tails: tuple[int, ...],
    selected: int,
    quotient_phase: Q,
    expected_sheet: int,
    expected_phase: Q,
    expected_clearance: Q,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    require(len(tails) == d, ("tail-count", d, tails))
    require(len(set(tails)) == d, ("distinct-tails", d, tails))
    require(all(gcd(tail, d) == 1 for tail in tails), ("unit-tails", d, tails))
    require(selected in tails, ("selected-tail", d, selected))
    require(clearance(core, quotient_phase) >= DELTA, ("core-safe", d))

    masks = tuple((tail, bad_sheets(d, tail, quotient_phase)) for tail in tails)
    for tail, mask in masks:
        require(len(mask) <= 1, ("one-sheet-cap", d, tail, mask))
        require(bool(mask) == is_active(d, tail, quotient_phase), ("activity-iff", d, tail))
    require(not dict(masks)[selected], ("selected-inactive", d, selected))

    burned = {sheet for _, mask in masks for sheet in mask}
    require(expected_sheet not in burned, ("free-sheet", d, expected_sheet, burned))
    physical_phase = (quotient_phase + expected_sheet) / d
    physical_row = tuple(d * speed for speed in core) + tails
    literal_clearance = clearance(physical_row, physical_phase)
    require(physical_phase == expected_phase, ("physical-phase", d, physical_phase))
    require(literal_clearance == expected_clearance, ("literal-clearance", d))
    require(literal_clearance >= DELTA, ("literal-safe", d))
    return masks


def main() -> None:
    source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(source_tree)), "no-assert")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(source_tree)
        ),
        "no-float",
    )

    # Initial segments C={1,...,13-d}.  Their first component and both
    # endpoint denominators are closed-form. Address data lowers the least
    # selected-tail d-unit threshold for d=2,3,4,6. For d=6 this does not
    # improve the whole row: six distinct positive units force max(T)>=17.
    expected_address_units = (17, 17, 17, 16, 13, 12)
    expected_raw_units = (23, 20, 19, 16, 17, 12)
    quotient_phases = (Q(6, 77), Q(23, 280), Q(3, 34), Q(3, 32), Q(3, 26), Q(1, 8))
    expected_sheets = (1, 1, 1, 1, 3, 2)
    expected_phases = (Q(83, 154), Q(101, 280), Q(37, 136), Q(7, 32), Q(27, 52), Q(17, 56))
    expected_clearances = (Q(6, 77), Q(23, 280), Q(3, 34), Q(3, 32), Q(3, 26), Q(5, 56))
    initial_rows: list[str] = []
    raw_improvements: list[int] = []

    for offset, d in enumerate(range(2, 8)):
        m = 13 - d
        core = tuple(range(1, m + 1))
        component = (Q(1, 14), Q(13, 14 * m))
        components = global_wall_components(core)
        require(component in components, ("initial-component", d, component))
        left, right = component
        width = right - left
        require(width == Q(d, 14 * m), ("initial-width", d, width))
        require((left.denominator, right.denominator) == (14, 14 * m), ("initial-Q", d))
        require(all(endpoint.denominator % 14 == 0 for endpoint in component), ("Q-divisibility", d))

        raw_ratio, simple_ratio, additive_ratio = certificate_ratios(d, component)
        raw_cutoff = least_unit_at_least(ceil_fraction(raw_ratio), d)
        simple_cutoff = least_unit_at_least(ceil_fraction(simple_ratio), d)
        additive_cutoff = least_unit_at_least(strict_integer_cutoff(additive_ratio), d)
        require(raw_cutoff == expected_raw_units[offset], ("raw-cutoff", d, raw_cutoff))
        require(simple_cutoff == expected_address_units[offset], ("simple-cutoff", d, simple_cutoff))
        require(additive_cutoff == expected_address_units[offset], ("additive-cutoff", d, additive_cutoff))
        if simple_cutoff < raw_cutoff:
            raw_improvements.append(d)

        selected = expected_address_units[offset]
        tails = [selected]
        candidate = 1
        while len(tails) < d:
            if gcd(candidate, d) == 1 and candidate not in tails:
                tails.append(candidate)
            candidate += 1
        masks = audit_literal_witness(
            d,
            core,
            tuple(tails),
            selected,
            quotient_phases[offset],
            expected_sheets[offset],
            expected_phases[offset],
            expected_clearances[offset],
        )
        initial_rows.append(
            f"d={d}:m={m}:I={left},{right}:W={width}:raw_unit={raw_cutoff}:"
            f"address_unit={simple_cutoff}:selected={selected}:x={quotient_phases[offset]}:"
            f"free={expected_sheets[offset]}:t={expected_phases[offset]}:"
            f"clearance={expected_clearances[offset]}:masks={masks}"
        )

    require(tuple(raw_improvements) == (2, 3, 4, 6), ("strict-improvement-d", raw_improvements))

    # A p=5 equality-capacity row on the h=420 minority wall.  Exact
    # component geometry cuts the old scale-free selected-tail threshold
    # from 1512 to the least admissible unit 141.
    wall_core = (1, 3, 5, 7, 9, 11, 13, 168)
    wall_tails = (1, 3, 7, 9, 141)
    wall_component = (Q(281, 2352), Q(293, 2352))
    wall_components = global_wall_components(wall_core)
    require(len(wall_components) == 58, ("wall-component-count", len(wall_components)))
    require(wall_component in wall_components, "wall-component-address")
    require(wall_component[1] - wall_component[0] == Q(1, 196), "wall-component-width")
    require(
        tuple(speed for speed in wall_core if distance(speed * wall_component[0]) == DELTA) == (168,),
        "wall-left-owner",
    )
    require(
        tuple(speed for speed in wall_core if distance(speed * wall_component[1]) == DELTA) == (168,),
        "wall-right-owner",
    )
    wall_raw, wall_simple, wall_additive = certificate_ratios(5, wall_component)
    require(wall_raw == 140, ("wall-raw-ratio", wall_raw))
    require(wall_simple == Q(1679, 12), ("wall-simple-ratio", wall_simple))
    require(wall_additive == Q(839, 6), ("wall-additive-ratio", wall_additive))
    require(least_unit_at_least(ceil_fraction(wall_simple), 5) == 141, "wall-unit-cutoff")
    wall_masks = audit_literal_witness(
        5,
        wall_core,
        wall_tails,
        141,
        Q(81, 658),
        3,
        Q(411, 658),
        DELTA,
    )
    wall_row = tuple(5 * speed for speed in wall_core) + wall_tails
    require(len(wall_row) == 13 and len(set(wall_row)) == 13, "wall-row-typing")
    require(840 in wall_row, "minority-anchor-speed")
    require(all(any(speed % modulus == 0 for speed in wall_row) for modulus in range(2, 15)), "wall-denominator-gates")

    # The additive inequality really must be strict.  This d=3 component
    # lies strictly in one open effective tooth, while attaining equality in
    # the sum of the two quantized endpoint gaps.
    equality_core = (9,)
    equality_component = (Q(19, 42), Q(23, 42))
    equality_d = 3
    equality_w = 4
    require(equality_component in global_wall_components(equality_core), "equality-component")
    eq_left, eq_right = equality_component
    eq_width = eq_right - eq_left
    tooth = (Q(25, 56), Q(31, 56))
    require(tooth[0] < eq_left < eq_right < tooth[1], "strict-tooth-containment")
    require(eq_left - tooth[0] == Q(1, equality_w * eq_left.denominator), "left-gap")
    require(tooth[1] - eq_right == Q(1, equality_w * eq_right.denominator), "right-gap")
    require(
        equality_w * eq_width
        == Q(equality_d, 7) - Q(1, eq_left.denominator) - Q(1, eq_right.denominator),
        "additive-equality",
    )
    require(all(is_active(equality_d, equality_w, x) for x in equality_component), "endpoint-activity")

    print("LRC14 UNIT COMPONENT-ADDRESS ESCAPE EXACT AUDIT")
    print("STATUS=FINITE-EXACT_CERTIFICATE_AUDIT_NOT_LRC14")
    print("LEMMA=failure forces every selected d-unit tail to contain each positive core component in one open effective tooth")
    print("FAILURE_BOUNDS=wW<d/7-1/Q_L and wW<d/7-1/Q_R and wW<=d/7-1/Q_L-1/Q_R")
    print("INITIAL_SEGMENT_ROWS")
    for row in initial_rows:
        print(row)
    print("SELECTED_TAIL_THRESHOLD_STRICTLY_IMPROVED_AT=d=2,3,4,6")
    print("D6_ROW_LEVEL_CAVEAT=six_distinct_positive_6_units_force_max_tail_at_least_17")
    print(
        "MINORITY_WALL_P5="
        f"C={wall_core}:T={wall_tails}:components={len(wall_components)}:"
        f"I={wall_component[0]},{wall_component[1]}:W={wall_component[1]-wall_component[0]}:"
        f"raw_ratio={wall_raw}:simple_ratio={wall_simple}:additive_ratio={wall_additive}:"
        f"unit_cutoff=141:witness=411/658:clearance=1/14:masks={wall_masks}"
    )
    print(
        "STRICT_EQUALITY_CONTROL="
        f"d={equality_d}:C={equality_core}:I={eq_left},{eq_right}:w={equality_w}:"
        f"tooth={tooth[0]},{tooth[1]}:gaps={eq_left-tooth[0]},{tooth[1]-eq_right}:"
        f"wW={equality_w*eq_width}"
    )
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
