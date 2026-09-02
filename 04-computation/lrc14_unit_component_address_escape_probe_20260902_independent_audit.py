#!/usr/bin/env python3
"""Independent interval-intersection audit of unit component-address escape.

Unlike the primary global-wall program, this script constructs each runner's
closed safe intervals and intersects interval unions successively.  It also
checks the d-lift activity equivalence over a finite rational universe and
audits every actual effective-tooth containment in the displayed core bank.
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


def runner_safe_intervals(speed: int) -> tuple[tuple[Q, Q], ...]:
    return tuple(
        (Q(14 * k + 1, 14 * speed), Q(14 * k + 13, 14 * speed))
        for k in range(speed)
    )


def intersect_unions(
    left_union: tuple[tuple[Q, Q], ...],
    right_union: tuple[tuple[Q, Q], ...],
) -> tuple[tuple[Q, Q], ...]:
    raw: list[tuple[Q, Q]] = []
    left_index = 0
    right_index = 0
    while left_index < len(left_union) and right_index < len(right_union):
        left = max(left_union[left_index][0], right_union[right_index][0])
        right = min(left_union[left_index][1], right_union[right_index][1])
        if left <= right:
            raw.append((left, right))
        if left_union[left_index][1] < right_union[right_index][1]:
            left_index += 1
        else:
            right_index += 1

    merged: list[tuple[Q, Q]] = []
    for left, right in raw:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def successive_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    current = ((Q(0), Q(1)),)
    for speed in speeds:
        current = intersect_unions(current, runner_safe_intervals(speed))
    return tuple((left, right) for left, right in current if left < right)


def bad_sheets(d: int, speed: int, quotient_phase: Q) -> tuple[int, ...]:
    return tuple(
        sheet
        for sheet in range(d)
        if distance(Q(speed) * (quotient_phase + sheet) / d) < DELTA
    )


def containing_tooth(d: int, speed: int, component: tuple[Q, Q]) -> tuple[Q, Q] | None:
    left, right = component
    for k in range(-1, speed + 2):
        tooth_left = Q(14 * k - d, 14 * speed)
        tooth_right = Q(14 * k + d, 14 * speed)
        if tooth_left < left < right < tooth_right:
            return tooth_left, tooth_right
    return None


def audit_containment(d: int, speed: int, component: tuple[Q, Q]) -> bool:
    tooth = containing_tooth(d, speed, component)
    if tooth is None:
        return False
    left, right = component
    tooth_left, tooth_right = tooth
    width = right - left
    require(left.denominator % 14 == 0, ("left-denominator", component))
    require(right.denominator % 14 == 0, ("right-denominator", component))
    require(left - tooth_left >= Q(1, speed * left.denominator), ("left-gap", d, speed, component))
    require(tooth_right - right >= Q(1, speed * right.denominator), ("right-gap", d, speed, component))
    require(speed * width < Q(d, 7) - Q(1, left.denominator), ("left-bound", d, speed, component))
    require(speed * width < Q(d, 7) - Q(1, right.denominator), ("right-bound", d, speed, component))
    require(
        speed * width <= Q(d, 7) - Q(1, left.denominator) - Q(1, right.denominator),
        ("additive-bound", d, speed, component),
    )
    return True


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

    # Direct finite audit of the lift fact used by the proof: a d-unit tail
    # owns at most one sheet, and owns one iff ||wx||<d/14.
    rational_universe = tuple(
        sorted({Q(numerator, denominator) for denominator in range(1, 43) for numerator in range(denominator + 1)})
    )
    activity_checks = 0
    for d in range(2, 8):
        for speed in range(1, 41):
            if gcd(d, speed) != 1:
                continue
            for phase in rational_universe:
                mask = bad_sheets(d, speed, phase)
                require(len(mask) <= 1, ("one-sheet", d, speed, phase, mask))
                require(bool(mask) == (distance(speed * phase) < Q(d, 14)), ("active-iff", d, speed, phase))
                activity_checks += 2

    initial_components: list[str] = []
    containment_count = 0
    expected_counts = (14, 20, 20, 20, 18, 12)
    for offset, d in enumerate(range(2, 8)):
        core = tuple(range(1, 14 - d))
        components = successive_components(core)
        expected = (Q(1, 14), Q(13, 14 * (13 - d)))
        require(len(components) == expected_counts[offset], ("component-count", d, len(components)))
        require(expected in components, ("first-component", d, expected))
        require(all(left.denominator % 14 == 0 and right.denominator % 14 == 0 for left, right in components), ("all-Q", d))
        for component in components:
            for speed in range(1, 81):
                if gcd(d, speed) == 1 and audit_containment(d, speed, component):
                    containment_count += 1
        initial_components.append(f"d={d}:count={len(components)}:first={expected[0]},{expected[1]}")

    wall_core = (1, 3, 5, 7, 9, 11, 13, 168)
    wall_component = (Q(281, 2352), Q(293, 2352))
    wall_components = successive_components(wall_core)
    require(len(wall_components) == 58, ("wall-count", len(wall_components)))
    require(wall_component in wall_components, "wall-address")
    require(all(left.denominator % 14 == 0 and right.denominator % 14 == 0 for left, right in wall_components), "wall-Q")

    wall_tails = (1, 3, 7, 9, 141)
    quotient_phase = Q(81, 658)
    masks = tuple((tail, bad_sheets(5, tail, quotient_phase)) for tail in wall_tails)
    require(masks == ((1, (0,)), (3, ()), (7, (2,)), (9, (1,)), (141, ())), ("wall-masks", masks))
    require(3 not in {sheet for _, mask in masks for sheet in mask}, "wall-free-sheet")
    physical_phase = (quotient_phase + 3) / 5
    physical_row = tuple(5 * speed for speed in wall_core) + wall_tails
    require(physical_phase == Q(411, 658), ("wall-phase", physical_phase))
    require(min(distance(speed * physical_phase) for speed in physical_row) == DELTA, "wall-clearance")

    equality_d = 3
    equality_speed = 4
    equality_component = (Q(19, 42), Q(23, 42))
    require(equality_component in successive_components((9,)), "equality-address")
    equality_tooth = containing_tooth(equality_d, equality_speed, equality_component)
    require(equality_tooth == (Q(25, 56), Q(31, 56)), ("equality-tooth", equality_tooth))
    require(audit_containment(equality_d, equality_speed, equality_component), "equality-containment")
    left, right = equality_component
    require(
        equality_speed * (right - left)
        == Q(equality_d, 7) - Q(1, left.denominator) - Q(1, right.denominator),
        "strictness-equality",
    )

    print("LRC14 UNIT COMPONENT-ADDRESS ESCAPE INDEPENDENT AUDIT")
    print("STATUS=FINITE-EXACT_INDEPENDENT_INTERVAL_INTERSECTION_NOT_LRC14")
    print(f"ACTIVITY_CHECKS={activity_checks}")
    print("INITIAL_COMPONENTS=" + ";".join(initial_components))
    print(f"AUDITED_EFFECTIVE_TOOTH_CONTAINMENTS={containment_count}")
    print(f"MINORITY_WALL_COMPONENTS={len(wall_components)}:witness={physical_phase}:clearance={DELTA}")
    print(
        "STRICT_EQUALITY="
        f"d={equality_d}:w={equality_speed}:I={left},{right}:"
        f"tooth={equality_tooth[0]},{equality_tooth[1]}:wW={equality_speed*(right-left)}"
    )
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
