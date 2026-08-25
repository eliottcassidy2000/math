#!/usr/bin/env python3
"""Exact primary audit for THM-4136.

The proof has one infinite and one finite part.  The infinite part bounds every
component of the two-sheet bad quotient for a primitive odd tail ratio.  The
finite part checks the 68 coprime ratios below the sharp large-ratio gate
against three named clocks.  Selected full rows are then solved directly on
their complete rational wall grids, independently of the proof carrier.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

U = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
EVEN_BODY = tuple(2 * speed for speed in U)
DELTA = Q(1, 14)
J = (Q(33, 70), Q(27, 56))
J_LENGTH = J[1] - J[0]
UNIVERSAL_BETA = Q(2, 63)
LOW_RATIOS = ((1, 3), (1, 5), (3, 5), (1, 7), (3, 7), (5, 7))
CLOCKS = (Q(89, 1176), Q(181, 4704), Q(431, 4480))
DIRECT_ROWS = (
    (1, 3),
    (1, 13),
    (13, 25),
    (1, 27),
    (3, 9),
    (3, 27),
    (101, 103),
    (2001, 18009),
)
EXPECTED_SEMANTIC = "4e62dbd9b4cd87a7044687d48668449b7129ac11adb29512803ae7d01be61dbf"

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def fmt(value: Q) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(circle_distance(speed * phase) for speed in speeds)


def closed_safe(speeds: tuple[int, ...], phase: Q) -> bool:
    return clearance(speeds, phase) >= DELTA


def strict_danger_walls(speeds: tuple[int, ...]) -> tuple[Q, ...]:
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for integer in range(speed):
            for sign in (-1, 1):
                walls.add(((Q(integer) + sign * DELTA) / speed) % 1)
    return tuple(sorted(walls))


def positive_safe_components(
    speeds: tuple[int, ...],
) -> tuple[tuple[Q, ...], tuple[tuple[Q, Q], ...]]:
    walls = strict_danger_walls(speeds)
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(walls, walls[1:]):
        middle = (left + right) / 2
        if not closed_safe(speeds, middle):
            continue
        require(closed_safe(speeds, left), ("safe left endpoint", left, right))
        require(closed_safe(speeds, right), ("safe right endpoint", left, right))
        cells.append((left, right))
    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left and closed_safe(speeds, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return walls, tuple(components)


def body_interval_audit() -> dict[str, object]:
    walls, components = positive_safe_components(U)
    maximum = max(right - left for left, right in components)
    maximizers = tuple(component for component in components if component[1] - component[0] == maximum)
    require(len(walls) == 244, "body wall count")
    require(len(components) == 28, "body component count")
    require(maximum == J_LENGTH == Q(3, 280), "body maximum length")
    require(maximizers == (J, (Q(29, 56), Q(37, 70))), "body maximizers")
    left_owners = tuple(speed for speed in U if circle_distance(speed * J[0]) == DELTA)
    right_owners = tuple(speed for speed in U if circle_distance(speed * J[1]) == DELTA)
    require(left_owners == (15,), "left owner")
    require(right_owners == (4,), "right owner")
    return {
        "wall_count": len(walls),
        "component_count": len(components),
        "J": (pair(J[0]), pair(J[1]), pair(J_LENGTH)),
        "owners": (left_owners, right_owners),
    }


def speed_bad(speed: int, phase: Q) -> bool:
    return circle_distance(speed * phase) < DELTA


def both_sheets_bad(p: int, q: int, w: Q) -> bool:
    z = w / 2
    sheet_zero = speed_bad(p, z) or speed_bad(q, z)
    sheet_one = speed_bad(p, z + Q(1, 2)) or speed_bad(q, z + Q(1, 2))
    return sheet_zero and sheet_one


def two_sheet_walls(p: int, q: int) -> tuple[Q, ...]:
    walls = {Q(0), Q(1)}
    for sheet_shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    phase_wall = (Q(integer) + sign * DELTA) / speed
                    walls.add((2 * (phase_wall - sheet_shift)) % 1)
    return tuple(sorted(walls))


def two_sheet_components(p: int, q: int) -> tuple[tuple[Q, Q], ...]:
    require(0 < p < q and p % 2 == q % 2 == 1, ("odd ratio", p, q))
    walls = two_sheet_walls(p, q)
    cells = tuple(
        (left, right)
        for left, right in zip(walls, walls[1:])
        if both_sheets_bad(p, q, (left + right) / 2)
    )
    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left and both_sheets_bad(p, q, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    require(not both_sheets_bad(p, q, Q(0)), ("zero outside quotient", p, q))
    for left, right in components:
        require(not both_sheets_bad(p, q, left), ("open quotient left", p, q, left))
        require(not both_sheets_bad(p, q, right), ("open quotient right", p, q, right))
    return tuple(components)


def quotient_audit() -> dict[str, object]:
    expected_low = {
        (1, 3): (),
        (1, 5): (),
        (3, 5): ((Q(13, 35), Q(8, 21)), (Q(13, 21), Q(22, 35))),
        (1, 7): ((Q(6, 49), Q(1, 7)), (Q(6, 7), Q(43, 49))),
        (3, 7): ((Q(2, 7), Q(15, 49)), (Q(34, 49), Q(5, 7))),
        (5, 7): ((Q(20, 49), Q(3, 7)), (Q(4, 7), Q(29, 49))),
    }
    low_ledger = []
    for p, q in LOW_RATIOS:
        components = two_sheet_components(p, q)
        require(components == expected_low[(p, q)], ("low quotient", p, q))
        beta = max((right - left for left, right in components), default=Q(0))
        require(beta <= UNIVERSAL_BETA, ("low quotient bound", p, q, beta))
        low_ledger.append(
            (
                (p, q),
                tuple((pair(left), pair(right)) for left, right in components),
                pair(beta),
            )
        )

    scan_max = Q(0)
    scan_maximizers: list[tuple[int, int]] = []
    scanned = 0
    # The analytic proof covers every q >= 9.  This finite sweep is a hostile
    # control, not a dependency, so q <= 101 is already deliberately far past
    # the only residual gate q <= 25.
    for q in range(3, 102, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scanned += 1
            components = two_sheet_components(p, q)
            beta = max((right - left for left, right in components), default=Q(0))
            if q >= 9:
                require(beta <= Q(2, 7 * q), ("analytic interval-width control", p, q, beta))
            require(beta <= UNIVERSAL_BETA, ("universal quotient control", p, q, beta))
            if beta > scan_max:
                scan_max = beta
                scan_maximizers = [(p, q)]
            elif beta == scan_max:
                scan_maximizers.append((p, q))
    require(scan_max == UNIVERSAL_BETA, "scan maximum")
    require(scan_maximizers == [(1, 9)], "scan maximum ratio")
    return {
        "low": tuple(low_ledger),
        "large_ratio_rule": "q>=9: every cross-component lies in one q-danger interval of width 2/(7q)",
        "universal_beta": pair(UNIVERSAL_BETA),
        "scan": (scanned, 101, pair(scan_max), tuple(scan_maximizers)),
    }


def clock_index(p: int, q: int) -> int:
    require(0 < p < q <= 25 and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
            ("finite ratio", p, q))
    if 13 not in (p, q):
        return 0
    if (p, q) in ((1, 13), (13, 25)):
        return 2
    return 1


def finite_clock_audit() -> dict[str, object]:
    pairs = tuple(
        (p, q)
        for q in range(3, 26, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1
    )
    require(len(pairs) == 68, "finite ratio count")
    expected_body_gaps = (Q(9, 98), Q(15, 196), Q(17, 224))
    category_pairs: list[list[tuple[int, int]]] = [[], [], []]
    category_minima = [Q(1), Q(1), Q(1)]
    for p, q in pairs:
        index = clock_index(p, q)
        phase = CLOCKS[index]
        body_gap = clearance(EVEN_BODY, phase)
        row_gap = clearance(EVEN_BODY + (p, q), phase)
        require(body_gap == expected_body_gaps[index], ("clock body gap", index, p, q))
        require(row_gap > DELTA, ("clock closes ratio", p, q, phase, row_gap))
        category_pairs[index].append((p, q))
        category_minima[index] = min(category_minima[index], row_gap)
    require(tuple(len(group) for group in category_pairs) == (56, 10, 2), "clock category counts")
    require(tuple(category_minima) == (Q(89, 1176), Q(15, 196), Q(17, 224)),
            "clock category minima")
    require(tuple(category_pairs[2]) == ((1, 13), (13, 25)), "exception clock pairs")
    return {
        "pair_count": len(pairs),
        "clocks": tuple(pair(phase) for phase in CLOCKS),
        "category_counts": tuple(len(group) for group in category_pairs),
        "category_minima": tuple(pair(value) for value in category_minima),
        "exception_pairs": tuple(category_pairs[2]),
        "pair_digest": digest(tuple(pairs)),
    }


def direct_full_row_control(a: int, b: int) -> tuple[Q, Q, int]:
    require(0 < a < b and a % 2 == b % 2 == 1, ("direct row tails", a, b))
    speeds = EVEN_BODY + (a, b)
    require(len(speeds) == len(set(speeds)) == 13, ("direct row distinct", a, b))
    require(reduce(gcd, speeds) == 1, ("direct row primitive", a, b))
    walls = strict_danger_walls(speeds)
    for left, right in zip(walls, walls[1:]):
        phase = (left + right) / 2
        gap = clearance(speeds, phase)
        if gap > DELTA:
            return phase, gap, len(walls)
    raise RuntimeError(f"no direct full-row witness for {(a, b)}")


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    require(len(U) == len(set(U)) == 11, "body universe")
    require(all(speed % 2 == 0 for speed in EVEN_BODY), "even body")
    require(reduce(gcd, EVEN_BODY) == 2, "even body gcd")

    body = body_interval_audit()
    quotient = quotient_audit()
    require(3 * J_LENGTH - UNIVERSAL_BETA == Q(1, 2520) > 0, "gcd scale gate")
    require(J_LENGTH - Q(2, 7 * 27) == Q(1, 7560) > 0, "primitive ratio gate")
    finite = finite_clock_audit()

    hostile_w = Q(1, 9)
    hostile_z = hostile_w / 2
    hostile_sheet_zero = (circle_distance(hostile_z), circle_distance(9 * hostile_z))
    hostile_sheet_one = (
        circle_distance(hostile_z + Q(1, 2)),
        circle_distance(9 * (hostile_z + Q(1, 2))),
    )
    require(hostile_sheet_zero == (Q(1, 18), Q(1, 2)), "fixed-sheet hostile zero")
    require(hostile_sheet_one == (Q(4, 9), Q(0)), "fixed-sheet hostile one")
    require(both_sheets_bad(1, 9, hostile_w), "hostile quotient membership")

    direct_rows = []
    for a, b in DIRECT_ROWS:
        phase, gap, wall_count = direct_full_row_control(a, b)
        direct_rows.append(((a, b), pair(phase), pair(gap), wall_count))
    direct_rows = tuple(direct_rows)

    ledger = {
        "theorem": "THM-4136",
        "statement": "for all distinct positive odd a,b, 2U union {a,b} is 1/14-safe",
        "U": U,
        "normalization": "a=pt,b=qt with odd t=gcd(a,b) and coprime odd p<q",
        "body": body,
        "quotient": quotient,
        "scale_gates": {
            "t_at_least_3_surplus": pair(3 * J_LENGTH - UNIVERSAL_BETA),
            "primitive_q_at_least_27_surplus": pair(J_LENGTH - Q(2, 7 * 27)),
        },
        "finite_clock_bank": finite,
        "fixed_sheet_hostile": {
            "w": pair(hostile_w),
            "sheet_zero": tuple(pair(value) for value in hostile_sheet_zero),
            "sheet_one": tuple(pair(value) for value in hostile_sheet_one),
        },
        "direct_full_row_controls": direct_rows,
        "scope": "fixed U only; arbitrary eleven-speed bodies, physical entry, and LRC14 remain open",
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    low_text = tuple(
        (ratio, tuple((fmt(Q(*left)), fmt(Q(*right))) for left, right in components), fmt(Q(*beta)))
        for ratio, components, beta in quotient["low"]
    )
    direct_text = tuple(
        (tails, f"{phase[0]}/{phase[1]}", f"{gap[0]}/{gap[1]}", wall_count)
        for tails, phase, gap, wall_count in direct_rows
    )
    print("THM4136_FIXED_BODY_UNIVERSAL_ODD_TAIL_COMPLETION_20260825")
    print("statement=every distinct positive odd a,b makes 2U union {a,b} 1/14-safe")
    print(f"U={U};J=[{fmt(J[0])},{fmt(J[1])}];length={fmt(J_LENGTH)};owners={body['owners']}")
    print(f"low_two_sheet_quotients={low_text}")
    print(
        f"universal_quotient_beta={fmt(UNIVERSAL_BETA)};"
        f"scan={quotient['scan']};large_ratio_rule=2/(7q) for q>=9"
    )
    print(
        f"scale_gates=t>=3 surplus:{fmt(3 * J_LENGTH - UNIVERSAL_BETA)};"
        f"primitive q>=27 surplus:{fmt(J_LENGTH - Q(2, 7 * 27))}"
    )
    print(
        f"finite_bank=pairs:{finite['pair_count']};clocks:{tuple(fmt(x) for x in CLOCKS)};"
        f"categories:{finite['category_counts']};minima:{tuple(fmt(Q(*x)) for x in finite['category_minima'])};"
        f"exceptions:{finite['exception_pairs']};pair_sha256:{finite['pair_digest']}"
    )
    print(
        f"fixed_sheet_hostile=w:{fmt(hostile_w)};"
        f"sheet0:{tuple(fmt(x) for x in hostile_sheet_zero)};"
        f"sheet1:{tuple(fmt(x) for x in hostile_sheet_one)}"
    )
    print(f"direct_full_row_controls={direct_text}")
    print("scope=fixed U only; arbitrary bodies, physical entry, and LRC14 remain open")
    print(f"semantic_sha256={semantic}")
    print(f"checks={CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
