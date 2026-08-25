#!/usr/bin/env python3
"""Clean-room exact referee for THM-4132.

This audit imports no primary implementation.  It reconstructs the fixed
body's closed safe components as the complement of all strict danger cells,
builds the two-sheet obstruction from its literal inequalities, and checks
the compact-to-open length gate.  Selected full rows are then audited on
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
DELTA = Q(1, 14)
EXPECTED_J = (Q(33, 70), Q(27, 56))
EXPECTED_C = ((Q(2, 21), Q(8, 63)), (Q(55, 63), Q(19, 21)))
LITERAL_SCALES = (1, 3, 5, 21, 23, 93, 95, 2001)
EXPECTED_SEMANTIC = "a1a708792f72909ea1ee517586dcc30f8dc16f46781005b4e3fc22edcf875f88"

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def fmt(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


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
    """All equality walls bounding ||v*x|| < 1/14, cut at zero."""
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for integer in range(speed):
            for sign in (-1, 1):
                walls.add(((Q(integer) + sign * DELTA) / speed) % 1)
    return tuple(sorted(walls))


def positive_closed_safe_components(
    speeds: tuple[int, ...],
) -> tuple[tuple[Q, ...], tuple[tuple[Q, Q], ...]]:
    """Complement strict danger cells, preserving closed equality walls."""
    walls = strict_danger_walls(speeds)
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        if not closed_safe(speeds, midpoint):
            continue
        require(closed_safe(speeds, left), ("safe-cell left endpoint", left, right))
        require(closed_safe(speeds, right), ("safe-cell right endpoint", left, right))
        cells.append((left, right))

    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left and closed_safe(speeds, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return walls, tuple(components)


def reconstruct_body_interval() -> dict[str, object]:
    walls, components = positive_closed_safe_components(U)
    measure = sum((right - left for left, right in components), Q(0))
    maximum = max(right - left for left, right in components)
    maximizers = tuple(component for component in components if component[1] - component[0] == maximum)
    left_maximizers = tuple(component for component in maximizers if component[0] < Q(1, 2))
    require(len(walls) == 244, "body wall count")
    require(len(components) == 28, "body positive-component count")
    require(measure == Q(5497, 64680), "body positive-safe measure")
    require(maximum == Q(3, 280), "body maximum component length")
    require(maximizers == (EXPECTED_J, (Q(29, 56), Q(37, 70))), "body maximum components")
    require(left_maximizers == (EXPECTED_J,), "reconstructed J")

    left, right = left_maximizers[0]
    left_owners = tuple(speed for speed in U if circle_distance(speed * left) == DELTA)
    right_owners = tuple(speed for speed in U if circle_distance(speed * right) == DELTA)
    require(left_owners == (15,), "J left owner")
    require(right_owners == (4,), "J right owner")
    require(closed_safe(U, left) and closed_safe(U, right), "J closed endpoints")
    require(closed_safe(U, (left + right) / 2), "J interior")
    return {
        "wall_count": len(walls),
        "component_count": len(components),
        "measure": pair(measure),
        "component_digest": digest(tuple((pair(a), pair(b)) for a, b in components)),
        "J": (pair(left), pair(right), pair(right - left)),
        "owners": (left_owners, right_owners),
        "maximizers": tuple((pair(a), pair(b)) for a, b in maximizers),
    }


def normalized_pair_bad(phase: Q) -> bool:
    return circle_distance(phase) < DELTA or circle_distance(9 * phase) < DELTA


def both_sheets_bad(w: Q) -> bool:
    z = w / 2
    return normalized_pair_bad(z) and normalized_pair_bad(z + Q(1, 2))


def two_sheet_walls() -> tuple[Q, ...]:
    """Solve ||v*(w/2+s/2)||=1/14 for v=1,9 and sheets s=0,1."""
    walls = {Q(0), Q(1)}
    for sheet_shift in (Q(0), Q(1, 2)):
        for speed in (1, 9):
            for integer in range(speed):
                for sign in (-1, 1):
                    phase_wall = (Q(integer) + sign * DELTA) / speed
                    walls.add((2 * (phase_wall - sheet_shift)) % 1)
    return tuple(sorted(walls))


def reconstruct_two_sheet_obstruction() -> dict[str, object]:
    walls = two_sheet_walls()
    active_cells = tuple(
        (left, right)
        for left, right in zip(walls, walls[1:])
        if both_sheets_bad((left + right) / 2)
    )
    require(not both_sheets_bad(Q(0)), "quotient nonwrap at zero")

    components: list[tuple[Q, Q]] = []
    for left, right in active_cells:
        if components and components[-1][1] == left and both_sheets_bad(left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    obstruction = tuple(components)
    require(len(walls) == 22, "two-sheet wall count")
    require(obstruction == EXPECTED_C, "two-sheet obstruction")
    for left, right in obstruction:
        require(not both_sheets_bad(left), ("open obstruction left", left))
        require(not both_sheets_bad(right), ("open obstruction right", right))
        require(both_sheets_bad((left + right) / 2), ("obstruction interior", left, right))

    widths = tuple(right - left for left, right in obstruction)
    beta = max(widths)
    measure = sum(widths, Q(0))
    require(widths == (Q(2, 63), Q(2, 63)), "obstruction widths")
    require(beta == Q(2, 63) and measure == Q(4, 63), "obstruction beta and measure")
    return {
        "wall_count": len(walls),
        "components": tuple((pair(a), pair(b)) for a, b in obstruction),
        "beta": pair(beta),
        "measure": pair(measure),
    }


def ceiling(value: Q) -> int:
    return -((-value.numerator) // value.denominator)


def compact_to_open_audit(j_length: Q, beta: Q) -> dict[str, object]:
    """A compact source arc fits in an open target arc only at strict length."""
    ratio = beta / j_length
    first_integer_gate = ceiling(ratio)
    first_surplus = first_integer_gate * j_length - beta
    wrap_threshold = ceiling(1 / j_length)
    require(ratio == Q(80, 27), "critical real scale")
    require(first_integer_gate == 3, "first integer compact-to-open gate")
    require(first_surplus == Q(1, 2520) > 0, "first integer surplus")
    require(2 * j_length < beta < 3 * j_length, "sharp neighboring integer controls")
    require(not (beta < beta), "equal compact/open lengths do not fit")
    for left, right in EXPECTED_C:
        require(not both_sheets_bad(left) and not both_sheets_bad(right),
                ("equal-length endpoint witness", left, right))
    require(wrap_threshold == 94, "first wrapping scale")
    return {
        "critical_real_scale": pair(ratio),
        "first_integer_gate": first_integer_gate,
        "first_surplus": pair(first_surplus),
        "equality_closes": True,
        "wrap_threshold": wrap_threshold,
    }


def direct_full_row_control(t: int) -> tuple[Q, Q, int]:
    """Find a strict witness from the full row's wall complement, not J or C."""
    require(t > 0 and t % 2 == 1, ("control scale", t))
    body = tuple(2 * speed for speed in U)
    speeds = body + (t, 9 * t)
    require(len(set(speeds)) == 13, ("distinct full row", t))
    require(reduce(gcd, speeds) == 1, ("primitive full row", t))
    walls = strict_danger_walls(speeds)
    for left, right in zip(walls, walls[1:]):
        phase = (left + right) / 2
        gap = clearance(speeds, phase)
        if gap > DELTA:
            return phase, gap, len(walls)
    raise RuntimeError(f"no direct full-row witness at t={t}")


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    require(len(U) == len(set(U)) == 11 and min(U) == 1 and max(U) == 22,
            "fixed body universe")
    even_body = tuple(2 * speed for speed in U)
    require(all(speed % 2 == 0 for speed in even_body), "even body")
    require(reduce(gcd, even_body) == 2, "body gcd two")

    body = reconstruct_body_interval()
    quotient = reconstruct_two_sheet_obstruction()
    j_length = Q(*body["J"][2])
    beta = Q(*quotient["beta"])
    gate = compact_to_open_audit(j_length, beta)

    t1_phase = Q(1, 13)
    t1_speeds = even_body + (1, 9)
    t1_gap = clearance(t1_speeds, t1_phase)
    require(t1_gap == Q(1, 13) > DELTA, "t=1 named clock")

    hostile_w = Q(1, 9)
    hostile_z = hostile_w / 2
    sheet_zero = (circle_distance(hostile_z), circle_distance(9 * hostile_z))
    sheet_one = (
        circle_distance(hostile_z + Q(1, 2)),
        circle_distance(9 * (hostile_z + Q(1, 2))),
    )
    require(sheet_zero == (Q(1, 18), Q(1, 2)), "hostile first sheet")
    require(sheet_one == (Q(4, 9), Q(0)), "hostile second sheet")
    require(sheet_zero[0] < DELTA and sheet_zero[1] > DELTA, "first sheet owner one")
    require(sheet_one[0] > DELTA and sheet_one[1] < DELTA, "second sheet owner nine")
    require(both_sheets_bad(hostile_w), "hostile belongs to obstruction")

    literal_rows = []
    for t in LITERAL_SCALES:
        phase, gap, wall_count = direct_full_row_control(t)
        literal_rows.append((t, pair(phase), pair(gap), wall_count))
    literal_rows = tuple(literal_rows)
    require(literal_rows[0][1] == (89, 1176), "independent t=1 first strict cell")
    require(literal_rows[1][1] == (233, 4032), "independent t=3 first strict cell")

    ledger = {
        "theorem": "THM-4132",
        "statement": "for every positive odd t, 2U union {t,9t} is 1/14-safe",
        "U": U,
        "normalization": {
            "even_body_gcd": 2,
            "odd_tails_disjoint_by_parity": True,
            "full_row_primitive_for_odd_t": True,
            "tail_ratio": (1, 9),
        },
        "body_strict_danger_complement": body,
        "two_sheet_obstruction": quotient,
        "compact_to_open": gate,
        "t1_clock": (pair(t1_phase), pair(t1_gap)),
        "fixed_sheet_hostile": {
            "w": pair(hostile_w),
            "z": pair(hostile_z),
            "sheet_zero_gaps": tuple(pair(value) for value in sheet_zero),
            "sheet_one_gaps": tuple(pair(value) for value in sheet_one),
        },
        "direct_full_row_controls": literal_rows,
        "scope": (
            "fixed U-body only; no arbitrary-body scale-two closure, physical entry, "
            "even-t normalization, or LRC14"
        ),
    }
    semantic = digest(ledger)
    require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    literal_text = tuple(
        (t, f"{phase[0]}/{phase[1]}", f"{gap[0]}/{gap[1]}", wall_count)
        for t, phase, gap, wall_count in literal_rows
    )
    print("THM4132_FIXED_BODY_SCALE_TWO_COMPLETION_INDEPENDENT_AUDIT_20260825")
    print("route=strict-danger complement + literal two-sheet inequalities + full-row wall controls")
    print(
        "body_complement="
        f"walls:{body['wall_count']};positive_components:{body['component_count']};"
        f"measure:{fmt(Q(*body['measure']))};component_sha256:{body['component_digest']}"
    )
    print(
        f"J=[{fmt(Q(*body['J'][0]))},{fmt(Q(*body['J'][1]))}];"
        f"length={fmt(j_length)};owners={body['owners']}"
    )
    quotient_text = tuple((fmt(Q(*left)), fmt(Q(*right))) for left, right in quotient["components"])
    print(
        f"two_sheet_quotient=walls:{quotient['wall_count']};components:{quotient_text};"
        f"beta:{fmt(beta)};measure:{fmt(Q(*quotient['measure']))};endpoints:open"
    )
    print(
        f"compact_to_open=critical_scale:{fmt(Q(*gate['critical_real_scale']))};"
        f"first_integer:{gate['first_integer_gate']};surplus:{fmt(Q(*gate['first_surplus']))};"
        f"equality_closes:{int(gate['equality_closes'])};wrap_threshold:{gate['wrap_threshold']}"
    )
    print(f"t1_clock={fmt(t1_phase)};clearance={fmt(t1_gap)}")
    print(
        f"fixed_sheet_hostile=w:{fmt(hostile_w)};z:{fmt(hostile_z)};"
        f"sheet0:{tuple(fmt(value) for value in sheet_zero)};"
        f"sheet1:{tuple(fmt(value) for value in sheet_one)}"
    )
    print(f"direct_full_row_controls={literal_text}")
    print("normalization=odd t gives 13 distinct speeds and primitive full row")
    print("scope=fixed U-body only; arbitrary bodies, physical entry, even-t normalization, and LRC14 remain open")
    print(f"semantic_sha256={semantic}")
    print(f"checks={CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
