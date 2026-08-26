#!/usr/bin/env python3
"""Clean-room exact audit for THM-4156.

This referee does not import or inspect the primary implementation.  It
constructs the pool safe set by successively intersecting closed safe-tooth
unions, one speed at a time.  The global endpoint set is used only for the
separate wall-count checksum, never to classify midpoint cells.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = F(1, 14)
HAAR_THRESHOLD = F(4, 63)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
ANCHORS = (120, 126, 143)
EXPECTED_MEASURE = F(298133356159, 4560289854120)
EXPECTED_SURPLUS = F(8591143199, 4560289854120)
EXPECTED_MAX_WIDTH = F(37, 25520)
EXPECTED_MAXIMIZERS = (
    (F(393, 1232), F(1301, 4060)),
    (F(2759, 4060), F(839, 1232)),
)
EXPECTED_SEMANTIC = "97954c5f353ca0977c87d1649abb646a43af4cdb636bcc211861a12129b1c415"

CHECKS = 0


def verify(predicate: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(f"independent audit failed: {label}")


def qpair(value: F) -> tuple[int, int]:
    return value.numerator, value.denominator


def circle_gap(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: F) -> F:
    verify(bool(speeds), ("nonempty clearance packet", phase))
    return min(circle_gap(speed * phase) for speed in speeds)


def merge_closed(
    intervals: list[tuple[F, F]] | tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    merged: list[tuple[F, F]] = []
    for left, right in sorted(intervals):
        if left > right:
            continue
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def intersect_closed(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    pieces: list[tuple[F, F]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left <= right:
            pieces.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge_closed(pieces)


def safe_teeth(speed: int) -> tuple[tuple[F, F], ...]:
    verify(speed > 0, ("positive speed", speed))
    return tuple(
        ((F(integer) + DELTA) / speed,
         (F(integer + 1) - DELTA) / speed)
        for integer in range(speed)
    )


def interval_measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def point_in_union(point: F, intervals: tuple[tuple[F, F], ...]) -> bool:
    return any(left <= point <= right for left, right in intervals)


def coverage_audit() -> dict[str, object]:
    verify(len(POOL) == len(set(POOL)) == 30, "pool cardinality")
    verify(tuple(sorted(POOL)) == POOL, "pool sorted")
    verify(len(ANCHORS) == len(set(ANCHORS)) == 3, "anchor cardinality")
    verify(set(ANCHORS) < set(POOL), "anchors are a proper pool subset")

    coverage = tuple(
        (divisor, tuple(value for value in ANCHORS if value % divisor == 0))
        for divisor in range(2, 15)
    )
    expected = (
        (2, (120, 126)),
        (3, (120, 126)),
        (4, (120,)),
        (5, (120,)),
        (6, (120, 126)),
        (7, (126,)),
        (8, (120,)),
        (9, (126,)),
        (10, (120,)),
        (11, (143,)),
        (12, (120,)),
        (13, (143,)),
        (14, (126,)),
    )
    verify(coverage == expected, "anchor divisor coverage 2 through 14")
    verify(all(owners for _, owners in coverage), "no missing divisor")

    free = tuple(value for value in POOL if value not in ANCHORS)
    verify(len(free) == 27, "free-label cardinality")
    verify(comb(len(free), 8) == 2_220_075, "hereditary family universe")
    return {
        "pool": POOL,
        "anchors": ANCHORS,
        "coverage": coverage,
        "free": free,
    }


def geometry_audit() -> tuple[dict[str, object], tuple[tuple[F, F], ...]]:
    current = ((F(0), F(1)),)
    previous_measure = F(1)
    processed: list[int] = []
    trace: list[tuple[object, ...]] = []

    for speed in POOL:
        teeth = safe_teeth(speed)
        verify(len(teeth) == speed, ("tooth count", speed))
        verify(interval_measure(teeth) == F(6, 7),
               ("one-speed safe measure", speed))
        current = intersect_closed(current, teeth)
        processed.append(speed)
        measure = interval_measure(current)
        verify(measure <= previous_measure, ("monotone measure", speed))
        verify(all(left <= right for left, right in current),
               ("ordered endpoints", speed))
        verify(all(current[index][1] < current[index + 1][0]
                   for index in range(len(current) - 1)),
               ("merged disjoint components", speed))
        for left, right in current:
            verify(clearance(tuple(processed), left) >= DELTA,
                   ("closed left endpoint", speed, left))
            verify(clearance(tuple(processed), right) >= DELTA,
                   ("closed right endpoint", speed, right))
        trace.append((speed, len(current), qpair(measure)))
        previous_measure = measure

    verify(all(left < right for left, right in current),
           "no zero-length final components")
    measure = interval_measure(current)
    maximum = max(right - left for left, right in current)
    maximizers = tuple(
        component for component in current
        if component[1] - component[0] == maximum
    )
    reflection = tuple(sorted((1 - right, 1 - left) for left, right in current))
    verify(reflection == current, "reflection symmetry")
    verify(len(current) == 150, "final component count")
    verify(measure == EXPECTED_MEASURE, "final Haar measure")
    verify(measure - HAAR_THRESHOLD == EXPECTED_SURPLUS,
           "Haar-threshold surplus")
    verify(maximum == EXPECTED_MAX_WIDTH, "maximum component width")
    verify(maximizers == EXPECTED_MAXIMIZERS,
           "maximum reflected component pair")

    # This endpoint census is deliberately separate from the iterative
    # geometry construction above.  No midpoint classification uses it.
    walls = {F(0), F(1)}
    for speed in POOL:
        for left, right in safe_teeth(speed):
            walls.add(left)
            walls.add(right)
    verify(len(walls) == 7_134, "distinct safe-wall count")

    return ({
        "trace": tuple(trace),
        "walls": len(walls),
        "components": len(current),
        "measure": qpair(measure),
        "threshold": qpair(HAAR_THRESHOLD),
        "surplus": qpair(measure - HAAR_THRESHOLD),
        "max_width": qpair(maximum),
        "maximizers": tuple(
            (qpair(left), qpair(right)) for left, right in maximizers
        ),
    }, current)


def control_audit(
    safe: tuple[tuple[F, F], ...],
) -> dict[str, object]:
    body_phase = F(1, 112)
    physical_phase = F(113, 224)
    tails = (1, 3)
    verify((2 * physical_phase) % 1 == body_phase,
           "positive-control physical lift")
    verify(point_in_union(body_phase, safe), "positive body phase in safe set")

    body_speeds = tuple(2 * value for value in POOL)
    body_gap = clearance(body_speeds, physical_phase)
    body_owners = tuple(
        value for value in POOL
        if circle_gap(2 * value * physical_phase) == body_gap
    )
    tail_gaps = tuple(circle_gap(tail * physical_phase) for tail in tails)
    full_gap = min(body_gap, *tail_gaps)
    verify(body_gap == DELTA, "positive-control body gap")
    verify(body_owners == (8, 120), "positive-control body owners")
    verify(tail_gaps == (F(111, 224), F(109, 224)),
           "positive-control tail gaps")
    verify(full_gap == DELTA, "positive-control full-row gap")

    hostile_phase = F(1, 12)
    hostile_zeros = tuple(
        value for value in POOL
        if circle_gap(2 * value * hostile_phase) == 0
    )
    hostile_anchor_zeros = tuple(
        value for value in ANCHORS
        if circle_gap(2 * value * hostile_phase) == 0
    )
    verify(hostile_zeros ==
           (30, 42, 60, 84, 120, 126, 132, 168, 240, 252, 264),
           "x=1/12 zero-owner list")
    verify(hostile_anchor_zeros == (120, 126),
           "x=1/12 mandatory zero anchors")
    verify(clearance(body_speeds, hostile_phase) == 0,
           "x=1/12 hostile body gap")

    return {
        "positive": (
            qpair(body_phase), qpair(physical_phase), tails,
            qpair(body_gap), body_owners,
            tuple(qpair(value) for value in tail_gaps), qpair(full_gap),
        ),
        "hostile_mod12": (
            qpair(hostile_phase), hostile_zeros, hostile_anchor_zeros,
            (0, 1),
        ),
    }


def choose(total: int, selected: int) -> int:
    if selected < 0 or selected > total:
        return 0
    return comb(total, selected)


def width_gate(minimum: int, maximum: int) -> bool:
    return 27 * (13 * minimum - maximum) >= 4 * minimum * maximum


def thm4151_affine_gate(minimum: int, maximum: int) -> bool:
    # Exact predicate of the proved THM-4151 affine gate.
    return (
        minimum >= 3
        and (12 * minimum + 1) * (13 * minimum - maximum)
        >= 4 * minimum * maximum
    )


def family_count_audit(free: tuple[int, ...]) -> dict[str, object]:
    direct_histogram: dict[tuple[int, int], int] = {}
    direct = {
        "total": 0,
        "thm4148": 0,
        "thm4151_affine_gate": 0,
        "both": 0,
        "thm4148_only": 0,
        "thm4151_affine_gate_only": 0,
        "neither": 0,
    }
    width_best: int | None = None
    width_witnesses: set[tuple[int, int]] = set()
    affine_gate_pass_min: int | None = None
    affine_gate_pass_witnesses: set[tuple[int, int]] = set()
    affine_gate_fail_max: int | None = None
    affine_gate_fail_witnesses: set[tuple[int, int]] = set()

    for selected in combinations(free, 8):
        minimum = min(ANCHORS[0], selected[0])
        maximum = max(ANCHORS[-1], selected[-1])
        direct_histogram[(minimum, maximum)] = (
            direct_histogram.get((minimum, maximum), 0) + 1
        )
        first = width_gate(minimum, maximum)
        second = thm4151_affine_gate(minimum, maximum)
        direct["total"] += 1
        direct["thm4148"] += int(first)
        direct["thm4151_affine_gate"] += int(second)
        direct["both"] += int(first and second)
        direct["thm4148_only"] += int(first and not second)
        direct["thm4151_affine_gate_only"] += int(second and not first)
        direct["neither"] += int(not first and not second)

        width_margin = 27 * (13 * minimum - maximum) - 4 * minimum * maximum
        if width_best is None or width_margin > width_best:
            width_best = width_margin
            width_witnesses = {(minimum, maximum)}
        elif width_margin == width_best:
            width_witnesses.add((minimum, maximum))

        affine_gate_margin = (
            (12 * minimum + 1) * (13 * minimum - maximum)
            - 4 * minimum * maximum
        )
        if second:
            if (affine_gate_pass_min is None
                    or affine_gate_margin < affine_gate_pass_min):
                affine_gate_pass_min = affine_gate_margin
                affine_gate_pass_witnesses = {(minimum, maximum)}
            elif affine_gate_margin == affine_gate_pass_min:
                affine_gate_pass_witnesses.add((minimum, maximum))
        else:
            if (affine_gate_fail_max is None
                    or affine_gate_margin > affine_gate_fail_max):
                affine_gate_fail_max = affine_gate_margin
                affine_gate_fail_witnesses = {(minimum, maximum)}
            elif affine_gate_margin == affine_gate_fail_max:
                affine_gate_fail_witnesses.add((minimum, maximum))

    # A second count uses only the two endpoint choices and binomial interior
    # multiplicities.  The anchors force min<=120 and max>=143.
    low_endpoints = tuple(value for value in free if value < ANCHORS[0]) + (
        ANCHORS[0],
    )
    high_endpoints = (ANCHORS[-1],) + tuple(
        value for value in free if value > ANCHORS[-1]
    )
    formula_histogram: dict[tuple[int, int], int] = {}
    formula = {key: 0 for key in direct}
    for minimum in low_endpoints:
        for maximum in high_endpoints:
            required = int(minimum < ANCHORS[0]) + int(maximum > ANCHORS[-1])
            interior = sum(1 for value in free if minimum < value < maximum)
            ways = choose(interior, 8 - required)
            if ways == 0:
                continue
            formula_histogram[(minimum, maximum)] = ways
            first = width_gate(minimum, maximum)
            second = thm4151_affine_gate(minimum, maximum)
            formula["total"] += ways
            formula["thm4148"] += ways * int(first)
            formula["thm4151_affine_gate"] += ways * int(second)
            formula["both"] += ways * int(first and second)
            formula["thm4148_only"] += ways * int(first and not second)
            formula["thm4151_affine_gate_only"] += ways * int(second and not first)
            formula["neither"] += ways * int(not first and not second)

    verify(direct_histogram == formula_histogram,
           "direct/formula min-max histogram")
    verify(direct == formula, "direct/formula gate counts")
    expected = {
        "total": 2_220_075,
        "thm4148": 0,
        "thm4151_affine_gate": 344_366,
        "both": 0,
        "thm4148_only": 0,
        "thm4151_affine_gate_only": 344_366,
        "neither": 1_875_709,
    }
    verify(direct == expected, "complete family gate census")
    verify(width_best == -5_629, "sharp THM-4148 hostile margin")
    verify(tuple(sorted(width_witnesses)) == ((8, 143),),
           "sharp THM-4148 hostile endpoint pair")
    verify(affine_gate_pass_min == 350,
           "THM-4151 affine gate smallest passing margin")
    verify(tuple(sorted(affine_gate_pass_witnesses)) == ((15, 145),),
           "THM-4151 affine gate smallest passing endpoint pair")
    verify(affine_gate_fail_max == -3_032,
           "THM-4151 affine gate closest failing margin")
    verify(tuple(sorted(affine_gate_fail_witnesses)) == ((16, 168),),
           "THM-4151 affine gate closest failing endpoint pair")

    return {
        "counts": tuple(sorted(direct.items())),
        "outside": (
            ("thm4148", direct["total"] - direct["thm4148"]),
            ("thm4151_affine_gate",
             direct["total"] - direct["thm4151_affine_gate"]),
        ),
        "histogram": tuple(
            (minimum, maximum, ways)
            for (minimum, maximum), ways in sorted(direct_histogram.items())
        ),
        "margins": (
            ("thm4148_best", width_best, tuple(sorted(width_witnesses))),
            ("thm4151_affine_gate_pass_min", affine_gate_pass_min,
             tuple(sorted(affine_gate_pass_witnesses))),
            ("thm4151_affine_gate_fail_max", affine_gate_fail_max,
             tuple(sorted(affine_gate_fail_witnesses))),
        ),
    }


def main() -> None:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    verify(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
           "no assert statements")
    verify(not any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "float"
        for node in ast.walk(tree)
    ), "no float calls")

    coverage = coverage_audit()
    geometry, safe = geometry_audit()
    controls = control_audit(safe)
    families = family_count_audit(coverage["free"])
    semantic = {
        "delta": qpair(DELTA),
        "coverage": coverage,
        "geometry": geometry,
        "controls": controls,
        "families": families,
    }
    semantic_hash = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC != "PENDING":
        verify(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    coverage_rows = coverage["coverage"]
    counts = dict(families["counts"])
    print("LRC14_DIVISOR_COMPLETE_ANCHOR_POOL_HAAR_THM4156_INDEPENDENT_20260825")
    print("status=ACCEPT;scope=clean-room iterative closed-interval audit;LRC14=OPEN")
    print(f"pool_size={len(POOL)};anchors={ANCHORS};free={len(coverage['free'])}")
    print(f"divisor_coverage={coverage_rows}")
    print(
        "geometry="
        f"(walls={geometry['walls']},components={geometry['components']},"
        f"measure={EXPECTED_MEASURE},surplus={EXPECTED_SURPLUS},"
        f"max_width={EXPECTED_MAX_WIDTH})"
    )
    print(f"maximizers={geometry['maximizers']}")
    print(f"iterative_trace={geometry['trace']}")
    print(f"positive_control={controls['positive']}")
    print(f"x_1_over_12_hostile={controls['hostile_mod12']}")
    print(
        "family_counts="
        f"(total={counts['total']},thm4148_proved={counts['thm4148']},"
        f"outside_thm4148={counts['total']-counts['thm4148']},"
        "thm4151_affine_gate="
        f"{counts['thm4151_affine_gate']},"
        "outside_thm4151_affine_gate="
        f"{counts['total']-counts['thm4151_affine_gate']},"
        "thm4151_affine_gate_only="
        f"{counts['thm4151_affine_gate_only']},both={counts['both']},"
        f"neither={counts['neither']})"
    )
    print(f"gate_margins={families['margins']}")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
