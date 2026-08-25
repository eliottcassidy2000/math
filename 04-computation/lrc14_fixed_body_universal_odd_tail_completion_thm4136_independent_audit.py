#!/usr/bin/env python3
"""Clean-room exact audit for THM-4136.

This companion imports no project code.  It independently checks the proposed
fixed-body theorem for

    U = (1,4,6,8,10,12,14,15,16,18,22).

For odd tails a=pt and b=qt, with t=gcd(a,b), the two physical lifts of a
body-safe quotient phase fail together on an open set C_(p,q).  The proof
uses three disjoint pieces of evidence:

* an exact wall-cell reconstruction of C_(p,q), with strict wall activity;
* the symbolic q-tooth containment bound beta(C_(p,q)) <= 2/(7q) for q>=9,
  stress-tested on a broad primitive odd bank; and
* an exact 68-pair low-ratio census with three literal physical clocks.

The computation verifies the consequence (all thirteen clearances at each
clock), not only the quotient model.  It uses only Python's standard library.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
BODY_ARC = (Q(33, 70), Q(27, 56))
BODY_ARC_LENGTH = Q(3, 280)

CLOCK_NON13 = Q(89, 1176)
CLOCK_13_MIDDLE = Q(181, 4704)
CLOCK_13_EDGE = Q(431, 4480)

EXPECTED_SEMANTIC_SHA256 = "92b303d4737b146a0ff31de4c13fed920e70531d486495e31793128e5ca4482e"

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def fmt(value: Q) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(circle_distance(speed * phase) for speed in speeds)


def pair_bad(z: Q, p: int, q: int) -> bool:
    return circle_distance(p * z) < DELTA or circle_distance(q * z) < DELTA


def both_lifts_bad(w: Q, p: int, q: int) -> bool:
    """Predicate on w=2z for failure of both physical two-sheet lifts."""
    z = (w % 1) / 2
    return pair_bad(z, p, q) and pair_bad(z + Q(1, 2), p, q)


def quotient_walls(p: int, q: int) -> tuple[Q, ...]:
    """All possible strict-danger boundaries after the quotient w=2z."""
    return tuple(sorted(
        {Q(0)}
        | {
            (2 * (Q(k) + sign * DELTA) / speed) % 1
            for speed in (p, q)
            for k in range(speed)
            for sign in (-1, 1)
        }
    ))


def quotient_components(p: int, q: int) -> tuple[tuple[Q, Q, Q], ...]:
    """Open circular components of C_(p,q), retaining inactive punctures.

    A component is returned as (one wall, the next wall after traversing the
    component, metric length).  The endpoint labels are diagnostic only; the
    metric length is the theorem input.  Adjacent active cells are joined only
    when their common wall is itself strictly active.
    """
    walls = quotient_walls(p, q)
    count = len(walls)
    cell_active: list[bool] = []
    for index, left in enumerate(walls):
        right = walls[(index + 1) % count]
        right_lift = right if right > left else right + 1
        midpoint = ((left + right_lift) / 2) % 1
        cell_active.append(both_lifts_bad(midpoint, p, q))
    wall_active = [both_lifts_bad(wall, p, q) for wall in walls]

    unseen = {index for index, active in enumerate(cell_active) if active}
    components: list[tuple[Q, Q, Q]] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        cells: set[int] = set()
        while stack:
            index = stack.pop()
            cells.add(index)
            neighbors = (
                ((index - 1) % count, index),
                ((index + 1) % count, (index + 1) % count),
            )
            for neighbor, shared_wall in neighbors:
                if (neighbor in unseen and cell_active[neighbor]
                        and wall_active[shared_wall]):
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        length = sum(
            ((walls[(index + 1) % count] - walls[index]) % 1
             for index in cells),
            Q(0),
        )
        first = min(cells)
        # Endpoint labels are made deterministic even for a wrapping component.
        components.append((walls[first], (walls[first] + length) % 1, length))
    return tuple(sorted(components, key=lambda row: (row[2], row[0], row[1])))


def quotient_width(p: int, q: int) -> Q:
    components = quotient_components(p, q)
    return max((length for _, _, length in components), default=Q(0))


def primitive_odd_pairs(max_q: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (p, q)
        for q in range(3, max_q + 1, 2)
        for p in range(1, q, 2)
        if math.gcd(p, q) == 1
    )


def body_arc_controls() -> dict[str, object]:
    left, right = BODY_ARC
    require(right - left == BODY_ARC_LENGTH, "body arc length")
    require(clearance(BODY, left) == DELTA, "body left endpoint safe")
    require(clearance(BODY, right) == DELTA, "body right endpoint safe")
    require(
        tuple(speed for speed in BODY if circle_distance(speed * left) == DELTA)
        == (15,),
        "body left owner",
    )
    require(
        tuple(speed for speed in BODY if circle_distance(speed * right) == DELTA)
        == (4,),
        "body right owner",
    )
    for numerator in range(101):
        point = left + (right - left) * Q(numerator, 100)
        require(clearance(BODY, point) >= DELTA, f"body arc sample {numerator}")
    return {
        "body": list(BODY),
        "arc": [fmt(left), fmt(right)],
        "length": fmt(BODY_ARC_LENGTH),
        "owners": [15, 4],
    }


def low_width_controls() -> dict[str, object]:
    expected = {
        (1, 3): (),
        (1, 5): (),
        (3, 5): (Q(1, 105), Q(1, 105)),
        (1, 7): (Q(1, 49), Q(1, 49)),
        (3, 7): (Q(1, 49), Q(1, 49)),
        (5, 7): (Q(1, 49), Q(1, 49)),
    }
    observed: dict[tuple[int, int], tuple[Q, ...]] = {}
    for pair, target in expected.items():
        lengths = tuple(length for _, _, length in quotient_components(*pair))
        require(lengths == target, f"low quotient topology {pair}")
        observed[pair] = lengths
    low_max = max((max(row, default=Q(0)) for row in observed.values()))
    require(low_max == Q(1, 49) < Q(2, 63), "low-pair width ceiling")

    exceptional = quotient_components(1, 9)
    exceptional_lengths = tuple(length for _, _, length in exceptional)
    require(exceptional_lengths == (Q(2, 63), Q(2, 63)), "(1,9) widths")
    endpoints = sorted({endpoint for left, right, _ in exceptional for endpoint in (left, right)})
    require(len(endpoints) == 4, "(1,9) four open endpoints")
    require(all(not both_lifts_bad(endpoint, 1, 9) for endpoint in endpoints),
            "strict endpoint equality controls")
    require(all(
        both_lifts_bad(((left + length / 2) % 1), 1, 9)
        for left, _, length in exceptional
    ), "(1,9) interior hostile controls")

    return {
        "low_components": {
            f"{p},{q}": [fmt(length) for length in lengths]
            for (p, q), lengths in sorted(observed.items())
        },
        "low_max": fmt(low_max),
        "exceptional_1_9": [fmt(length) for length in exceptional_lengths],
        "strict_endpoints": [fmt(value) for value in endpoints],
    }


def q_tooth_bound_controls() -> dict[str, object]:
    # For odd q, unshifted and half-shifted q-tooth centres have circular
    # separation 1/(2q), while their two open radii total only 1/(7q).
    # The cross terms X=D_p intersect (D_q-1/2) and X+1/2 are antipodal.
    # Doubling identifies those copies; it does not concatenate them.  Thus a
    # quotient component lies in the doubled image of one q-tooth.
    require(Q(1, 2) - Q(1, 7) == Q(5, 14), "q-tooth separation surplus")

    exhaustive = primitive_odd_pairs(101)
    scanned: list[tuple[int, int, Q, Q]] = []
    transcript = sha256()
    for p, q in exhaustive:
        if q < 9:
            continue
        width = quotient_width(p, q)
        bound = Q(2, 7 * q)
        require(width <= bound, f"q-tooth bound {(p, q)}")
        scanned.append((p, q, width, bound))
        transcript.update(f"{p},{q},{fmt(width)},{fmt(bound)}\n".encode("ascii"))

    # Sparse large-q hostiles exercise highly unbalanced, near-diagonal, and
    # intermediate primitive ratios without making replay needlessly slow.
    large_qs = (127, 151, 181, 211, 251, 301, 401, 503, 701, 1001)
    sparse: list[tuple[int, int, Q, Q]] = []
    for q in large_qs:
        candidates = {
            1,
            3,
            5,
            q - 2,
            q - 4,
            max(1, 2 * (q // 6) + 1),
            max(1, 2 * (q // 4) + 1),
        }
        for p in sorted(candidates):
            if not (0 < p < q and p % 2 == 1 and math.gcd(p, q) == 1):
                continue
            width = quotient_width(p, q)
            bound = Q(2, 7 * q)
            require(width <= bound, f"large q-tooth hostile {(p, q)}")
            sparse.append((p, q, width, bound))
            transcript.update(f"{p},{q},{fmt(width)},{fmt(bound)}\n".encode("ascii"))

    equality_rows = tuple(
        (p, q) for p, q, width, bound in scanned + sparse if width == bound
    )
    require((1, 9) in equality_rows, "sharp q-tooth row")
    global_max = max(width for _, _, width, _ in scanned)
    global_max_rows = tuple(
        (p, q) for p, q, width, _ in scanned + sparse if width == global_max
    )
    require(global_max == Q(2, 63), "scanned global width")
    require(global_max_rows == ((1, 9),), "unique scanned global width row")
    require(all(width <= Q(2, 63) for _, _, width, _ in scanned + sparse),
            "universal 2/63 stress ceiling")

    width_histogram = Counter(fmt(width) for _, _, width, _ in scanned)
    return {
        "symbolic": {
            "q_tooth_width_z": "1/(7q)",
            "opposite_grid_center_gap_z": "1/(2q)",
            "strict_grid_gap_z": "5/(14q)",
            "doubled_component_bound_w": "2/(7q)",
        },
        "exhaustive_q_max": 101,
        "exhaustive_rows": len(scanned),
        "large_q_max": max(large_qs),
        "large_hostile_rows": len(sparse),
        "bound_equality_count": len(equality_rows),
        "bound_equality_first": [list(pair) for pair in equality_rows[:8]],
        "bound_equality_last": [list(pair) for pair in equality_rows[-8:]],
        "global_max": fmt(global_max),
        "global_max_rows": [list(pair) for pair in global_max_rows],
        "width_histogram": sorted(width_histogram.items()),
        "scan_transcript_sha256": transcript.hexdigest(),
    }


def residual_clock_controls() -> dict[str, object]:
    residual = primitive_odd_pairs(25)
    require(len(residual) == 68 and len(set(residual)) == 68,
            "68 residual primitive odd pairs")

    rows: list[tuple[str, int, int, Q, Q]] = []
    for p, q in residual:
        if 13 not in (p, q):
            category = "non13"
            phase = CLOCK_NON13
        else:
            other = p if q == 13 else q
            if other in (1, 25):
                category = "13_edge"
                phase = CLOCK_13_EDGE
            else:
                category = "13_middle"
                phase = CLOCK_13_MIDDLE
        speeds = tuple(2 * speed for speed in BODY) + (p, q)
        require(len(speeds) == 13 and len(set(speeds)) == 13,
                f"distinct residual row {(p, q)}")
        gap = clearance(speeds, phase)
        require(gap >= DELTA, f"literal clock clearance {(p, q)}")
        rows.append((category, p, q, phase, gap))

    expected_counts = {"non13": 56, "13_middle": 10, "13_edge": 2}
    counts = Counter(category for category, _, _, _, _ in rows)
    require(dict(counts) == expected_counts, "clock category split 56/10/2")

    expected_minima = {
        "non13": Q(89, 1176),
        "13_middle": Q(15, 196),
        "13_edge": Q(17, 224),
    }
    minima: dict[str, Q] = {}
    minimizers: dict[str, tuple[tuple[int, int], ...]] = {}
    for category in expected_counts:
        minimum = min(gap for key, _, _, _, gap in rows if key == category)
        require(minimum == expected_minima[category], f"clock minimum {category}")
        require(minimum > DELTA, f"strict clock surplus {category}")
        minima[category] = minimum
        minimizers[category] = tuple(
            (p, q) for key, p, q, _, gap in rows
            if key == category and gap == minimum
        )

    return {
        "residual_q_max": 25,
        "residual_pairs": [list(pair) for pair in residual],
        "counts": expected_counts,
        "clocks": {
            "non13": fmt(CLOCK_NON13),
            "13_middle": fmt(CLOCK_13_MIDDLE),
            "13_edge": fmt(CLOCK_13_EDGE),
        },
        "minima": {key: fmt(value) for key, value in minima.items()},
        "minimizers": {
            key: [list(pair) for pair in pairs]
            for key, pairs in minimizers.items()
        },
    }


def scale_split_controls() -> dict[str, str]:
    first_scale_surplus = 3 * BODY_ARC_LENGTH - Q(2, 63)
    large_ratio_surplus = BODY_ARC_LENGTH - Q(2, 7 * 27)
    require(first_scale_surplus == Q(1, 2520) > 0,
            "first t>=3 compact/open surplus")
    require(large_ratio_surplus == Q(1, 7560) > 0,
            "first q>=27 compact/open surplus")
    return {
        "t_ge_3_first_surplus": fmt(first_scale_surplus),
        "t_1_q_ge_27_first_surplus": fmt(large_ratio_surplus),
        "strict_rule": "compact circular arc containment in an open component forces strict length inequality",
    }


def main() -> None:
    semantic = {
        "universe": {
            "body": list(BODY),
            "tail_scope": "two distinct positive odd integers",
            "normalization": "a=pt,b=qt,t=gcd(a,b),p<q odd coprime",
            "danger": "strict ||v*x||<1/14; equality is safe",
        },
        "body_arc": body_arc_controls(),
        "low_widths": low_width_controls(),
        "q_tooth_bound": q_tooth_bound_controls(),
        "residual_clocks": residual_clock_controls(),
        "scale_split": scale_split_controls(),
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    clock_data = semantic["residual_clocks"]
    width_data = semantic["q_tooth_bound"]
    low_data = semantic["low_widths"]
    scale_data = semantic["scale_split"]

    print("LRC14_FIXED_BODY_UNIVERSAL_ODD_TAIL_THM4136_INDEPENDENT_AUDIT")
    print("scope=U_fixed;two_distinct_positive_odd_tails;LRC14_global=OPEN")
    print("implementation=clean_room_stdlib_only_exact_Fraction_no_project_imports")
    print("body=1,4,6,8,10,12,14,15,16,18,22")
    print("closed_body_arc=[33/70,27/56];length=3/280;owners=15,4")
    print("low_quotient_components=" + repr(low_data["low_components"]))
    print("exceptional_1_9_widths=" + repr(tuple(low_data["exceptional_1_9"])))
    print("q_tooth_bound=beta_pq<=2/(7q)_for_odd_q>=9")
    print(
        f"q_tooth_scan=primitive_odd_q<=101:{width_data['exhaustive_rows']};"
        f"large_hostiles_to_q={width_data['large_q_max']}:{width_data['large_hostile_rows']};"
        f"bound_equalities={width_data['bound_equality_count']};"
        f"global_max={width_data['global_max']};"
        f"global_max_rows={tuple(tuple(row) for row in width_data['global_max_rows'])}"
    )
    print("q_tooth_scan_transcript_sha256=" + width_data["scan_transcript_sha256"])
    print(
        "scale_surpluses="
        f"t>=3_first:{scale_data['t_ge_3_first_surplus']};"
        f"t=1,q>=27_first:{scale_data['t_1_q_ge_27_first_surplus']}"
    )
    print(
        f"residual_pairs={len(clock_data['residual_pairs'])};"
        "clock_split=56,10,2;"
        f"clocks={clock_data['clocks']}"
    )
    print("clock_minima=" + repr(clock_data["minima"]))
    print("consequence_check=all_13_literal_clearances_verified_in_each_residual_row")
    print("hostiles=strict_endpoints;active_wall_punctures;unbalanced_and_near_diagonal_large_q")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=ACCEPT;THM4136_CANDIDATE_VERIFIED;LRC14=OPEN")


if __name__ == "__main__":
    main()
