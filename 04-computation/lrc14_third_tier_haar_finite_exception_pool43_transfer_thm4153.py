#!/usr/bin/env python3
"""Exact primary audit for THM-4153's 4/91 Haar/clock transfer.

The proof has three independent finite objects:

* the Bernoulli cross-comb formula and the complete pq<160 cap census;
* the positive-length safe components of the explicit 43-label pool; and
* the residual primitive-ratio/odd-scale bank discharged by one clock.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
THRESHOLD = Q(4, 91)
POOL40 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 67, 69, 71, 73,
    75, 76, 80, 82, 86, 89, 93, 95, 141,
)
POOL43 = tuple(sorted(POOL40 + (111, 159, 285)))
BODY_PHASE = Q(57, 742)
PHYSICAL_CLOCK = Q(799, 1484)
COMPONENT_WIDTH = Q(1, 399)
EXPECTED_SEMANTIC = "8fc29d9da1111f0b2cf5e13c6cd6dbca6e68072acd7f4b95f5455f7e412245f7"

# (primitive ratio): (cross measure, largest open-component width, residual
# odd scales for a safe component of length 1/399).
EXCEPTION_DATA = {
    (1, 9): (Q(4, 63), Q(2, 63), (1, 3, 5, 7, 9, 11)),
    (1, 11): (Q(4, 77), Q(2, 77), (1, 3, 5, 7, 9)),
    (3, 11): (Q(4, 77), Q(2, 77), (1, 3, 5, 7, 9)),
    (1, 23): (Q(8, 161), Q(2, 161), (1, 3)),
    (5, 11): (Q(18, 385), Q(9, 385), (1, 3, 5, 7, 9)),
    (1, 37): (Q(12, 259), Q(2, 259), (1, 3)),
    (1, 25): (Q(8, 175), Q(2, 175), (1, 3)),
    (3, 25): (Q(8, 175), Q(2, 175), (1, 3)),
    (3, 23): (Q(22, 483), Q(2, 161), (1, 3)),
    (1, 51): (Q(16, 357), Q(2, 357), (1,)),
    (5, 9): (Q(2, 45), Q(1, 45), (1, 3, 5, 7)),
}
EXPECTED_ABOVE = tuple(sorted(EXCEPTION_DATA))
EXPECTED_EQUAL = ((1, 13), (1, 39), (1, 65), (3, 13), (5, 13))

CHECKS = 0


def require(predicate: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def digest(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(raw).hexdigest()


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(circle_distance(speed * phase) for speed in speeds)


def b2(value: Q) -> Q:
    value %= 1
    return value * value - value + Q(1, 6)


def cross_measure_formula(p: int, q: int) -> Q:
    require(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
            ("primitive ratio", p, q))
    u_minus = (Q(1, 2) + Q(q - p, 14)) % 1
    u_plus = (Q(1, 2) + Q(q + p, 14)) % 1
    return 2 * (Q(1, 49) + (b2(u_minus) - b2(u_plus)) / (p * q))


def speed_bad(speed: int, phase: Q) -> bool:
    return circle_distance(speed * phase) < DELTA


def both_sheets_bad(p: int, q: int, quotient_phase: Q) -> bool:
    lower = quotient_phase / 2
    upper = lower + Q(1, 2)
    return (
        (speed_bad(p, lower) or speed_bad(q, lower))
        and (speed_bad(p, upper) or speed_bad(q, upper))
    )


def cross_geometry(p: int, q: int) -> tuple[Q, int, Q]:
    walls = {Q(0), Q(1)}
    for shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    lifted_wall = (Q(integer) + sign * DELTA) / speed - shift
                    walls.add((2 * lifted_wall) % 1)
    ordered = tuple(sorted(walls))
    cells = [
        (left, right)
        for left, right in zip(ordered, ordered[1:])
        if both_sheets_bad(p, q, (left + right) / 2)
    ]
    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left and both_sheets_bad(p, q, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    total = sum((right - left for left, right in components), Q(0))
    maximum = max((right - left for left, right in components), default=Q(0))
    return total, len(components), maximum


def safe_walls(speeds: tuple[int, ...]) -> tuple[Q, ...]:
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for integer in range(speed):
            walls.add((Q(integer) + DELTA) / speed)
            walls.add((Q(integer + 1) - DELTA) / speed)
    return tuple(sorted(walls))


def safe_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    walls = safe_walls(speeds)
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(walls, walls[1:]):
        middle = (left + right) / 2
        if clearance(speeds, middle) < DELTA:
            continue
        require(clearance(speeds, left) >= DELTA, ("safe left", left, right))
        require(clearance(speeds, right) >= DELTA, ("safe right", left, right))
        cells.append((left, right))
    result: list[tuple[Q, Q]] = []
    for left, right in cells:
        if result and result[-1][1] == left:
            result[-1] = (result[-1][0], right)
        else:
            result.append((left, right))
    return tuple(result)


def threshold_audit() -> dict[str, object]:
    # Osc(B_2)=1/4 gives mu(C_pq)<=2/49+1/(2pq).  Since
    # 4/91-2/49=2/637, every odd product pq>=161 is automatic.
    require(THRESHOLD - Q(2, 49) == Q(2, 637), "threshold surplus")
    require(Q(1, 2 * 161) < Q(2, 637), "product-161 automatic gate")

    small_pairs = tuple(
        (p, q)
        for q in range(3, 160, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1 and p * q < 160
    )
    above = tuple(sorted(
        ratio for ratio in small_pairs
        if cross_measure_formula(*ratio) > THRESHOLD
    ))
    equal = tuple(sorted(
        ratio for ratio in small_pairs
        if cross_measure_formula(*ratio) == THRESHOLD
    ))
    require(len(small_pairs) == 117, "finite product census")
    require(above == EXPECTED_ABOVE, "complete above-threshold list")
    require(equal == EXPECTED_EQUAL, "complete equality list")

    scanned = 0
    maximum_nonexception = Q(0)
    maximizers: list[tuple[int, int]] = []
    for q in range(3, 102, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scanned += 1
            formula = cross_measure_formula(p, q)
            direct, _, _ = cross_geometry(p, q)
            require(formula == direct, ("formula/direct wall measure", p, q))
            if (p, q) in EXCEPTION_DATA:
                continue
            require(formula <= THRESHOLD, ("third-tier cap", p, q))
            if formula > maximum_nonexception:
                maximum_nonexception = formula
                maximizers = [(p, q)]
            elif formula == maximum_nonexception:
                maximizers.append((p, q))
    require(scanned == 1053, "direct cross scan")
    require(maximum_nonexception == THRESHOLD, "sharp third-tier cap")
    require(tuple(sorted(maximizers)) == EXPECTED_EQUAL,
            "sharp equality ratios")

    geometry = []
    for ratio in EXPECTED_ABOVE:
        expected_measure, expected_width, expected_scales = EXCEPTION_DATA[ratio]
        measure, count, width = cross_geometry(*ratio)
        scales = tuple(
            scale for scale in range(1, 100, 2)
            if scale * COMPONENT_WIDTH < width
        )
        require(measure == expected_measure, ("exception measure", ratio))
        require(width == expected_width, ("exception width", ratio))
        require(scales == expected_scales, ("exception scales", ratio))
        geometry.append((ratio, pair(measure), count, pair(width), scales))
    return {
        "finite_product_count": len(small_pairs),
        "above": above,
        "equal": equal,
        "scan": (scanned, pair(maximum_nonexception), tuple(sorted(maximizers))),
        "geometry": tuple(geometry),
    }


def pool_audit() -> dict[str, object]:
    require(len(POOL43) == len(set(POOL43)) == 43, "pool cardinality")
    require(set(POOL40) < set(POOL43), "strict P40 inheritance")
    require(tuple(value for value in POOL43 if value not in POOL40) ==
            (111, 159, 285), "three new labels")

    walls = safe_walls(POOL43)
    components = safe_components(POOL43)
    measure = sum((right - left for left, right in components), Q(0))
    largest = max(right - left for left, right in components)
    maximizers = tuple(
        component for component in components
        if component[1] - component[0] == largest
    )
    require(len(walls) == 4776, "pool wall count")
    require(len(components) == 46, "pool component count")
    require(measure ==
            Q(10080921463555906580413, 211196778145191767531400),
            "pool Haar measure")
    require(measure - THRESHOLD ==
            Q(10368105800402918384569, 2745558115887492977908200),
            "pool threshold surplus")
    require(largest == COMPONENT_WIDTH, "largest component width")
    require(maximizers ==
            ((Q(911, 3990), Q(307, 1330)),
             (Q(1023, 1330), Q(3079, 3990))),
            "largest reflection pair")

    require(PHYSICAL_CLOCK == (BODY_PHASE + 1) / 2, "upper physical lift")
    body_gap = clearance(POOL43, BODY_PHASE)
    physical_body_gap = clearance(tuple(2 * value for value in POOL43),
                                  PHYSICAL_CLOCK)
    require(body_gap == physical_body_gap == DELTA, "clock body clearance")
    require(tuple(value for value in POOL43
                  if circle_distance(value * BODY_PHASE) == DELTA) == (53,),
            "unique body owner")

    clock_rows = []
    global_tail_gap = Q(1, 2)
    for ratio in EXPECTED_ABOVE:
        p, q = ratio
        rows = []
        for scale in EXCEPTION_DATA[ratio][2]:
            gaps = (
                circle_distance(p * scale * PHYSICAL_CLOCK),
                circle_distance(q * scale * PHYSICAL_CLOCK),
            )
            require(min(gaps) > DELTA, ("clock tail gap", ratio, scale, gaps))
            require(clearance(
                tuple(2 * value for value in POOL43) + (p * scale, q * scale),
                PHYSICAL_CLOCK,
            ) == DELTA, ("full clock row", ratio, scale))
            global_tail_gap = min(global_tail_gap, *gaps)
            rows.append((scale, pair(gaps[0]), pair(gaps[1])))
        clock_rows.append((ratio, tuple(rows)))
    require(global_tail_gap == Q(113, 1484), "sharp finite-bank tail gap")
    require(global_tail_gap - DELTA == Q(1, 212), "clock strict surplus")

    hostile_measure, _, hostile_width = cross_geometry(5, 13)
    hostile_scales = tuple(
        scale for scale in range(1, 100, 2)
        if scale * COMPONENT_WIDTH < hostile_width
    )
    hostile_gaps = tuple(
        (
            scale,
            pair(circle_distance(5 * scale * PHYSICAL_CLOCK)),
            pair(circle_distance(13 * scale * PHYSICAL_CLOCK)),
        )
        for scale in hostile_scales
    )
    require(hostile_measure == THRESHOLD, "next-tier hostile measure")
    require(hostile_width == Q(2, 91), "next-tier hostile width")
    require(hostile_scales == (1, 3, 5, 7), "next-tier hostile scale bank")
    require(all(
        min(Q(*first), Q(*second)) < DELTA
        for _, first, second in hostile_gaps
    ), "current clock fails every hostile scale")

    total = comb(43, 11)
    inherited = comb(40, 11)
    width_covered = 0
    for i, minimum in enumerate(POOL43):
        for j in range(i + 10, len(POOL43)):
            maximum = POOL43[j]
            if 27 * (13 * minimum - maximum) >= 4 * minimum * maximum:
                width_covered += comb(j - i - 1, 9)
    require(total == 5_752_004_349, "family count")
    require(inherited == 2_311_801_440, "P40 inherited count")
    require(total - inherited == 3_440_202_909, "strict P40 increment")
    require(width_covered == 253_100, "width-gate overlap")
    require(total - width_covered == 5_751_751_249, "beyond-width families")

    return {
        "pool": POOL43,
        "geometry": (
            len(walls), len(components), pair(measure), pair(largest),
            tuple((pair(left), pair(right)) for left, right in maximizers),
        ),
        "surplus": pair(measure - THRESHOLD),
        "clock": (
            pair(BODY_PHASE), pair(PHYSICAL_CLOCK), pair(body_gap),
            pair(global_tail_gap), tuple(clock_rows),
        ),
        "next_hostile": (
            (5, 13), pair(hostile_measure), pair(hostile_width),
            hostile_scales, hostile_gaps,
        ),
        "families": (
            total, inherited, total - inherited,
            width_covered, total - width_covered,
        ),
    }


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "no assertions")
    require(b2(0) - b2(Q(1, 2)) == Q(1, 4), "Bernoulli oscillation")
    threshold = threshold_audit()
    pool = pool_audit()
    semantic = {
        "delta": pair(DELTA),
        "threshold": pair(THRESHOLD),
        "component_width": pair(COMPONENT_WIDTH),
        "threshold_audit": threshold,
        "pool_audit": pool,
    }
    semantic_hash = digest(semantic)
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        require(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    print("LRC14_THIRD_TIER_HAAR_FINITE_EXCEPTION_POOL43_THM4153_20260825")
    print("status=PASS;scope=exact theorem audit;LRC14=OPEN")
    print(f"threshold=4/91;above={threshold['above']};equal={threshold['equal']}")
    print(f"cross_scan={threshold['scan']};finite_product_count={threshold['finite_product_count']}")
    print(f"exception_geometry={threshold['geometry']}")
    print(f"pool_size={len(POOL43)};pool={POOL43}")
    print(f"pool_geometry={pool['geometry']};surplus={pool['surplus']}")
    print(f"clock_summary={pool['clock'][:4]};clock_rows={pool['clock'][4]}")
    print(f"next_hostile={pool['next_hostile']}")
    print(f"families={pool['families'][0]};P40={pool['families'][1]};increment={pool['families'][2]}")
    print(f"width_covered={pool['families'][3]};beyond_width={pool['families'][4]}")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
