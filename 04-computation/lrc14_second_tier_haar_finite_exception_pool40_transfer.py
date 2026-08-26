#!/usr/bin/env python3
"""Exact primary audit for the second-tier Haar/finite-exception transfer.

Candidate theorem (identifier deliberately not reserved): a full-safe-set
measure of at least 4/77 excludes every primitive odd cross-comb except
(1,9).  A largest safe component reduces the remaining common scales to a
finite bank.  The explicit 40-label pool below closes that bank at one clock.
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
FIRST_TIER = Q(4, 63)
SECOND_TIER = Q(4, 77)
COMPONENT_CAP_19 = Q(2, 63)
WIDTH_GATE = Q(2, 189)
POOL33 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 69, 71, 73, 75, 76, 80,
)
POOL40 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 67, 69, 71, 73,
    75, 76, 80, 82, 86, 89, 93, 95, 141,
)
MAX_COMPONENT = (Q(299, 658), Q(321, 700))
CLOCK_BODY_PHASE = Q(57, 742)
CLOCK_PHASE = Q(799, 1484)
FINITE_SCALES = (1, 3, 5, 7)
EXPECTED_SEMANTIC = "f13989721e4bddfbfea0b48abdaba9201e06daed15e7cbf7ad139084ce00168c"

CHECKS = 0


def require(predicate: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def semantic_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(circle_distance(speed * phase) for speed in speeds)


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


def both_sheets_bad(p: int, q: int, w: Q) -> bool:
    z = w / 2
    return (
        (speed_bad(p, z) or speed_bad(q, z))
        and (speed_bad(p, z + Q(1, 2)) or speed_bad(q, z + Q(1, 2)))
    )


def cross_wall_geometry(p: int, q: int) -> tuple[Q, int, Q]:
    walls = {Q(0), Q(1)}
    for shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    z_wall = (Q(integer) + sign * DELTA) / speed - shift
                    walls.add((2 * z_wall) % 1)
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
    width = max((right - left for left, right in components), default=Q(0))
    return total, len(components), width


def second_tier_audit() -> dict[str, object]:
    # Since B_2 has oscillation 1/4, mu(C)<=2/49+1/(2pq).
    # Product 45 is the first automatic value for the 4/77 target.
    require(Q(1, 2 * 45) < SECOND_TIER - Q(2, 49), "product-45 gate")
    small_pairs = tuple(
        (p, q)
        for q in range(3, 45, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1 and p * q < 45
    )
    require(len(small_pairs) == 26, "finite product universe")
    above = tuple(pair_ for pair_ in small_pairs if cross_measure_formula(*pair_) > SECOND_TIER)
    equal = tuple(pair_ for pair_ in small_pairs if cross_measure_formula(*pair_) == SECOND_TIER)
    require(above == ((1, 9),), "unique exceptional primitive ratio")
    require(equal == ((1, 11), (3, 11)), "second-tier equality ratios")
    require(cross_measure_formula(1, 9) == FIRST_TIER, "exception measure")

    scanned = 0
    maximum_nonexception = Q(0)
    maximizers: list[tuple[int, int]] = []
    for q in range(3, 102, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scanned += 1
            formula = cross_measure_formula(p, q)
            wall_measure, _, _ = cross_wall_geometry(p, q)
            require(formula == wall_measure, ("formula/wall", p, q))
            if (p, q) == (1, 9):
                continue
            require(formula <= SECOND_TIER, ("second-tier cap", p, q))
            if formula > maximum_nonexception:
                maximum_nonexception = formula
                maximizers = [(p, q)]
            elif formula == maximum_nonexception:
                maximizers.append((p, q))
    require(scanned == 1053, "scan size")
    require(maximum_nonexception == SECOND_TIER, "sharp nonexception cap")
    require(maximizers == [(1, 11), (3, 11)], "sharp nonexception ratios")
    return {
        "small_pair_count": len(small_pairs),
        "above": above,
        "equal": equal,
        "scan": (scanned, 101, pair(maximum_nonexception), tuple(maximizers)),
    }


def pool_audit() -> dict[str, object]:
    require(len(POOL40) == len(set(POOL40)) == 40, "pool40 cardinality")
    require(set(POOL33) < set(POOL40), "strict P33 inheritance")
    require(tuple(value for value in POOL40 if value not in POOL33) ==
            (67, 82, 86, 89, 93, 95, 141), "seven new labels")

    walls = safe_walls(POOL40)
    components = safe_components(POOL40)
    measure = sum((right - left for left, right in components), Q(0))
    largest = max(right - left for left, right in components)
    maximizers = tuple(component for component in components if component[1] - component[0] == largest)
    require(len(walls) == 3744, "pool wall count")
    require(len(components) == 44, "pool component count")
    require(measure == Q(23518182747542658511201, 441420293060220631236800),
            "pool measure")
    require(measure - SECOND_TIER ==
            Q(6459842759986025773611, 4855623223662426943604800),
            "pool second-tier surplus")
    require(measure < FIRST_TIER, "first-tier hostile")
    require(largest == Q(137, 32900), "largest component width")
    require(maximizers == (MAX_COMPONENT, (Q(379, 700), Q(359, 658))),
            "largest reflection pair")
    require(9 * largest - COMPONENT_CAP_19 == Q(1697, 296100) > 0,
            "scale-nine escape")
    require(7 * largest < COMPONENT_CAP_19, "scale-seven width hostile")
    finite_scales = tuple(t for t in range(1, 100, 2) if t * largest < COMPONENT_CAP_19)
    require(finite_scales == FINITE_SCALES, "finite exceptional scale bank")

    require(CLOCK_PHASE == (CLOCK_BODY_PHASE + 1) / 2, "physical upper lift")
    body_gap = clearance(tuple(2 * speed for speed in POOL40), CLOCK_PHASE)
    require(body_gap == DELTA, "clock body clearance")
    require(tuple(speed for speed in POOL40
                  if circle_distance(2 * speed * CLOCK_PHASE) == DELTA) == (53,),
            "clock body owner")
    tail_gaps = []
    for scale in FINITE_SCALES:
        gaps = (
            circle_distance(scale * CLOCK_PHASE),
            circle_distance(9 * scale * CLOCK_PHASE),
        )
        require(min(gaps) > DELTA, ("finite scale clock", scale, gaps))
        require(clearance(tuple(2 * speed for speed in POOL40) + (scale, 9 * scale),
                          CLOCK_PHASE) == DELTA,
                ("complete finite row", scale))
        tail_gaps.append((scale, pair(gaps[0]), pair(gaps[1])))

    total = comb(len(POOL40), 11)
    old_total = comb(len(POOL33), 11)
    width_covered = 0
    for i, minimum in enumerate(POOL40):
        for j in range(i + 10, len(POOL40)):
            maximum = POOL40[j]
            if 27 * (13 * minimum - maximum) >= 4 * minimum * maximum:
                width_covered += comb(j - i - 1, 9)
    require(total == 2_311_801_440, "pool family total")
    require(old_total == 193_536_720, "P33 inherited family total")
    require(total - old_total == 2_118_264_720, "strict P33 increment")
    require(width_covered == 253_100, "current width-gate overlap")
    require(total - width_covered == 2_311_548_340, "beyond-width families")

    return {
        "pool": POOL40,
        "geometry": (len(walls), len(components), pair(measure), pair(largest)),
        "surplus": pair(measure - SECOND_TIER),
        "maximizers": tuple((pair(left), pair(right)) for left, right in maximizers),
        "finite_scales": finite_scales,
        "clock": (pair(CLOCK_BODY_PHASE), pair(CLOCK_PHASE), pair(body_gap), tuple(tail_gaps)),
        "families": (total, old_total, total - old_total, width_covered, total - width_covered),
    }


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    require(b2(0) - b2(Q(1, 2)) == Q(1, 4), "Bernoulli oscillation")
    tier = second_tier_audit()
    pool = pool_audit()
    semantic = {
        "delta": pair(DELTA),
        "second_tier": pair(SECOND_TIER),
        "component_cap": pair(COMPONENT_CAP_19),
        "tier": tier,
        "pool": pool,
    }
    semantic_hash = semantic_digest(semantic)
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        require(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    print("LRC14_SECOND_TIER_HAAR_FINITE_EXCEPTION_POOL40_20260825")
    print("status=PASS;scope=candidate exact theorem audit;LRC14=OPEN")
    print("second_tier=4/77;only_above=(1,9);equal=((1,11),(3,11))")
    print(f"cross_scan={tier['scan']}")
    print(f"pool_size={len(POOL40)};pool={POOL40}")
    print(f"pool_geometry={pool['geometry']};surplus={pool['surplus']}")
    print(f"max_components={pool['maximizers']}")
    print(f"finite_scales={pool['finite_scales']};clock={pool['clock']}")
    print(f"families={pool['families'][0]};P33={pool['families'][1]};increment={pool['families'][2]}")
    print(f"width_covered={pool['families'][3]};beyond_width={pool['families'][4]}")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
