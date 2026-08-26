#!/usr/bin/env python3
"""Exact primary audit for THM-4150.

This script keeps the two proof coordinates separate:

* the cross-comb measure is reconstructed both from strict rational wall
  cells and from the Bernoulli-polynomial formula in the proof;
* body-safe Haar measure is reconstructed from the complete rational wall
  arrangement, without selecting a preferred component or phase.

All arithmetic affecting a stated value is exact ``Fraction`` arithmetic.
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
THRESHOLD = Q(4, 63)
WIDTH_GATE = Q(2, 189)
POOL = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 69, 71, 73, 75, 76, 80,
)
FRAGMENTED_BODY = (1, 17, 31, 32, 41, 47, 50, 51, 58, 62, 71)
THM4142_POOL = {
    1, 3, 4, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18,
    22, 24, 28, 30, 32, 35, 37, 41, 43, 45, 49, 60, 64,
}
EXPECTED_SEMANTIC = "9b08f96a0b6771ddadc366232f0cd7e3d07b208d3328668890b6a6cc24a5390f"

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def digest(value: object) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(data).hexdigest()


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
    """Positive-length components of the closed body-safe set."""
    walls = safe_walls(speeds)
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(walls, walls[1:]):
        middle = (left + right) / 2
        if clearance(speeds, middle) < DELTA:
            continue
        require(clearance(speeds, left) >= DELTA, ("safe left", left, right))
        require(clearance(speeds, right) >= DELTA, ("safe right", left, right))
        cells.append((left, right))
    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left:
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return tuple(components)


def bernoulli_two(value: Q) -> Q:
    value %= 1
    return value * value - value + Q(1, 6)


def cross_measure_formula(p: int, q: int) -> Q:
    """Haar measure of C_(p,q), using the proved Fourier/Bernoulli formula."""
    require(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
            ("primitive odd ratio", p, q))
    minus = (Q(1, 2) + Q(q - p, 14)) % 1
    plus = (Q(1, 2) + Q(q + p, 14)) % 1
    b_measure = Q(1, 49) + (bernoulli_two(minus) - bernoulli_two(plus)) / (p * q)
    return 2 * b_measure


def speed_bad(speed: int, phase: Q) -> bool:
    return circle_distance(speed * phase) < DELTA


def both_sheets_bad(p: int, q: int, w: Q) -> bool:
    z = w / 2
    return (
        (speed_bad(p, z) or speed_bad(q, z))
        and (speed_bad(p, z + Q(1, 2)) or speed_bad(q, z + Q(1, 2)))
    )


def cross_measure_walls(p: int, q: int) -> tuple[Q, int, Q]:
    """Independent-in-script strict-wall reconstruction of C_(p,q)."""
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
    measure = sum((right - left for left, right in components), Q(0))
    width = max((right - left for left, right in components), default=Q(0))
    return measure, len(components), width


def cross_comb_audit() -> dict[str, object]:
    small_pairs = tuple(
        (p, q)
        for q in range(3, 23, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1 and p * q <= 22
    )
    expected_small = {
        (1, 3): Q(0),
        (1, 5): Q(0),
        (1, 7): Q(2, 49),
        (1, 9): Q(4, 63),
        (1, 11): Q(4, 77),
        (1, 13): Q(4, 91),
        (1, 15): Q(4, 105),
        (1, 17): Q(4, 119),
        (1, 19): Q(4, 133),
        (1, 21): Q(2, 49),
        (3, 5): Q(2, 105),
        (3, 7): Q(2, 49),
    }
    require(set(small_pairs) == set(expected_small), "complete pq<=22 universe")
    for ratio, expected in expected_small.items():
        require(cross_measure_formula(*ratio) == expected, ("small formula", ratio))
        require(cross_measure_walls(*ratio)[0] == expected, ("small walls", ratio))
    maximizers = tuple(ratio for ratio, value in expected_small.items() if value == THRESHOLD)
    require(maximizers == ((1, 9),), "unique small equality")

    # B_2 takes values in [-1/12,1/6], hence the correction in mu(C)
    # is at most 1/(2pq).  At pq>=23 this is strictly below the equality
    # correction 10/441.
    require(Q(1, 2 * 23) < THRESHOLD - Q(2, 49), "large-product analytic gate")

    scanned = 0
    scan_max = Q(0)
    scan_maximizers: list[tuple[int, int]] = []
    for q in range(3, 102, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scanned += 1
            formula = cross_measure_formula(p, q)
            walls, _, _ = cross_measure_walls(p, q)
            require(formula == walls, ("formula/wall agreement", p, q))
            require(formula <= THRESHOLD, ("sharp measure bound", p, q))
            if formula > scan_max:
                scan_max = formula
                scan_maximizers = [(p, q)]
            elif formula == scan_max:
                scan_maximizers.append((p, q))
    require(scanned == 1053, "hostile scan count")
    require(scan_max == THRESHOLD and scan_maximizers == [(1, 9)],
            "hostile scan unique maximum")
    return {
        "small": tuple((ratio, pair(value)) for ratio, value in expected_small.items()),
        "large_product_start": 23,
        "scan": (scanned, 101, pair(scan_max), tuple(scan_maximizers)),
    }


def body_audit() -> dict[str, object]:
    require(len(POOL) == len(set(POOL)) == 33, "pool cardinality")
    components = safe_components(POOL)
    measure = sum((right - left for left, right in components), Q(0))
    maximum_width = max(right - left for left, right in components)
    require(len(safe_walls(POOL)) == 2472, "pool wall count")
    require(len(components) == 46, "pool component count")
    require(measure == Q(110551382435042260737, 1702610555154297252800),
            "pool safe measure")
    require(measure - THRESHOLD ==
            Q(3148874954907620719, 2189070713769810753600),
            "pool measure surplus")
    require(maximum_width == Q(3, 700) < WIDTH_GATE, "pool fragmented width hostile")
    require(tuple((1 - right, 1 - left) for left, right in reversed(components)) == components,
            "pool reflection symmetry")

    family_count = comb(len(POOL), 11)
    require(family_count == 193_536_720, "pool eleven-body count")
    require(len(set(POOL) & THM4142_POOL) == 9, "not enough labels for THM4142 family")

    width_covered = 0
    for i, minimum in enumerate(POOL):
        for j in range(i + 10, len(POOL)):
            maximum = POOL[j]
            choices = comb(j - i - 1, 9)
            if 27 * (13 * minimum - maximum) >= 4 * minimum * maximum:
                width_covered += choices
    require(width_covered == 208_000, "THM4148-covered subfamily")
    new_beyond_width = family_count - width_covered
    require(new_beyond_width == 193_328_720, "new beyond-width subfamily")

    fragmented = safe_components(FRAGMENTED_BODY)
    fragmented_measure = sum((right - left for left, right in fragmented), Q(0))
    fragmented_width = max(right - left for left, right in fragmented)
    require(len(fragmented) == 100, "fragmented body components")
    require(fragmented_measure == Q(1142176622583, 5854727790800),
            "fragmented body measure")
    require(fragmented_width == Q(21, 4100) < WIDTH_GATE,
            "fragmented body has no wide component")
    require(Q(13, 14 * max(FRAGMENTED_BODY)) - Q(1, 14) < 0,
            "fragmented body first window empty")

    # A below-threshold safe set is not a negative certificate.  This control
    # has measure <4/63, yet the particular completed row is visibly safe.
    low_body = tuple(range(1, 12))
    low_components = safe_components(low_body)
    low_measure = sum((right - left for left, right in low_components), Q(0))
    require(low_measure == Q(10931, 194040) < THRESHOLD, "criterion is sufficient only")
    require(clearance(tuple(2 * h for h in low_body) + (1, 3), Q(1, 13)) == Q(1, 13),
            "below-threshold positive row")

    return {
        "pool": POOL,
        "pool_geometry": (
            len(safe_walls(POOL)), len(components), pair(measure), pair(maximum_width)
        ),
        "surplus": pair(measure - THRESHOLD),
        "families": (family_count, width_covered, new_beyond_width),
        "fragmented_body": (
            FRAGMENTED_BODY, len(fragmented), pair(fragmented_measure), pair(fragmented_width)
        ),
        "below_threshold_control": (pair(low_measure), pair(Q(1, 13))),
    }


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    require(bernoulli_two(Q(0)) == Q(1, 6), "B2 maximum")
    require(bernoulli_two(Q(1, 2)) == -Q(1, 12), "B2 minimum")
    cross = cross_comb_audit()
    body = body_audit()

    semantic = {
        "delta": pair(DELTA),
        "threshold": pair(THRESHOLD),
        "cross": cross,
        "body": body,
    }
    semantic_hash = digest(semantic)
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        require(semantic_hash == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("THM4150_SAFE_SET_HAAR_ODD_TAIL_TRANSFER_20260825")
    print("status=PASS;scope=exact cross-comb measure and body-safe-set audit")
    print("cross_measure_bound=4/63;unique_equality=(1,9);large_product_start=23")
    print(f"cross_scan={cross['scan']}")
    print(f"pool_size={len(POOL)};pool={POOL}")
    print(f"pool_wall_count={body['pool_geometry'][0]};pool_components={body['pool_geometry'][1]}")
    print(f"pool_measure={body['pool_geometry'][2][0]}/{body['pool_geometry'][2][1]}")
    print(f"pool_surplus={body['surplus'][0]}/{body['surplus'][1]}")
    print(f"pool_max_component={body['pool_geometry'][3][0]}/{body['pool_geometry'][3][1]}")
    print(f"families={body['families'][0]};width_covered={body['families'][1]};beyond_width={body['families'][2]}")
    print(f"fragmented_body={FRAGMENTED_BODY}")
    print(f"fragmented_geometry=(components={body['fragmented_body'][1]},measure={body['fragmented_body'][2][0]}/{body['fragmented_body'][2][1]},max_width={body['fragmented_body'][3][0]}/{body['fragmented_body'][3][1]})")
    print(f"below_threshold_control={body['below_threshold_control']}")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
