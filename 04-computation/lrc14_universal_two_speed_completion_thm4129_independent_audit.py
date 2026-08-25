#!/usr/bin/env python3
"""Clean-room cumulative-measure and lattice audit for THM-4129."""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import ceil, floor, gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
BODY_SET = frozenset(BODY)
DELTA = Fraction(1, 14)
WINDOW = (Fraction(33, 70), Fraction(27, 56))
LOW_MISSING = (2, 3, 5, 7, 9, 11, 13, 17, 19, 20, 21)
LIVE_SCALE_ONE_TYPES = (
    (1, 3), (1, 4), (1, 9), (1, 10), (2, 3), (2, 9),
    (3, 7), (3, 8), (4, 7), (5, 6), (5, 12), (6, 11),
    (7, 10), (8, 9), (8, 21), (9, 11),
)
EXPECTED_SEMANTIC = "a717248ff237c9cbadf1a15dbaeb2902c0b7b7da10ebca0e29d797ea09bbdbe5"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def rat(value: Fraction) -> tuple[int, int]:
    return (value.numerator, value.denominator)


def semantic_digest(value: object) -> str:
    raw = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(raw).hexdigest()


def fractional_part(value: Fraction) -> Fraction:
    return value - floor(value)


def torus_distance(value: Fraction) -> Fraction:
    residue = fractional_part(value)
    return min(residue, 1 - residue)


def row_clearance(row: tuple[int, ...], phase: Fraction) -> Fraction:
    return min(torus_distance(speed * phase) for speed in row)


def cumulative_danger(real: Fraction) -> Fraction:
    """Measure of {y in [0,real]: ||y||<1/14} for real >= 0."""
    integer = floor(real)
    residue = real - integer
    if residue <= DELTA:
        tail = residue
    elif residue <= 1 - DELTA:
        tail = DELTA
    else:
        tail = residue - Fraction(6, 7)
    return Fraction(integer, 7) + tail


def antiderivative_measure(speed: int, interval: tuple[Fraction, Fraction]) -> Fraction:
    left, right = interval
    return (cumulative_danger(speed * right) - cumulative_danger(speed * left)) / speed


def endpoint_certificate(
    row: tuple[int, ...], interval: tuple[Fraction, Fraction]
) -> tuple[int, int]:
    left, right = interval
    require(left < right, ("positive interval", interval))
    left_values = []
    right_values = []
    for speed in row:
        require((speed * left).numerator // (speed * left).denominator
                == (speed * right).numerator // (speed * right).denominator,
                ("no integer crossing", speed, interval))
        left_values.append((torus_distance(speed * left), speed))
        right_values.append((torus_distance(speed * right), speed))
    left_min = min(left_values)
    right_min = min(right_values)
    require(left_min[0] >= DELTA and right_min[0] >= DELTA,
            ("closed endpoint safety", row, interval))
    require(left_min[0] == DELTA and right_min[0] == DELTA,
            ("sharp closed endpoints", row, interval))
    left_owners = tuple(speed for speed in row
                        if torus_distance(speed * left) == DELTA)
    right_owners = tuple(speed for speed in row
                         if torus_distance(speed * right) == DELTA)
    require(len(left_owners) == len(right_owners) == 1,
            ("unique closed endpoint owners", row, interval))
    return (left_owners[0], right_owners[0])


def common_grid_witness(first: int, second: int) -> Fraction | None:
    """Search the exact grid containing every wall endpoint in WINDOW."""
    denominator = 280 * first * second
    left_scaled = WINDOW[0] * denominator
    right_scaled = WINDOW[1] * denominator
    lower = ceil(left_scaled)
    upper = floor(right_scaled)
    for numerator in range(lower, upper + 1):
        phase = Fraction(numerator, denominator)
        if (torus_distance(first * phase) >= DELTA
                and torus_distance(second * phase) >= DELTA):
            return phase
    return None


def smallest_prime_factors(limit: int) -> list[int]:
    factors = list(range(limit + 1))
    for prime in range(2, int(limit ** 0.5) + 1):
        if factors[prime] == prime:
            for multiple in range(prime * prime, limit + 1, prime):
                if factors[multiple] == multiple:
                    factors[multiple] = prime
    return factors


def admissible_sum_sieve(value: int, factors: list[int]) -> bool:
    if value <= 1:
        return False
    while value > 1:
        prime = factors[value]
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        if prime % 3 != 2 or exponent > 2:
            return False
    return True


def atlas_sieve() -> tuple[tuple[int, int], ...]:
    factors = smallest_prime_factors(356)
    result = []
    for total in range(3, 357):
        if not admissible_sum_sieve(total, factors):
            continue
        for first in range(1, (total + 1) // 2):
            second = total - first
            if first < second and gcd(first, total) == 1:
                result.append((first, second))
    return tuple(sorted(result))


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "no assert statements")
    require(len(BODY_SET) == 11 and min(BODY) == 1 and max(BODY) == 22,
            "body cardinality and primitive origin")

    body_owners = endpoint_certificate(BODY, WINDOW)
    length = WINDOW[1] - WINDOW[0]
    require(body_owners == (15, 4) and length == Fraction(3, 280),
            "window geometry")
    body_endpoint_table = tuple(
        (speed, rat(torus_distance(speed * WINDOW[0])),
         rat(torus_distance(speed * WINDOW[1])))
        for speed in BODY
    )

    expected_high_table = (
        (23, (1, 161)), (24, (0, 1)), (25, (1, 200)), (26, (0, 1)),
        (27, (5, 1512)), (28, (0, 1)), (29, (3, 1624)), (30, (0, 1)),
        (31, (1, 1736)), (32, (0, 1)), (33, (0, 1)), (34, (3, 2380)),
        (35, (0, 1)), (36, (1, 360)), (37, (0, 1)), (38, (1, 266)),
        (39, (0, 1)), (40, (1, 280)), (41, (0, 1)),
    )
    cumulative_table = tuple(
        (speed, rat(antiderivative_measure(speed, WINDOW)))
        for speed in range(23, 42)
    )
    require(cumulative_table == expected_high_table,
            "independent cumulative high-tail table")
    majorant_at_27 = length / 7 + Fraction(6, 49 * 27)
    require(majorant_at_27 == Fraction(107, 17640) < Fraction(1, 161),
            "independent one-tail ceiling")
    cutoff_upper = Fraction(1, 161) + length / 7 + Fraction(6, 49 * 42)
    require(cutoff_upper == Fraction(3363, 315560), "cutoff upper value")
    require(length - cutoff_upper == Fraction(9, 157780) > 0,
            "cutoff strictness")

    grid_pairs = tuple(combinations(range(23, 42), 2))
    require(len(grid_pairs) == 171, "small high-tail pair count")
    grid_witnesses = {}
    for first, second in grid_pairs:
        witness = common_grid_witness(first, second)
        require(witness is not None, ("common-grid pair witness", first, second))
        require(row_clearance(BODY + (first, second), witness) >= DELTA,
                ("common-grid full-row consequence", first, second, witness))
        grid_witnesses[(first, second)] = witness
    sum_exceptions = tuple(
        (first, second)
        for first, second in grid_pairs
        if antiderivative_measure(first, WINDOW)
        + antiderivative_measure(second, WINDOW) >= length
    )
    require(sum_exceptions == ((23, 25),), "unique measure-sum exception")
    exceptional_phase = Fraction(381, 805)
    exceptional_clearance = row_clearance(BODY + (23, 25), exceptional_phase)
    require(exceptional_clearance == Fraction(16, 161),
            "independent strict exceptional witness")

    low_intervals = {
        2: (Fraction(29, 196), Fraction(13, 84)),
        3: WINDOW,
        5: WINDOW,
        7: WINDOW,
        9: WINDOW,
        11: WINDOW,
        13: WINDOW,
        17: (Fraction(113, 238), Fraction(27, 56)),
        19: (Fraction(29, 84), Fraction(69, 196)),
        20: WINDOW,
        21: (Fraction(29, 196), Fraction(13, 84)),
    }
    low_interval_table = []
    for added in LOW_MISSING:
        interval = low_intervals[added]
        owners = endpoint_certificate(BODY + (added,), interval)
        interval_length = interval[1] - interval[0]
        low_interval_table.append((
            added,
            rat(interval[0]),
            rat(interval[1]),
            rat(interval_length),
            owners,
            ceil(Fraction(1, 7) / interval_length),
        ))
    expected_low_interval_table = (
        (2, (29, 196), (13, 84), (1, 147), (14, 6), 21),
        (3, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (5, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (7, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (9, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (11, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (13, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (17, (113, 238), (27, 56), (1, 136), (17, 4), 20),
        (19, (29, 84), (69, 196), (1, 147), (6, 14), 21),
        (20, (33, 70), (27, 56), (3, 280), (15, 4), 14),
        (21, (29, 196), (13, 84), (1, 147), (14, 6), 21),
    )
    require(tuple(low_interval_table) == expected_low_interval_table,
            "independent low interval table")

    residual_pairs = tuple(
        (first, second)
        for first, second in combinations(LOW_MISSING, 2)
        if second * (low_intervals[first][1] - low_intervals[first][0])
        < Fraction(1, 7)
    )
    at_thirteen = tuple(sorted(
        tuple(combinations((2, 3, 5, 7, 9, 11), 2))
        + ((2, 17), (2, 19), (2, 20), (17, 19), (19, 20))
    ))
    at_twenty = ((2, 13),)
    at_nineteen = ((3, 13), (5, 13), (7, 13), (9, 13), (11, 13))
    require(tuple(sorted(at_thirteen + at_twenty + at_nineteen)) == residual_pairs,
            "independent residual partition")
    low_witness_groups = (
        ((1, 13), (1, 13), at_thirteen),
        ((7, 20), (1, 10), at_twenty),
        ((9, 19), (2, 19), at_nineteen),
    )
    for phase_data, clearance_data, pairs in low_witness_groups:
        phase = Fraction(*phase_data)
        for first, second in pairs:
            require(row_clearance(BODY + (first, second), phase)
                    == Fraction(*clearance_data),
                    ("independent named low witness", first, second, phase))

    external = tuple(value for value in range(1, 501) if value not in BODY_SET)
    routed = 0
    for first, second in combinations(external, 2):
        if first > 22:
            if second < 42:
                require((first, second) in grid_witnesses,
                        ("small high pair routed by grid", first, second))
            else:
                first_ceiling = (
                    antiderivative_measure(first, WINDOW)
                    if first <= 26 else majorant_at_27
                )
                second_ceiling = length / 7 + Fraction(6, 49 * second)
                require(first_ceiling + second_ceiling < length,
                        ("analytic high pair route", first, second))
        else:
            interval_length = low_intervals[first][1] - low_intervals[first][0]
            if second * interval_length < Fraction(1, 7):
                require((first, second) in residual_pairs,
                        ("low finite route", first, second))
            else:
                require(second * interval_length >= Fraction(1, 7),
                        ("compact-to-open route", first, second))
        routed += 1
    require(routed == 119316, "independent hostile rectangle routing")

    atlas = atlas_sieve()
    require(len(atlas) == 5855, "independent atlas count")
    require(all(
        first < second
        and 23 * first > 22
        and 23 * second > 22
        and (
            23 * second >= 42
            or (23 * first, 23 * second) in grid_witnesses
        )
        for first, second in atlas
    ), "every atlas type enters the universal high-tail proof")
    require(all(pair in atlas for pair in LIVE_SCALE_ONE_TYPES),
            "independent live-type membership")
    atlas_hash = semantic_digest(atlas)

    dilation_controls = []
    for scale in (1, 2, 3, 7, 19, 64):
        phase = exceptional_phase / scale
        actual = row_clearance(
            tuple(scale * speed for speed in BODY + (23, 25)), phase
        )
        require(actual == exceptional_clearance,
                ("independent dilation control", scale))
        dilation_controls.append((scale, rat(phase), rat(actual)))

    ledger = {
        "body": BODY,
        "threshold": (1, 14),
        "universal_completion": "distinct positive a,b outside body",
        "body_interval": (rat(WINDOW[0]), rat(WINDOW[1]), rat(length)),
        "body_endpoint_owners": body_owners,
        "body_endpoint_table": body_endpoint_table,
        "danger_arc_duty": (1, 7),
        "discrepancy_excess": (6, 49),
        "high_tail_max_measure": (1, 161),
        "high_second_cutoff": 42,
        "high_second_cutoff_upper": rat(cutoff_upper),
        "high_second_cutoff_slack": rat(length - cutoff_upper),
        "high_small_measure_table": expected_high_table,
        "high_small_pair_count": 171,
        "high_sum_exception": (23, 25),
        "high_exception_witness": (rat(exceptional_phase), rat(exceptional_clearance)),
        "missing_low": LOW_MISSING,
        "low_interval_table": expected_low_interval_table,
        "low_residual_pairs": residual_pairs,
        "low_witness_groups": low_witness_groups,
        "hostile_rectangle": (500, routed, "all_closed"),
        "atlas_type_count": len(atlas),
        "atlas_sha256": atlas_hash,
        "live_scale_one_types": LIVE_SCALE_ONE_TYPES,
        "scale_two_not_claimed": (2, 1, 9),
        "dilation_controls": tuple(dilation_controls),
    }
    semantic = semantic_digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "shared frozen semantic digest")

    print("status=PASS")
    print("implementation=periodic_antiderivative_plus_exact_common_wall_grid")
    print(f"body_interval={rat(WINDOW[0])},{rat(WINDOW[1])};length={rat(length)};owners={body_owners}")
    print(f"cumulative_high_measure_table={expected_high_table}")
    print(f"common_grid_high_pairs={len(grid_pairs)};measure_sum_exception={(23,25)}")
    print(f"exceptional_strict_witness={rat(exceptional_phase)};clearance={rat(exceptional_clearance)}")
    print(f"low_interval_table={expected_low_interval_table}")
    print(f"low_residual_pairs={residual_pairs}")
    print(f"low_witness_groups={low_witness_groups}")
    print(f"hostile_rectangle_routes=external_speeds_le_500;pairs={routed}")
    print(f"thm3878_scale_one_atlas_types={len(atlas)};atlas_sha256={atlas_hash}")
    print(f"live_scale_one_types={LIVE_SCALE_ONE_TYPES};scale_two_not_claimed={(2,1,9)}")
    print(f"dilation_controls={tuple(dilation_controls)}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
