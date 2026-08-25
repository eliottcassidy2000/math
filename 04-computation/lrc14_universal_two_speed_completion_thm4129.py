#!/usr/bin/env python3
"""Exact interval-comb referee for THM-4129."""

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
THRESHOLD = Fraction(1, 14)
J = (Fraction(33, 70), Fraction(27, 56))
MISSING_LOW = (2, 3, 5, 7, 9, 11, 13, 17, 19, 20, 21)
LIVE_SCALE_ONE_TYPES = (
    (1, 3), (1, 4), (1, 9), (1, 10), (2, 3), (2, 9),
    (3, 7), (3, 8), (4, 7), (5, 6), (5, 12), (6, 11),
    (7, 10), (8, 9), (8, 21), (9, 11),
)
EXPECTED_SEMANTIC = "a717248ff237c9cbadf1a15dbaeb2902c0b7b7da10ebca0e29d797ea09bbdbe5"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def rational(value: Fraction) -> tuple[int, int]:
    return (value.numerator, value.denominator)


def digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(values: tuple[int, ...], phase: Fraction) -> Fraction:
    return min(distance(value * phase) for value in values)


def certify_safe_interval(
    values: tuple[int, ...], interval: tuple[Fraction, Fraction]
) -> tuple[int, int]:
    """Certify a closed safe interval by one-cell concavity."""
    left, right = interval
    require(left < right, ("positive interval", interval))
    for value in values:
        require(floor(value * left) == floor(value * right),
                ("one integer cell", value, interval))
        require(distance(value * left) >= THRESHOLD,
                ("left endpoint safe", value, interval))
        require(distance(value * right) >= THRESHOLD,
                ("right endpoint safe", value, interval))
    left_min = min((distance(value * left), value) for value in values)
    right_min = min((distance(value * right), value) for value in values)
    require(left_min[0] == THRESHOLD and right_min[0] == THRESHOLD,
            ("sharp endpoint clearances", values, interval))
    left_owners = tuple(value for value in values
                        if distance(value * left) == THRESHOLD)
    right_owners = tuple(value for value in values
                         if distance(value * right) == THRESHOLD)
    require(len(left_owners) == len(right_owners) == 1,
            ("unique endpoint owners", values, interval))
    return (left_owners[0], right_owners[0])


def danger_intervals(
    speed: int, interval: tuple[Fraction, Fraction]
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Strict danger components clipped to a real interval; endpoints have zero mass."""
    left, right = interval
    pieces = []
    for integer in range(floor(speed * left) - 1, ceil(speed * right) + 2):
        lower = max(left, Fraction(14 * integer - 1, 14 * speed))
        upper = min(right, Fraction(14 * integer + 1, 14 * speed))
        if lower < upper:
            pieces.append((lower, upper))
    return tuple(pieces)


def danger_measure(speed: int, interval: tuple[Fraction, Fraction]) -> Fraction:
    return sum((right - left for left, right in danger_intervals(speed, interval)),
               Fraction(0))


def discrepancy_majorant(speed: int, length: Fraction) -> Fraction:
    return length / 7 + Fraction(6, 49 * speed)


def interval_witness(
    interval: tuple[Fraction, Fraction], first: int, second: int
) -> Fraction | None:
    """Find a positive complementary interval to two strict danger combs."""
    left, right = interval
    current = left
    pieces = sorted(danger_intervals(first, interval) + danger_intervals(second, interval))
    for lower, upper in pieces:
        if lower > current:
            return (current + lower) / 2
        if upper > current:
            current = upper
    if right > current:
        return (current + right) / 2
    return None


def factor_trial(value: int) -> tuple[tuple[int, int], ...]:
    factors = []
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            exponent = 0
            while value % prime == 0:
                value //= prime
                exponent += 1
            factors.append((prime, exponent))
        prime += 1
    if value > 1:
        factors.append((value, 1))
    return tuple(factors)


def admissible_atlas_sum(value: int) -> bool:
    factors = factor_trial(value)
    return bool(factors) and all(prime % 3 == 2 and exponent <= 2
                                 for prime, exponent in factors)


def atlas_trial() -> tuple[tuple[int, int], ...]:
    return tuple(
        (first, second)
        for first in range(1, 356)
        for second in range(first + 1, 357 - first)
        if first + second <= 356
        and gcd(first, second) == 1
        and admissible_atlas_sum(first + second)
    )


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "no assert statements")
    require(len(BODY) == 11 and len(BODY_SET) == 11 and 1 in BODY_SET,
            "primitive eleven-speed body")

    body_owners = certify_safe_interval(BODY, J)
    require(body_owners == (15, 4), "sharp J endpoint owners")
    require(J[1] - J[0] == Fraction(3, 280), "J length")
    body_endpoint_table = tuple(
        (value, rational(distance(value * J[0])), rational(distance(value * J[1])))
        for value in BODY
    )

    expected_high_table = (
        (23, (1, 161)), (24, (0, 1)), (25, (1, 200)), (26, (0, 1)),
        (27, (5, 1512)), (28, (0, 1)), (29, (3, 1624)), (30, (0, 1)),
        (31, (1, 1736)), (32, (0, 1)), (33, (0, 1)), (34, (3, 2380)),
        (35, (0, 1)), (36, (1, 360)), (37, (0, 1)), (38, (1, 266)),
        (39, (0, 1)), (40, (1, 280)), (41, (0, 1)),
    )
    actual_high_table = tuple(
        (speed, rational(danger_measure(speed, J))) for speed in range(23, 42)
    )
    require(actual_high_table == expected_high_table, "exact low end of high-tail table")

    length = J[1] - J[0]
    require(discrepancy_majorant(27, length) == Fraction(107, 17640),
            "v=27 discrepancy value")
    require(discrepancy_majorant(27, length) < Fraction(1, 161),
            "universal one-tail measure ceiling after v=26")
    require(tuple(rational(danger_measure(speed, J)) for speed in range(23, 27))
            == ((1, 161), (0, 1), (1, 200), (0, 1)),
            "four initial one-tail measures")
    high_pair_exceptions = tuple(
        (first, second)
        for first, second in combinations(range(23, 42), 2)
        if danger_measure(first, J) + danger_measure(second, J) >= length
    )
    require(high_pair_exceptions == ((23, 25),), "unique small-tail sum exception")
    cutoff_upper = Fraction(1, 161) + discrepancy_majorant(42, length)
    require(cutoff_upper == Fraction(3363, 315560), "high second-tail cutoff value")
    require(length - cutoff_upper == Fraction(9, 157780) > 0,
            "strict high second-tail cutoff slack")
    exceptional_phase = Fraction(381, 805)
    exceptional_clearance = clearance(BODY + (23, 25), exceptional_phase)
    require(J[0] < exceptional_phase < J[1], "exceptional witness inside J")
    require(exceptional_clearance == Fraction(16, 161) > THRESHOLD,
            "strict exceptional-pair witness")

    low_intervals = {
        2: (Fraction(29, 196), Fraction(13, 84)),
        3: J,
        5: J,
        7: J,
        9: J,
        11: J,
        13: J,
        17: (Fraction(113, 238), Fraction(27, 56)),
        19: (Fraction(29, 84), Fraction(69, 196)),
        20: J,
        21: (Fraction(29, 196), Fraction(13, 84)),
    }
    require(tuple(sorted(low_intervals)) == MISSING_LOW, "complete low complement")
    low_interval_table = []
    for added in MISSING_LOW:
        interval = low_intervals[added]
        owners = certify_safe_interval(BODY + (added,), interval)
        interval_length = interval[1] - interval[0]
        cutoff = ceil(Fraction(1, 7) / interval_length)
        low_interval_table.append((
            added,
            rational(interval[0]),
            rational(interval[1]),
            rational(interval_length),
            owners,
            cutoff,
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
            "exact low safe-interval table")

    residual_pairs = tuple(
        (first, second)
        for first, second in combinations(MISSING_LOW, 2)
        if second * (low_intervals[first][1] - low_intervals[first][0])
        < Fraction(1, 7)
    )
    at_thirteen = tuple(sorted(
        tuple(combinations((2, 3, 5, 7, 9, 11), 2))
        + ((2, 17), (2, 19), (2, 20), (17, 19), (19, 20))
    ))
    at_twenty = ((2, 13),)
    at_nineteen = ((3, 13), (5, 13), (7, 13), (9, 13), (11, 13))
    require(len(residual_pairs) == 26, "low residual count")
    require(tuple(sorted(at_thirteen + at_twenty + at_nineteen)) == residual_pairs,
            "three-clock low residual partition")
    low_witness_groups = (
        ((1, 13), (1, 13), at_thirteen),
        ((7, 20), (1, 10), at_twenty),
        ((9, 19), (2, 19), at_nineteen),
    )
    for phase_pair, expected_clearance, pairs in low_witness_groups:
        phase = Fraction(*phase_pair)
        for first, second in pairs:
            require(clearance(BODY + (first, second), phase)
                    == Fraction(*expected_clearance),
                    ("named low-pair witness", first, second, phase))

    # Direct hostile rectangle: print the theorem's consequence, not only its bounds.
    external = tuple(value for value in range(1, 501) if value not in BODY_SET)
    hostile_checked = 0
    minimum_record = None
    for first, second in combinations(external, 2):
        if first > 22:
            phase = interval_witness(J, first, second)
        else:
            interval = low_intervals[first]
            if second * (interval[1] - interval[0]) >= Fraction(1, 7):
                pieces = danger_intervals(second, interval)
                current = interval[0]
                phase = None
                for lower, upper in pieces:
                    if lower > current:
                        phase = (current + lower) / 2
                        break
                    if upper > current:
                        current = upper
                if phase is None and interval[1] > current:
                    phase = (current + interval[1]) / 2
            elif (first, second) in at_thirteen:
                phase = Fraction(1, 13)
            elif (first, second) in at_twenty:
                phase = Fraction(7, 20)
            elif (first, second) in at_nineteen:
                phase = Fraction(9, 19)
            else:
                raise RuntimeError(f"FAILED: unrouted hostile pair {(first, second)}")
        require(phase is not None, ("hostile rectangle witness exists", first, second))
        actual = clearance(BODY + (first, second), phase)
        require(actual >= THRESHOLD,
                ("hostile rectangle consequence", first, second, phase, actual))
        record = (actual, first, second, phase)
        if minimum_record is None or record < minimum_record:
            minimum_record = record
        hostile_checked += 1
    require(hostile_checked == 119316, "hostile rectangle pair count")
    require(minimum_record == (
        Fraction(995, 13916), 34, 497, Fraction(223651, 473144)
    ), "frozen hostile rectangle minimum witness")

    atlas = atlas_trial()
    require(len(atlas) == 5855, "all THM-3878 ratio types")
    require(all(
        first < second
        and 23 * first not in BODY_SET
        and 23 * second not in BODY_SET
        and interval_witness(J, 23 * first, 23 * second) is not None
        for first, second in atlas
    ), "canonical representative of every scale-one atlas type")
    require(all(pair in atlas for pair in LIVE_SCALE_ONE_TYPES),
            "sixteen live scale-one types occur in atlas")
    atlas_digest = digest(atlas)

    dilation_controls = []
    for scale in (1, 2, 3, 7, 19, 64):
        phase = exceptional_phase / scale
        actual = clearance(tuple(scale * value for value in BODY + (23, 25)), phase)
        require(actual == exceptional_clearance, ("dilation control", scale))
        dilation_controls.append((scale, rational(phase), rational(actual)))

    ledger = {
        "body": BODY,
        "threshold": (1, 14),
        "universal_completion": "distinct positive a,b outside body",
        "body_interval": (rational(J[0]), rational(J[1]), rational(length)),
        "body_endpoint_owners": body_owners,
        "body_endpoint_table": body_endpoint_table,
        "danger_arc_duty": (1, 7),
        "discrepancy_excess": (6, 49),
        "high_tail_max_measure": (1, 161),
        "high_second_cutoff": 42,
        "high_second_cutoff_upper": rational(cutoff_upper),
        "high_second_cutoff_slack": rational(length - cutoff_upper),
        "high_small_measure_table": expected_high_table,
        "high_small_pair_count": 171,
        "high_sum_exception": (23, 25),
        "high_exception_witness": (rational(exceptional_phase),
                                    rational(exceptional_clearance)),
        "missing_low": MISSING_LOW,
        "low_interval_table": expected_low_interval_table,
        "low_residual_pairs": residual_pairs,
        "low_witness_groups": low_witness_groups,
        "hostile_rectangle": (500, hostile_checked, "all_closed"),
        "atlas_type_count": len(atlas),
        "atlas_sha256": atlas_digest,
        "live_scale_one_types": LIVE_SCALE_ONE_TYPES,
        "scale_two_not_claimed": (2, 1, 9),
        "dilation_controls": tuple(dilation_controls),
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=PASS")
    print("theorem=every_distinct_two_speed_completion_of_U_is_1/14_lonely")
    print(f"body_interval={rational(J[0])},{rational(J[1])};length={rational(length)};owners={body_owners}")
    print(f"high_small_measure_table={expected_high_table}")
    print(f"high_pair_exception={(23,25)};witness={rational(exceptional_phase)};clearance={rational(exceptional_clearance)}")
    print(f"high_second_cutoff=42;strict_slack={rational(length-cutoff_upper)}")
    print(f"low_interval_table={expected_low_interval_table}")
    print(f"low_residual_pairs={residual_pairs}")
    print(f"low_witness_groups={low_witness_groups}")
    print("hostile_rectangle=external_speeds_le_500;"
          f"pairs={hostile_checked};minimum={(rational(minimum_record[0]), minimum_record[1], minimum_record[2], rational(minimum_record[3]))}")
    print(f"thm3878_scale_one_atlas_types={len(atlas)};atlas_sha256={atlas_digest}")
    print(f"live_scale_one_types={LIVE_SCALE_ONE_TYPES};scale_two_not_claimed={(2,1,9)}")
    print(f"dilation_controls={tuple(dilation_controls)}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
