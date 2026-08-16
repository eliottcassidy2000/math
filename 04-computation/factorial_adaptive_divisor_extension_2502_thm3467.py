#!/usr/bin/env python3
"""Dual exact adaptive-divisor audit of the factorial boundary row d=2502."""

import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path


D = 2502
PRIMARY_NAME = "factorial_seven_exit_first_flag_2400_thm3201.py"
PRIMARY_SHA256 = "1cf067eb3cc8da29e51ee8f560f4100fb66a4b9684a2da8b1f7ef14d380e4171"
INDEPENDENT_NAME = "factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py"
INDEPENDENT_SHA256 = "afe3ec56c7efe99ff4a09c1af856b4e144312d48ec45e729d063a5fd83b0f0b0"
EXPECTED_SEMANTIC_DIGEST = "c767afd684dbf63910db607a742c530ed9f03a8d61c97b6e2e4dd5cb22bddf98"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_pinned(filename, expected_hash, module_name):
    source = Path(__file__).resolve().with_name(filename)
    require(source.is_file(), ("missing dependency", source))
    normalized = source.read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(normalized).hexdigest()
    require(actual_hash == expected_hash, (source, actual_hash, expected_hash))
    tree = ast.parse(normalized.decode("utf-8"), filename=str(source))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            (source, "assert node"))
    spec = importlib.util.spec_from_file_location(module_name, source)
    require(spec is not None and spec.loader is not None, ("bad import spec", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalized_primary_profile(companion, row, prime):
    slopes = companion.slope_lengths(row, prime)
    return (
        companion.zero_order(row),
        tuple(
            (slope.numerator, slope.denominator, slopes[slope])
            for slope in sorted(slopes)
        ),
    )


def normalized_independent_profile(engine, row, prime):
    zero_order, slopes = engine.determinant_lower_profile(row, prime)
    return (
        zero_order,
        tuple(
            (slope[0], slope[1], slopes[slope])
            for slope in sorted(slopes, key=lambda item: Fraction(item[0], item[1]))
        ),
    )


def common_blocks(profiles):
    common = {entry[:2] for entry in profiles[0][1]}
    for _, entries in profiles[1:]:
        common &= {entry[:2] for entry in entries}
    capacities = tuple(
        {entry[:2]: entry[2] for entry in entries} for _, entries in profiles
    )
    return tuple(
        (slope[0], slope[1], min(capacity[slope] for capacity in capacities), slope[1])
        for slope in sorted(common, key=lambda item: Fraction(item[0], item[1]))
    )


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "thm3201_primary_for_d2502")
    independent = load_pinned(
        INDEPENDENT_NAME,
        INDEPENDENT_SHA256,
        "thm3201_independent_for_d2502",
    )
    primary_base = primary.load_engine()
    companion = primary_base.load_companion()
    independent_engine = independent.load_engine()

    planted, _, planted_zero, _ = companion.planted_control()
    require(1 in planted and 1 in planted_zero, ("primary planted controls", planted, planted_zero))
    independent_control = independent_engine.positive_controls()
    require(
        independent_control == "planted(v+1) and planted(v), degree 1 retained",
        independent_control,
    )

    require(primary.seven_exit_residual(primary_base, companion, D), "primary residual filter")
    require(independent.seven_exit_residual(independent_engine, D), "independent residual filter")
    invoice = (
        (D, ((2, 1), (3, 2), (139, 1))),
        (D - 1, ((41, 1), (61, 1))),
        (D - 2, ((2, 2), (5, 4))),
        (D - 3, ((3, 1), (7, 2), (17, 1))),
        (D - 4, ((2, 1), (1249, 1))),
        (D - 5, ((11, 1), (227, 1))),
        (D - 6, ((2, 6), (3, 1), (13, 1))),
    )
    require(
        tuple((n, tuple(companion.factor(n))) for n, _ in invoice) == invoice,
        "factorization invoice",
    )

    p_primary, q_primary = companion.pair(D)
    e_primary = companion.first_full_pseudoremainder(p_primary, q_primary, D)
    primary_rows = (p_primary, q_primary, e_primary)
    p_independent, q_independent = independent_engine.moments(D)
    e_independent = independent_engine.first_full_remainder(
        p_independent, q_independent, D
    )
    independent_rows = (p_independent, q_independent, e_independent)
    expected_degrees = (D - 2, D - 1, D - 3)
    require(tuple(row.degree() for row in primary_rows) == expected_degrees, "primary degrees")
    require(tuple(row.degree() for row in independent_rows) == expected_degrees,
            "independent degrees")
    require(primary_rows == independent_rows, "dual polynomial reconstruction mismatch")

    primes = (2, 3, 5, 41, 61)
    profiles = {}
    blocks = {}
    local_degrees = {}
    for prime in primes:
        primary_profiles = tuple(
            normalized_primary_profile(companion, row, prime) for row in primary_rows
        )
        independent_profiles = tuple(
            normalized_independent_profile(independent_engine, row, prime)
            for row in independent_rows
        )
        require(primary_profiles == independent_profiles,
                (prime, "raw profile mismatch", primary_profiles, independent_profiles))
        profiles[prime] = primary_profiles
        blocks[prime] = common_blocks(primary_profiles)
        primary_local, _ = companion.degree_barcode_many(primary_rows, prime)
        independent_local = independent_engine.allowed_degrees(independent_rows, prime)
        require(primary_local == independent_local,
                (prime, "local degree mismatch", primary_local, independent_local))
        local_degrees[prime] = tuple(sorted(primary_local))

    expected_profiles_41 = (
        (0, ((0, 1, 20), (2, 41, 779), (84, 1681, 1681), (1, 20, 20))),
        (0, ((2, 41, 820), (84, 1681, 1681))),
        (0, ((1, 21, 21), (2, 41, 779), (84, 1681, 1681), (1, 18, 18))),
    )
    expected_profiles_61 = (
        (0, ((0, 1, 30), (2, 61, 1830), (11, 320, 640))),
        (0, ((2, 61, 1830), (1, 30, 30), (22, 641, 641))),
        (0, ((1, 31, 31), (2, 61, 1830), (21, 610, 610), (1, 28, 28))),
    )
    require(profiles[41] == expected_profiles_41, (41, profiles[41]))
    require(profiles[61] == expected_profiles_61, (61, profiles[61]))
    require(blocks[41] == ((2, 41, 779, 41), (84, 1681, 1681, 1681)), blocks[41])
    require(blocks[61] == ((2, 61, 1830, 61),), blocks[61])

    positive_universe = set(range(1, D - 2))
    hostile_trace = []
    possible = positive_universe.copy()
    for prime in (61, 2, 3):
        possible &= set(local_degrees[prime])
        hostile_trace.append((prime, tuple(sorted(possible))))
    require(hostile_trace[-1] == (3, (61,)), hostile_trace)

    divisor_trace = []
    possible = positive_universe.copy()
    for prime in (61, 41):
        possible &= set(local_degrees[prime])
        divisor_trace.append((prime, tuple(sorted(possible))))
    require(not possible, divisor_trace)

    fixed_counts = []
    possible = positive_universe.copy()
    for prime in (2, 3, 5):
        possible &= set(local_degrees[prime])
        fixed_counts.append((prime, len(possible)))
    require(fixed_counts == [(2, 2499), (3, 223), (5, 0)], fixed_counts)

    semantic_record = {
        "d": D,
        "invoice": invoice,
        "degrees": expected_degrees,
        "profiles": profiles,
        "blocks": blocks,
        "local_degrees": local_degrees,
        "hostile_trace": hostile_trace,
        "divisor_trace": divisor_trace,
        "fixed_counts": fixed_counts,
    }
    semantic_digest = digest(semantic_record)
    require(
        EXPECTED_SEMANTIC_DIGEST == "TO_BE_FILLED"
        or semantic_digest == EXPECTED_SEMANTIC_DIGEST,
        (semantic_digest, EXPECTED_SEMANTIC_DIGEST),
    )

    print("THM-3467 ADAPTIVE DIVISOR EXTENSION d=2502 DUAL EXACT AUDIT")
    print("primary=%s sha256=%s" % (PRIMARY_NAME, PRIMARY_SHA256))
    print("independent=%s sha256=%s" % (INDEPENDENT_NAME, INDEPENDENT_SHA256))
    print("universe=single seven-exit residual d=2502; rows=(F,G,E); places=(2,3,5,41,61)")
    print("factorization_invoice=%s" % (invoice,))
    print("degrees=%s dual_polynomials_equal=True" % (expected_degrees,))
    print("p41_raw_profiles=%s" % (profiles[41],))
    print("p61_raw_profiles=%s" % (profiles[61],))
    print("p41_common_blocks=%s" % (blocks[41],))
    print("p61_common_blocks=%s" % (blocks[61],))
    print("p41_positive_degrees=%s" % (tuple(x for x in local_degrees[41] if x),))
    print("p61_positive_degrees=%s" % (tuple(x for x in local_degrees[61] if x),))
    print("hostile_reset_2_3_trace_sizes=(30,30,1) survivor=(61,)")
    print("adaptive_divisor_trace_sizes=(30,0) survivor=()")
    print("fixed_2_3_5_counts=%s" % (tuple(fixed_counts),))
    print("controls=dual raw profiles/local degrees; planted(v+1),planted(v) retain degree 1")
    print("semantic_sha256=%s" % semantic_digest)
    print("consequence=with THM-3201 and d=2501, every exact quadratic {0,1,2} window 1<=r<=2500 has a nonzero moment")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
