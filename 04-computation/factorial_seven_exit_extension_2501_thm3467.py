#!/usr/bin/env python3
"""Dual exact audit of the former THM-3201 boundary row d=2501."""

import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import prod
from pathlib import Path


D = 2501
PRIMARY_NAME = "factorial_seven_exit_first_flag_2400_thm3201.py"
PRIMARY_SHA256 = "1cf067eb3cc8da29e51ee8f560f4100fb66a4b9684a2da8b1f7ef14d380e4171"
INDEPENDENT_NAME = "factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py"
INDEPENDENT_SHA256 = "afe3ec56c7efe99ff4a09c1af856b4e144312d48ec45e729d063a5fd83b0f0b0"
EXPECTED_TRACE = (
    (2, (256, 2048, 2304)),
    (3, (256,)),
    (5, ()),
)
EXPECTED_TRACE_DIGEST = "365533925519a4d8d44db78394f0785e87be5f4cc03e0a98d759f93609fb09ee"
STRUCTURAL_PARAMETERS = (
    (5, 1),
    (5, 2),
    (5, 3),
    (7, 1),
    (7, 2),
    (11, 1),
    (11, 2),
    (13, 1),
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_pinned(filename, expected_hash, module_name):
    source = Path(__file__).resolve().with_name(filename)
    require(source.is_file(), ("missing dependency", source))
    actual_hash = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual_hash == expected_hash, (source, actual_hash, expected_hash))
    spec = importlib.util.spec_from_file_location(module_name, source)
    require(spec is not None and spec.loader is not None, ("bad import spec", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def primary_blocks(companion, rows, prime):
    _, blocks = companion.degree_barcode_many(rows, prime)
    return tuple(
        ((slope.numerator, slope.denominator), capacity, denominator)
        for slope, capacity, denominator in blocks
    )


def independent_blocks(engine, rows, prime):
    profiles = tuple(engine.determinant_lower_profile(row, prime) for row in rows)
    common = set(profiles[0][1])
    for _, capacities in profiles[1:]:
        common &= set(capacities)
    return tuple(
        (slope, min(capacities[slope] for _, capacities in profiles), slope[1])
        for slope in sorted(common, key=lambda item: Fraction(item[0], item[1]))
    )


def faces_with_slope(companion, polynomial, prime, slope):
    hull = companion.lower_hull(polynomial, prime)
    return tuple(
        (left, right)
        for left, right in zip(hull, hull[1:])
        if Fraction(right[1] - left[1], right[0] - left[0]) == slope
    )


def structural_carry_audit(companion):
    records = []
    for prime, depth in STRUCTURAL_PARAMETERS:
        height = prime**depth
        for multiplier in range(2, prime):
            size = multiplier * height
            d = size + 1
            midpoint = (height - 1) // 2
            carry = 2 * (height - 1) // (prime - 1)
            intercept = (height - 1) // (prime - 1) - depth
            extent = min(multiplier - 1, (prime - 1) // 2)
            g_extent = min(multiplier, (prime - 1) // 2)
            slope = Fraction(carry, height)
            f_row, g_row = companion.pair(d)
            e_row = companion.first_full_pseudoremainder(f_row, g_row, d)
            require(
                (f_row.degree(), g_row.degree(), e_row.degree())
                == (size - 1, size, size - 2),
                (prime, depth, multiplier, "degrees"),
            )
            f_faces = faces_with_slope(companion, f_row, prime, slope)
            g_faces = faces_with_slope(companion, g_row, prime, slope)
            e_faces = faces_with_slope(companion, e_row, prime, slope)
            expected_faces = (
                (((midpoint, intercept),
                  (midpoint + extent * height, intercept + extent * carry)),),
                (((0, 0), (g_extent * height, g_extent * carry)),),
                (((midpoint + 1, intercept + depth),
                  (midpoint + 1 + extent * height,
                   intercept + depth + extent * carry)),),
            )
            require(
                (f_faces, g_faces, e_faces) == expected_faces,
                (prime, depth, multiplier, "faces"),
            )
            f_anchors = tuple(
                (midpoint + t * height,
                 companion.vp(f_row[midpoint + t * height], prime))
                for t in range(extent + 1)
            )
            g_anchors = tuple(
                (t * height, companion.vp(g_row[t * height], prime))
                for t in range(g_extent + 1)
            )
            e_anchors = tuple(
                (midpoint + 1 + t * height,
                 companion.vp(e_row[midpoint + 1 + t * height], prime))
                for t in range(extent + 1)
            )
            require(
                f_anchors
                == tuple(
                    (midpoint + t * height, intercept + t * carry)
                    for t in range(extent + 1)
                ),
                (prime, depth, multiplier, "F anchors"),
            )
            require(
                g_anchors
                == tuple(
                    (t * height, t * carry) for t in range(g_extent + 1)
                ),
                (prime, depth, multiplier, "G anchors"),
            )
            require(
                e_anchors
                == tuple(
                    (midpoint + 1 + t * height,
                     intercept + depth + t * carry)
                    for t in range(extent + 1)
                ),
                (prime, depth, multiplier, "E anchors"),
            )
            local, blocks = companion.degree_barcode_many(
                (f_row, g_row, e_row), prime
            )
            require(
                blocks == ((slope, extent * height, height),),
                (prime, depth, multiplier, "blocks", blocks),
            )
            require(
                local == set(range(0, extent * height + 1, height)),
                (prime, depth, multiplier, "local", tuple(sorted(local))),
            )
            require(carry % prime == 2 and slope.denominator == height, (
                prime,
                depth,
                multiplier,
                "primitive slope",
            ))
            records.append((
                prime,
                depth,
                multiplier,
                d,
                str(slope),
                extent * height,
                height,
                f_faces,
                g_faces,
                e_faces,
                f_anchors,
                g_anchors,
                e_anchors,
                tuple((str(s), capacity, denominator)
                      for s, capacity, denominator in blocks),
                tuple(sorted(local)),
            ))
    require(len(records) == 48, ("structural cell count", len(records)))

    hostile_parameters = ((3, 3, 2), (5, 1, 6), (5, 2, 1))
    hostiles = []
    for prime, depth, multiplier in hostile_parameters:
        height = prime**depth
        d = multiplier * height + 1
        f_row, g_row = companion.pair(d)
        e_row = companion.first_full_pseudoremainder(f_row, g_row, d)
        local, blocks = companion.degree_barcode_many(
            (f_row, g_row, e_row), prime
        )
        entries = tuple(
            (
                str(slope),
                capacity,
                denominator,
                faces_with_slope(companion, f_row, prime, slope),
                faces_with_slope(companion, g_row, prime, slope),
                faces_with_slope(companion, e_row, prime, slope),
            )
            for slope, capacity, denominator in blocks
        )
        hostiles.append((
            prime,
            depth,
            multiplier,
            d,
            entries,
            tuple(sorted(local)),
        ))
    hostiles = tuple(hostiles)
    expected_hostiles = (
        (
            3,
            3,
            2,
            55,
            ((
                "26/27",
                27,
                27,
                (((13, 10), (40, 36)),),
                (((0, 0), (27, 26)),),
                (((13, 13), (40, 39)),),
            ),),
            (0, 27),
        ),
        (
            5,
            1,
            6,
            31,
            ((
                "12/25",
                25,
                25,
                (((2, 0), (27, 12)),),
                (((5, 2), (30, 14)),),
                (((3, 1), (28, 13)),),
            ),),
            (0, 25),
        ),
        (5, 2, 1, 26, (), (0,)),
    )
    require(hostiles == expected_hostiles, ("structural hostiles", hostiles))
    return tuple(records), hostiles


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "thm3201_primary_for_thm3467")
    independent = load_pinned(
        INDEPENDENT_NAME,
        INDEPENDENT_SHA256,
        "thm3201_independent_for_thm3467",
    )

    primary_engine = primary.load_engine()
    primary_companion = primary_engine.load_companion()
    independent_engine = independent.load_engine()

    structural_records, structural_hostiles = structural_carry_audit(
        primary_companion
    )

    primary_counts, primary_residuals = primary.exit_counts(
        primary_engine, primary_companion, D, D
    )
    independent_counts, independent_residuals = independent.exit_counts(
        independent_engine, D, D
    )
    expected_counts = (1, 1, 1, 1, 1, 1, 1)
    require(primary_counts == expected_counts, ("primary census", primary_counts))
    require(independent_counts == expected_counts, ("independent census", independent_counts))
    require(primary_residuals == independent_residuals == (D,), (
        primary_residuals,
        independent_residuals,
    ))

    p_primary, q_primary = primary_companion.pair(D)
    r_primary = primary_companion.first_full_pseudoremainder(p_primary, q_primary, D)
    primary_rows = (p_primary, q_primary, r_primary)
    require(
        tuple(row.degree() for row in primary_rows) == (D - 2, D - 1, D - 3),
        ("primary row degrees", tuple(row.degree() for row in primary_rows)),
    )
    primary_possible, primary_trace = primary_engine.flag_trace(
        primary_companion, primary_rows
    )

    p_independent, q_independent = independent_engine.moments(D)
    r_independent = independent_engine.first_full_remainder(
        p_independent, q_independent, D
    )
    independent_rows = (p_independent, q_independent, r_independent)
    require(
        tuple(row.degree() for row in independent_rows) == (D - 2, D - 1, D - 3),
        ("independent row degrees", tuple(row.degree() for row in independent_rows)),
    )
    independent_possible, independent_trace = independent_engine.flag_trace(
        independent_rows
    )

    require(not primary_possible and not independent_possible, (
        primary_possible,
        independent_possible,
    ))
    require(primary_trace == independent_trace == EXPECTED_TRACE, (
        primary_trace,
        independent_trace,
    ))
    records = ((D, primary_trace),)
    require(digest(records) == EXPECTED_TRACE_DIGEST, digest(records))

    block_ledgers = []
    degree_ledgers = []
    for prime in (2, 3, 5):
        primary_degrees, _ = primary_companion.degree_barcode_many(
            primary_rows, prime
        )
        independent_degrees = independent_engine.allowed_degrees(
            independent_rows, prime
        )
        require(primary_degrees == independent_degrees, (
            "degree ledger mismatch",
            prime,
            tuple(sorted(primary_degrees)),
            tuple(sorted(independent_degrees)),
        ))
        left = primary_blocks(primary_companion, primary_rows, prime)
        right = independent_blocks(independent_engine, independent_rows, prime)
        require(left == right, ("common-slope ledger mismatch", prime, left, right))
        block_ledgers.append((prime, left))
        degree_ledgers.append((prime, tuple(sorted(primary_degrees))))
    block_ledgers = tuple(block_ledgers)
    degree_ledgers = tuple(degree_ledgers)
    require(
        block_ledgers[-1] == (5, (((312, 625), 1250, 625),)),
        ("5-adic carry ledger", block_ledgers[-1]),
    )
    require(degree_ledgers[0] == (2, (0, 256, 2048, 2304)), degree_ledgers[0])
    require(degree_ledgers[-1] == (5, (0, 625, 1250)), degree_ledgers[-1])

    planted, _, planted_zero, _ = primary_companion.planted_control()
    require(1 in planted and 1 in planted_zero, (planted, planted_zero))
    independent_control = independent_engine.positive_controls()

    factorization_invoice = (
        (D, (41, 61)),
        (D - 1, (2, 2, 5, 5, 5, 5)),
        (D - 2, (3, 7, 7, 17)),
        (D - 3, (2, 1249)),
        (D - 4, (11, 227)),
        (D - 5, (2, 2, 2, 2, 2, 2, 3, 13)),
        (D - 6, (5, 499)),
    )
    require(
        all(number == prod(factors) for number, factors in factorization_invoice),
        factorization_invoice,
    )

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "assert node",
    )
    semantic = (
        expected_counts,
        structural_records,
        structural_hostiles,
        records,
        block_ledgers,
        degree_ledgers,
        factorization_invoice,
        independent_control,
    )

    print("THM-3467 FACTORIAL SEVEN-EXIT d=2501 DUAL EXACT AUDIT")
    print("prime_power_carry_cells=48 max_d=1211 complete_singleton_ledgers=PASS")
    print("prime_power_carry_hostiles=p3_shift;a_ge_p_depth_change;a1_empty")
    print("prime_power_carry_semantic_sha256=%s" % digest((
        structural_records,
        structural_hostiles,
    )))
    print("universe=d=2501 r=2499 seven_exit_census=%s" % (expected_counts,))
    print("row_degrees=(2499,2500,2498)")
    print("progressive_degree_trace=%s" % (primary_trace,))
    print("trace_digest=%s" % EXPECTED_TRACE_DIGEST)
    print("common_slope_ledgers=%s" % (block_ledgers,))
    print("raw_degree_ledger_p2=%s" % (degree_ledgers[0],))
    print("raw_degree_ledger_p5=%s" % (degree_ledgers[-1],))
    print("raw_degree_ledgers_sha256=%s" % digest(degree_ledgers))
    print("p5_carry=single_slope_312/625 cap=1250 denominator=625=5^4")
    print("positive_controls=primary(v+1,v); independent=%s" % independent_control)
    print("factorization_invoice=%s" % (factorization_invoice,))
    print("semantic_sha256=%s" % digest(semantic))
    print("STATUS=PASS closed=1 survivors=0 next_unaudited=d=2502/r=2500")


if __name__ == "__main__":
    main()
