#!/usr/bin/env python3
"""Exact dual-engine companion for THM-3474."""

import ast
import hashlib
import importlib.util
import json
from pathlib import Path


PRIMARY_NAME = "factorial_multiplace_newton_degree_barcode_thm3152.py"
PRIMARY_SHA256 = "f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0"
INDEPENDENT_NAME = "factorial_six_exit_first_flag_2200_independent_audit_thm3180.py"
INDEPENDENT_SHA256 = "5d00729ae39cdf1d085eca73d4612d4960a0788dd746f945ab2e8936ebbd466d"
BINARY_DIGEST = "888b5d72e72bfeae0ef53140420e728fc5c7a78bfab60e6ac29e694a0ee60edd"
RESET_DIGEST = "3512b41ae23e168e5c500a989f37b6120e2dee7670c700b1f40f32b2fe3c45e0"

RESET_CASES = (
    (5, 2, 1),
    (5, 2, 2),
    (5, 2, 3),
    (5, 4, 1),
    (5, 4, 2),
    (5, 4, 3),
    (7, 2, 1),
    (7, 2, 2),
    (11, 2, 1),
    (11, 4, 1),
    (11, 4, 2),
    (17, 2, 1),
    (17, 4, 1),
    (17, 8, 1),
    (17, 16, 1),
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":")).encode("ascii")
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


def predicted_binary_hull(even_size):
    prefix = 0
    count = 0
    answer = [(0, 0)]
    for exponent in range(even_size.bit_length()):
        if (even_size >> exponent) & 1:
            count += 1
            prefix += 1 << exponent
            answer.append((prefix, 2 * prefix - count))
    return tuple(answer)


def independent_hull(independent, polynomial, prime):
    answer = []
    for x_coordinate in range(polynomial.degree() + 1):
        y_coordinate = independent.valuation(polynomial[x_coordinate], prime)
        if y_coordinate is None:
            continue
        while len(answer) >= 2:
            x0, y0 = answer[-2]
            x1, y1 = answer[-1]
            if ((y1 - y0) * (x_coordinate - x1)
                    < (y_coordinate - y1) * (x1 - x0)):
                break
            answer.pop()
        answer.append((x_coordinate, y_coordinate))
    return tuple(answer)


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "thm3474_primary")
    independent = load_pinned(
        INDEPENDENT_NAME,
        INDEPENDENT_SHA256,
        "thm3474_independent",
    )

    binary_records = []
    for even_size in range(2, 121, 2):
        expected_hull = predicted_binary_hull(even_size)
        _, primary_g = primary.pair(even_size + 1)
        _, independent_g = independent.moments(even_size + 1)
        primary_hull = tuple(primary.lower_hull(primary_g, 2))
        referee_hull = independent_hull(independent, independent_g, 2)
        require(primary_hull == expected_hull, (even_size, "primary hull", primary_hull))
        require(referee_hull == expected_hull, (even_size, "independent hull", referee_hull))
        submasks = {
            degree
            for degree in range(even_size + 1)
            if degree & ~even_size == 0
        }
        primary_degrees, _ = primary.degree_barcode_many((primary_g,), 2)
        referee_degrees = independent.allowed_degrees((independent_g,), 2)
        require(primary_degrees == submasks, (
            even_size,
            "primary submasks",
            tuple(sorted(primary_degrees)),
        ))
        require(referee_degrees == submasks, (
            even_size,
            "independent submasks",
            tuple(sorted(referee_degrees)),
        ))
        binary_records.append((even_size, expected_hull))
    require(digest(binary_records) == BINARY_DIGEST, digest(binary_records))

    reset_records = []
    for prime, multiplier, depth in RESET_CASES:
        height = prime**depth
        size = multiplier * height
        d = size + 1
        primary_f, primary_g = primary.pair(d)
        primary_e = primary.first_full_pseudoremainder(primary_f, primary_g, d)
        primary_rows = (primary_f, primary_g, primary_e)
        primary_binary = primary.degree_barcode_many(primary_rows, 2)[0]
        primary_reset = primary.degree_barcode_many(primary_rows, prime)[0]

        referee_f, referee_g = independent.moments(d)
        referee_e = independent.first_full_remainder(referee_f, referee_g, d)
        referee_rows = (referee_f, referee_g, referee_e)
        referee_binary = independent.allowed_degrees(referee_rows, 2)
        referee_reset = independent.allowed_degrees(referee_rows, prime)
        require(primary_binary == referee_binary, (
            prime,
            multiplier,
            depth,
            "binary engine disagreement",
        ))
        require(primary_reset == referee_reset, (
            prime,
            multiplier,
            depth,
            "reset engine disagreement",
        ))
        positive_intersection = (primary_binary & primary_reset) - {0}
        require(not positive_intersection, (
            prime,
            multiplier,
            depth,
            tuple(sorted(positive_intersection)),
        ))
        reset_records.append((
            prime,
            multiplier,
            depth,
            d,
            tuple(sorted(primary_binary)),
            tuple(sorted(primary_reset)),
        ))
    require(digest(reset_records) == RESET_DIGEST, digest(reset_records))

    formal_odd_hull = predicted_binary_hull(5)
    _, primary_odd = primary.pair(6)
    _, referee_odd = independent.moments(6)
    actual_odd_hull = tuple(primary.lower_hull(primary_odd, 2))
    referee_odd_hull = independent_hull(independent, referee_odd, 2)
    require(actual_odd_hull == ((0, 3), (5, 8)), actual_odd_hull)
    require(referee_odd_hull == actual_odd_hull, referee_odd_hull)
    require(actual_odd_hull != formal_odd_hull, (formal_odd_hull, actual_odd_hull))

    planted, _, planted_zero, _ = primary.planted_control()
    require(1 in planted and 1 in planted_zero, (planted, planted_zero))
    referee_controls = independent.positive_controls()

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "assert node",
    )
    semantic = (
        tuple(binary_records),
        tuple(reset_records),
        formal_odd_hull,
        actual_odd_hull,
        referee_controls,
    )

    print("THM-3474 FACTORIAL BINARY-SUBMASK RESET FAMILIES EXACT AUDIT")
    print("dependencies=%s:%s;%s:%s" % (
        PRIMARY_NAME,
        PRIMARY_SHA256,
        INDEPENDENT_NAME,
        INDEPENDENT_SHA256,
    ))
    print("binary_universe=even_2<=N<=120 cells=60 dual_hulls_and_submasks=PASS")
    print("binary_semantic_sha256=%s" % BINARY_DIGEST)
    print("reset_cases=%s" % (RESET_CASES,))
    print("reset_cells=15 positive_intersections=0 dual_engines=PASS")
    print("reset_semantic_sha256=%s" % RESET_DIGEST)
    print("hostile_odd_N=5 formal=%s actual=%s" % (formal_odd_hull, actual_odd_hull))
    print("positive_controls=primary(v+1,v); independent=%s" % referee_controls)
    print("consequence=power_two_reset_families_verified_in_15_cells")
    print("semantic_sha256=%s" % digest(semantic))
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
