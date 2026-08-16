#!/usr/bin/env python3
"""Exact all-divisor digit-polygon and pair-ledger audit for THM-3475."""

import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path


PRIMARY_NAME = "factorial_multiplace_newton_degree_barcode_thm3152.py"
PRIMARY_SHA256 = "f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0"
DIGITAL_NAME = "factorial_predecessor_digital_skeleton_audit_thm3161.py"
DIGITAL_SHA256 = "0a9da30ffd1462bc2f6ca3d62ad4e7d71b82afe02099931542928105704fd318"
ODD_DIGEST = "1fe62de92a499ecaa9eefd4aaca90ae77330f4857fd33103f15f9b848c116142"
BINARY_DIGEST = "c1a56e649a584069678c9e275a8a28b69d80e359f27cbc64aef36f078daaeed8"
PAIR_DIGEST = "02915c06dcb6569bbd7a19193e38b964a4e35320ec3368ff6177958492d5c202"
CENSUS_DIGEST = "27c1c59959a58054a3420d0fed7944da7a3ce480ad91a00e9e753862af83efc7"
INDEPENDENT_CENSUS_DIGEST = (
    "e019eb61019620cfaffa8b5bb5769e8d171f08fb5ffeb2153f209c0128d42115"
)
SEMANTIC_DIGEST = "886f2bcd66711a44668b003717dd4f39643fc0e9cb0b694708b8670eeaf21499"

CENSUS_RESIDUALS = (
    2501, 2502, 2510, 2511, 2512, 2513, 2514, 2515, 2516, 2517,
    2518, 2519, 2520, 2528, 2529, 2530, 2538, 2564, 2565, 2566,
    2567, 2568, 2569, 2570, 2571, 2572, 2573, 2574, 2575, 2576,
    2577, 2578, 2586, 2587, 2588, 2589, 2590, 2600,
)

EXPECTED_CENSUS_SURVIVORS = (
    (2516, (503, 1006, 1509, 2012)),
    (2564, (466, 699, 1165, 1631, 1864, 2097, 2330)),
    (2571, (2056,)),
    (2576, (
        103, 206, 309, 412, 515, 618, 721, 824, 927, 1030, 1133, 1236,
        1339, 1442, 1545, 1648, 1751, 1854, 1957, 2060, 2163, 2266,
        2369, 2472,
    )),
    (2586, (
        47, 141, 188, 235, 282, 329, 2209, 2256, 2303, 2350, 2397,
        2444, 2491, 2538,
    )),
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
    normalized = source.read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(normalized).hexdigest()
    require(actual_hash == expected_hash, (source, actual_hash, expected_hash))
    spec = importlib.util.spec_from_file_location(module_name, source)
    require(spec is not None and spec.loader is not None, ("bad import spec", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def factorial_valuation(number, prime):
    answer = 0
    while number:
        number //= prime
        answer += number
    return answer


def integer_valuation(number, prime):
    require(number > 0, ("positive valuation input", number, prime))
    answer = 0
    while number % prime == 0:
        number //= prime
        answer += 1
    return answer


def f_weight(size, index, prime):
    """Raw Kummer--Legendre weight for A_(size-1)^(size+1)."""
    return (
        factorial_valuation(size - 1, prime)
        - factorial_valuation(index, prime)
        - factorial_valuation(size - 1 - index, prime)
        + factorial_valuation(2 * index, prime)
    )


def f_z_mod_prime(size, index, prime):
    """Residual factor Z for A_(size-1)^(size+1), reduced modulo prime."""
    remainder_length = size - 1 - index
    answer = 0
    rising = 1
    base = (size + 1) % prime
    for ell in range(remainder_length + 1):
        if ell:
            rising = rising * (2 * index + ell) % prime
        term = (
            comb(remainder_length, ell)
            * pow(base, remainder_length - ell, prime)
            * rising
        )
        answer = (answer - term if ell % 2 else answer + term) % prime
    return answer


def slope_ledger(points):
    answer = {}
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        slope = Fraction(y1 - y0, x1 - x0)
        require(slope not in answer, ("repeated maximal slope", slope, points))
        answer[slope] = x1 - x0
    return answer


def degrees_from_blocks(blocks):
    answer = {0}
    for slope, capacity, denominator in blocks:
        require(slope.denominator == denominator, (slope, denominator))
        require(capacity % denominator == 0, (slope, capacity, denominator))
        answer = {
            old + used
            for old in answer
            for used in range(0, capacity + 1, denominator)
        }
    return answer


def predicted_pair_compiler(primary, digital, size, prime):
    g_hull = tuple(digital.hull(
        tuple((index, digital.weight(size, index, prime))
              for index in range(size + 1))
    ))
    if prime == 2:
        f_points = tuple(
            (index, f_weight(size, index, prime))
            for index in range(1, size, 2)
        )
    else:
        f_points = tuple(
            (index, f_weight(size, index, prime))
            for index in range(size)
        )
    f_hull = tuple(digital.hull(f_points))
    f_ledger = slope_ledger(f_hull)
    g_ledger = slope_ledger(g_hull)
    blocks = tuple(
        (slope, min(f_ledger[slope], g_ledger[slope]), slope.denominator)
        for slope in sorted(set(f_ledger) & set(g_ledger))
        if slope > 0
    )
    return f_hull, g_hull, blocks, degrees_from_blocks(blocks)


def odd_audit(primary, digital):
    records = []
    pair_records = []
    for size in range(3, 161):
        odd_primes = tuple(prime for prime, _ in primary.factor(size) if prime % 2)
        if not odd_primes:
            continue
        actual_f, actual_g = primary.pair(size + 1)
        for prime in odd_primes:
            half = (prime - 1) // 2
            quotient = size // prime
            raw_hull = tuple(digital.hull(
                tuple((index, f_weight(size, index, prime))
                      for index in range(size))
            ))
            raw_suffix = tuple(point for point in raw_hull if point[0] >= half)
            actual_hull = tuple(primary.lower_hull(actual_f, prime))
            actual_suffix = tuple(point for point in actual_hull if point[0] >= half)
            require(raw_suffix == actual_suffix, (
                size,
                prime,
                "odd suffix hull",
                raw_suffix,
                actual_suffix,
            ))

            for q in range(quotient):
                base_height = 2 * q + f_weight(quotient, q, prime)
                for residue in range(prime):
                    index = prime * q + residue
                    predicted = base_height
                    if residue > half:
                        predicted += 1 + integer_valuation(2 * q + 1, prime)
                    require(f_weight(size, index, prime) == predicted, (
                        size,
                        prime,
                        index,
                        "odd recurrence",
                    ))
                    require(
                        digital.weight(size - 1, index, prime)
                        == f_weight(size, index, prime),
                        (size, prime, index, "weight implementation disagreement"),
                    )
                for residue in (half, prime - 1):
                    index = prime * q + residue
                    require(f_z_mod_prime(size, index, prime) == 1, (
                        size,
                        prime,
                        index,
                        "odd anchor Z",
                    ))
                    require(primary.vp(actual_f[index], prime) == f_weight(
                        size,
                        index,
                        prime,
                    ), (size, prime, index, "odd anchor coefficient"))

            predicted_f, predicted_g, predicted_blocks, predicted_degrees = (
                predicted_pair_compiler(primary, digital, size, prime)
            )
            require(predicted_f == raw_hull, (size, prime, "predicted F hull"))
            actual_g_hull = tuple(primary.lower_hull(actual_g, prime))
            require(predicted_g == actual_g_hull, (
                size,
                prime,
                "inherited G hull",
                predicted_g,
                actual_g_hull,
            ))
            actual_degrees, actual_blocks = primary.degree_barcode_many(
                (actual_f, actual_g),
                prime,
            )
            require(predicted_blocks == actual_blocks, (
                size,
                prime,
                "odd pair blocks",
                predicted_blocks,
                actual_blocks,
            ))
            require(predicted_degrees == actual_degrees, (
                size,
                prime,
                "odd pair degrees",
            ))
            records.append((size, prime, raw_suffix, actual_suffix))
            pair_records.append((
                size,
                prime,
                tuple((str(slope), capacity, denominator)
                      for slope, capacity, denominator in predicted_blocks),
                tuple(sorted(predicted_degrees)),
            ))
    require(len(records) == 209, ("odd pair count", len(records)))
    require(digest(tuple(records)) == ODD_DIGEST, digest(tuple(records)))
    return tuple(records), tuple(pair_records)


def binary_audit(primary, digital):
    records = []
    pair_records = []
    zero_coefficients = []
    for size in range(2, 201, 2):
        actual_f, actual_g = primary.pair(size + 1)
        odd_anchors = tuple(
            (index, f_weight(size, index, 2))
            for index in range(1, size, 2)
        )
        predicted_suffix = tuple(digital.hull(odd_anchors))
        actual_suffix = tuple(
            point for point in primary.lower_hull(actual_f, 2) if point[0] >= 1
        )
        require(predicted_suffix == actual_suffix, (
            size,
            "binary suffix hull",
            predicted_suffix,
            actual_suffix,
        ))
        quotient = size // 2
        for q in range(quotient):
            base_height = 2 * q + f_weight(quotient, q, 2)
            require(f_weight(size, 2 * q, 2) == base_height, (
                size,
                q,
                "binary even recurrence",
            ))
            require(f_weight(size, 2 * q + 1, 2) == base_height + 1, (
                size,
                q,
                "binary odd recurrence",
            ))
        for index in range(size):
            actual_height = primary.vp(actual_f[index], 2)
            if index % 2:
                require(f_z_mod_prime(size, index, 2) == 1, (
                    size,
                    index,
                    "binary odd Z",
                ))
                require(actual_height == f_weight(size, index, 2), (
                    size,
                    index,
                    "binary odd coefficient",
                ))
            else:
                require(f_z_mod_prime(size, index, 2) == 0, (
                    size,
                    index,
                    "binary even Z",
                ))
                if actual_height is None:
                    zero_coefficients.append((size, index))
                else:
                    require(actual_height >= f_weight(size, index, 2) + 1, (
                        size,
                        index,
                        "binary even coefficient",
                        actual_height,
                    ))

        predicted_f, predicted_g, predicted_blocks, predicted_degrees = (
            predicted_pair_compiler(primary, digital, size, 2)
        )
        require(predicted_f == predicted_suffix, (size, "binary predicted F hull"))
        actual_g_hull = tuple(primary.lower_hull(actual_g, 2))
        require(predicted_g == actual_g_hull, (
            size,
            "binary inherited G hull",
            predicted_g,
            actual_g_hull,
        ))
        actual_degrees, actual_blocks = primary.degree_barcode_many(
            (actual_f, actual_g),
            2,
        )
        require(predicted_blocks == actual_blocks, (
            size,
            "binary pair blocks",
            predicted_blocks,
            actual_blocks,
        ))
        require(predicted_degrees == actual_degrees, (
            size,
            "binary pair degrees",
        ))
        records.append((size, predicted_suffix))
        pair_records.append((
            size,
            2,
            tuple((str(slope), capacity, denominator)
                  for slope, capacity, denominator in predicted_blocks),
            tuple(sorted(predicted_degrees)),
        ))
    require(len(records) == 100, ("binary cell count", len(records)))
    require(digest(tuple(records)) == BINARY_DIGEST, digest(tuple(records)))
    require(tuple(zero_coefficients) == ((4, 2),), zero_coefficients)
    return tuple(records), tuple(pair_records), tuple(zero_coefficients)


def divisor_census(primary, digital):
    """Digit-only divisor-prime pair census on the exact 38-row universe.

    The hash payload is the tuple of records

      (d,
       tuple((p, tuple((str(slope), capacity, denominator), ...),
              post_intersection_cardinality, post_intersection_preview), ...),
       tuple(sorted(final_positive_degrees))).

    The preview is the complete sorted tuple at cardinality at most 30, and
    otherwise ``(first_12, last_12)`` as a pair of tuples.

    It is serialized by ``json.dumps(payload, separators=(",", ":"),
    sort_keys=True).encode("ascii")`` before SHA-256.
    """
    records = []
    for d in CENSUS_RESIDUALS:
        size = d - 1
        possible = set(range(1, size))
        trace = []
        for prime, _ in primary.factor(size):
            _, _, blocks, local_degrees = predicted_pair_compiler(
                primary,
                digital,
                size,
                prime,
            )
            possible &= local_degrees
            serialized_blocks = tuple(
                (str(slope), capacity, denominator)
                for slope, capacity, denominator in blocks
            )
            ordered = tuple(sorted(possible))
            preview = ordered if len(ordered) <= 30 else (ordered[:12], ordered[-12:])
            trace.append((prime, serialized_blocks, len(possible), preview))
        records.append((d, tuple(trace), tuple(sorted(possible))))
    records = tuple(records)
    require(digest(records) == CENSUS_DIGEST, digest(records))
    survivors = tuple((d, degrees) for d, _, degrees in records if degrees)
    require(survivors == EXPECTED_CENSUS_SURVIVORS, survivors)
    closed = tuple(d for d, _, degrees in records if not degrees)
    require(len(closed) == 33, ("closed census rows", len(closed)))
    return records, closed, survivors


def hostile_controls(primary, digital):
    hostile_g_hull = tuple(digital.hull(
        tuple((index, digital.weight(6, index, 3)) for index in range(7))
    ))
    require(hostile_g_hull == ((0, 0), (3, 2), (6, 5)), hostile_g_hull)
    require(hostile_g_hull != ((0, 0), (6, 4)), hostile_g_hull)
    _, _, hostile_pair_blocks, _ = predicted_pair_compiler(primary, digital, 6, 3)
    require(
        hostile_pair_blocks == ((Fraction(2, 3), 3, 3),),
        hostile_pair_blocks,
    )
    hostile_pair_serialized = tuple(
        (str(slope), capacity, denominator)
        for slope, capacity, denominator in hostile_pair_blocks
    )

    nondivisor_actual = (
        factorial_valuation(4, 3)
        - factorial_valuation(2, 3)
        - factorial_valuation(2, 3)
    )
    nondivisor_false_prediction = (
        factorial_valuation(0, 3)
        - factorial_valuation(0, 3)
        - factorial_valuation(0, 3)
    )
    require(nondivisor_actual == 1, nondivisor_actual)
    require(nondivisor_false_prediction == 0, nondivisor_false_prediction)

    raw = tuple((index, f_weight(6, index, 3)) for index in range(6))
    raw_hull = tuple(digital.hull(raw))
    raised = tuple((index, 2 if index == 0 else height) for index, height in raw)
    raised_hull = tuple(digital.hull(raised))
    require(raw_hull != raised_hull, (raw_hull, raised_hull))
    require(
        tuple(point for point in raw_hull if point[0] >= 1)
        == tuple(point for point in raised_hull if point[0] >= 1),
        (raw_hull, raised_hull),
    )

    planted, _, planted_zero, _ = primary.planted_control()
    require(1 in planted and 1 in planted_zero, (planted, planted_zero))
    return (
        hostile_g_hull,
        hostile_pair_serialized,
        (nondivisor_actual, nondivisor_false_prediction),
        raw_hull,
        raised_hull,
        planted,
        planted_zero,
    )


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "thm3475_primary")
    digital = load_pinned(DIGITAL_NAME, DIGITAL_SHA256, "thm3475_digital")
    odd_records, odd_pair_records = odd_audit(primary, digital)
    binary_records, binary_pair_records, zero_coefficients = binary_audit(
        primary,
        digital,
    )
    census_records, census_closed, census_survivors = divisor_census(
        primary,
        digital,
    )
    hostiles = hostile_controls(primary, digital)
    pair_records = odd_pair_records + binary_pair_records
    pair_digest = digest(pair_records)
    require(pair_digest == PAIR_DIGEST, pair_digest)

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "assert node",
    )
    semantic = (
        odd_records,
        binary_records,
        pair_records,
        census_records,
        zero_coefficients,
        hostiles,
    )
    semantic_digest = digest(semantic)
    require(semantic_digest == SEMANTIC_DIGEST, semantic_digest)

    print("THM-3475 FACTORIAL ALL-DIVISOR DIGIT PAIR COMPILER EXACT AUDIT")
    print("dependencies=%s:%s;%s:%s" % (
        PRIMARY_NAME,
        PRIMARY_SHA256,
        DIGITAL_NAME,
        DIGITAL_SHA256,
    ))
    print("odd_universe=3<=N<=160 every odd prime p|N pairs=209")
    print("odd_recurrences=PASS anchor_Z_units=PASS actual_suffix_hulls=PASS")
    print("odd_semantic_sha256=%s" % ODD_DIGEST)
    print("binary_universe=even 2<=N<=200 cells=100")
    print("binary_parity_Z=PASS odd_anchor_hulls=PASS zero_coefficients=%s" % (
        zero_coefficients,
    ))
    print("binary_semantic_sha256=%s" % BINARY_DIGEST)
    print("pair_ledgers=309 digit_vs_actual_blocks_and_degrees=PASS")
    print("pair_semantic_sha256=%s" % pair_digest)
    print("census_universe=38 seven-exit residuals 2501<=d<=2600")
    print("census_schema=(d,((p,((str_slope,cap,den),...),post_count,preview),...),final_degrees)")
    print("census_preview=full_sorted_tuple_if_count<=30_else_(first12,last12)")
    print("census_serialization=json.dumps(tuple(records),separators=(',',':'),sort_keys=True).encode('ascii')")
    print("census_closed=33 survivors=%s" % (census_survivors,))
    print("census_semantic_sha256=%s independent_sha256=%s" % (
        CENSUS_DIGEST,
        INDEPENDENT_CENSUS_DIGEST,
    ))
    print("hostile_unclipped_N6_p3_G_hull=%s pair_blocks=%s" % (
        hostiles[0],
        hostiles[1],
    ))
    print("hostile_nondivisor_N5_p3_binomial=(actual,predicted)=%s" % (
        hostiles[2],
    ))
    print("hostile_raised_prefix_raw=%s raised=%s suffix_preserved=PASS" % (
        hostiles[3],
        hostiles[4],
    ))
    print("positive_controls=planted(v+1,v)_degree_one_retained")
    print("semantic_sha256=%s" % semantic_digest)
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
