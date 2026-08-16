#!/usr/bin/env python3
"""Dual exact fixed-bank replay and divisor-boundary audit for THM-3478.

The expensive default run independently reconstructs the resonant factorial
rows twice.  The primary route uses the canonical FLINT coefficient engine;
the independent route uses the separately pinned determinant/lower-hull
engine.  A cheap third stage compares the two proved THM-3475 digit-only
divisor compilers on the next six seven-exit rows.
"""

import argparse
import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path


START = 2501
END = 2600
CHUNKS = ((2501, 2516), (2517, 2565), (2566, 2575), (2576, 2600))
MAX_WORKERS = 4
BANK = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)

PRIMARY_NAME = "factorial_seven_exit_first_flag_2400_thm3201.py"
PRIMARY_SHA256 = "1cf067eb3cc8da29e51ee8f560f4100fb66a4b9684a2da8b1f7ef14d380e4171"
INDEPENDENT_NAME = (
    "factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py"
)
INDEPENDENT_SHA256 = (
    "afe3ec56c7efe99ff4a09c1af856b4e144312d48ec45e729d063a5fd83b0f0b0"
)
DIVISOR_PRIMARY_NAME = "factorial_all_divisor_digit_pair_compiler_thm3475.py"
DIVISOR_PRIMARY_SHA256 = (
    "834d0913eb5cd5b15684c7fb88af60e42d2a6ef36feb821e261c0498f55027ab"
)
DIVISOR_INDEPENDENT_NAME = (
    "factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py"
)
DIVISOR_INDEPENDENT_SHA256 = (
    "9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665"
)

EXPECTED_COUNTS = (89, 78, 69, 60, 52, 44, 38)
EXPECTED_RESIDUALS = (
    2501, 2502, 2510, 2511, 2512, 2513, 2514, 2515, 2516, 2517,
    2518, 2519, 2520, 2528, 2529, 2530, 2538, 2564, 2565, 2566,
    2567, 2568, 2569, 2570, 2571, 2572, 2573, 2574, 2575, 2576,
    2577, 2578, 2586, 2587, 2588, 2589, 2590, 2600,
)
EXPECTED_KILLERS = (
    (2501, 5), (2502, 5), (2510, 11), (2511, 7), (2512, 7),
    (2513, 3), (2514, 7), (2515, 3), (2516, 7), (2517, 3),
    (2518, 7), (2519, 3), (2520, 7), (2528, 7), (2529, 3),
    (2530, 7), (2538, 13), (2564, 13), (2565, 2), (2566, 13),
    (2567, 3), (2568, 17), (2569, 2), (2570, 13), (2571, 3),
    (2572, 17), (2573, 3), (2574, 13), (2575, 3), (2576, 13),
    (2577, 2), (2578, 13), (2586, 13), (2587, 3), (2588, 13),
    (2589, 3), (2590, 13), (2600, 7),
)
EXPECTED_HISTOGRAM = ((2, 3), (3, 11), (5, 2), (7, 9), (11, 1),
                      (13, 10), (17, 2))
EXPECTED_CHUNK_DIGESTS = (
    "d8ddf55e638c503d653927422899c4bddc381ac2bdcc0c53c9248bf161679da6",
    "f3575666ce9261b9186bbfe6bcc961dcd1600c484fe36e35e310084c45e99b32",
    "e58d91f4c2637b0555e079f5abc47703049f29e2c7e9258f792157d4df77e1cb",
    "8ec99dd199f66f764a163e8fb52d562b151a939a124e5a1efc8c4c6404fe95f2",
)
EXPECTED_GLOBAL_DIGEST = (
    "28bdc29a6dadfc941ccab7b8eddafd77bd6bbcec9a2574455d4ba3a6dd439b9f"
)
EXTENSION_ROWS = tuple(range(2601, 2607))
EXPECTED_EXTENSION_KILLERS = (
    (2601, 13), (2602, 17), (2603, 1301), (2604, 137), (2605, 7),
)
EXPECTED_EXTENSION_SURVIVORS = ((2606, (521, 1042, 1563, 2084)),)
EXPECTED_EXTENSION_DIGEST = (
    "602f0fac54c487114457683a3264d2a095a7f048f9a0b3769332d3ead0e61289"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def load_pinned(filename, expected_hash, module_name):
    source = Path(__file__).resolve().with_name(filename)
    require(source.is_file(), ("missing dependency", source))
    normalized = source.read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(normalized).hexdigest()
    require(actual_hash == expected_hash, (source, actual_hash, expected_hash))
    spec = importlib.util.spec_from_file_location(module_name, source)
    require(spec is not None and spec.loader is not None, ("bad import", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixed_module(label, suffix):
    if label == "primary":
        return load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "thm3478_primary_" + suffix)
    require(label == "independent", ("unknown engine", label))
    return load_pinned(
        INDEPENDENT_NAME,
        INDEPENDENT_SHA256,
        "thm3478_independent_" + suffix,
    )


def scan_job(job):
    label, chunk = job
    module = fixed_module(label, "%d_%d" % chunk)
    return module.scan_chunk(chunk)


def fixed_census(label):
    module = fixed_module(label, "controller")
    if label == "primary":
        engine = module.load_engine()
        companion = engine.load_companion()
        counts, residuals = module.exit_counts(engine, companion, START, END)
        planted, _, planted_zero, _ = companion.planted_control()
        require(1 in planted and 1 in planted_zero, ("primary controls", planted))
        control = "planted (v+1) and v factors retain degree one"
    else:
        engine = module.load_engine()
        counts, residuals = module.exit_counts(engine, START, END)
        control = engine.positive_controls()

    require(counts == EXPECTED_COUNTS, (label, counts))
    require(residuals == EXPECTED_RESIDUALS, (label, residuals))
    jobs = tuple((label, chunk) for chunk in CHUNKS)
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        results = tuple(pool.map(scan_job, jobs))

    records = tuple(record for result in results for record in result["records"])
    survivors = tuple(item for result in results for item in result["survivors"])
    chunk_digests = tuple(result["semantic_trace_digest"] for result in results)
    require(tuple(d for d, _ in records) == EXPECTED_RESIDUALS,
            (label, "row order"))
    require(not survivors, (label, "survivors", survivors))
    require(chunk_digests == EXPECTED_CHUNK_DIGESTS,
            (label, "chunk digests", chunk_digests))
    global_digest = digest(records)
    require(global_digest == EXPECTED_GLOBAL_DIGEST,
            (label, "global digest", global_digest))
    killers = tuple((d, trace[-1][0]) for d, trace in records)
    histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    require(killers == EXPECTED_KILLERS, (label, "killers", killers))
    require(histogram == EXPECTED_HISTOGRAM, (label, "histogram", histogram))
    return records, control


def load_divisor_compilers():
    primary = load_pinned(
        DIVISOR_PRIMARY_NAME,
        DIVISOR_PRIMARY_SHA256,
        "thm3478_divisor_primary",
    )
    independent = load_pinned(
        DIVISOR_INDEPENDENT_NAME,
        DIVISOR_INDEPENDENT_SHA256,
        "thm3478_divisor_independent",
    )
    coefficient = primary.load_pinned(
        primary.PRIMARY_NAME,
        primary.PRIMARY_SHA256,
        "thm3478_divisor_coefficient",
    )
    digital = primary.load_pinned(
        primary.DIGITAL_NAME,
        primary.DIGITAL_SHA256,
        "thm3478_divisor_digital",
    )
    return primary, independent, coefficient, digital


def primary_divisor_record(primary, coefficient, digital):
    records = []
    for d in EXTENSION_ROWS:
        size = d - 1
        possible = set(range(1, size))
        trace = []
        for prime, _ in coefficient.factor(size):
            _, _, blocks, local = primary.predicted_pair_compiler(
                coefficient, digital, size, prime
            )
            possible &= local
            normalized = tuple(
                (slope.numerator, slope.denominator, capacity)
                for slope, capacity, _ in blocks
            )
            trace.append((prime, normalized, tuple(sorted(possible))))
        records.append((d, tuple(trace), tuple(sorted(possible))))
    return tuple(records)


def independent_divisor_record(independent):
    records = []
    for d in EXTENSION_ROWS:
        size = d - 1
        possible = set(range(1, size))
        trace = []
        for prime in independent.prime_divisors(size):
            local, blocks, _, _ = independent.pair_barcode(size, prime)
            possible &= set(local)
            trace.append((prime, blocks, tuple(sorted(possible))))
        records.append((d, tuple(trace), tuple(sorted(possible))))
    return tuple(records)


def divisor_extension(fixed_primary, fixed_independent):
    primary, independent, coefficient, digital = load_divisor_compilers()
    primary_record = primary_divisor_record(primary, coefficient, digital)
    independent_record = independent_divisor_record(independent)
    require(primary_record == independent_record, (
        "divisor compiler disagreement",
        digest(primary_record),
        digest(independent_record),
    ))
    semantic_digest = digest(primary_record)
    require(semantic_digest == EXPECTED_EXTENSION_DIGEST,
            ("extension digest", semantic_digest))
    survivors = tuple((d, final) for d, _, final in primary_record if final)
    killers = tuple(
        (d, next(prime for prime, _, post in trace if not post))
        for d, trace, final in primary_record if not final
    )
    require(survivors == EXPECTED_EXTENSION_SURVIVORS, survivors)
    require(killers == EXPECTED_EXTENSION_KILLERS, killers)

    primary_engine = fixed_primary.load_engine()
    primary_coefficient = primary_engine.load_companion()
    independent_engine = fixed_independent.load_engine()
    seven_exit_checks = tuple(
        (
            d,
            fixed_primary.seven_exit_residual(
                primary_engine, primary_coefficient, d
            ),
            fixed_independent.seven_exit_residual(independent_engine, d),
        )
        for d in EXTENSION_ROWS
    )
    require(
        seven_exit_checks == tuple((d, True, True) for d in EXTENSION_ROWS),
        seven_exit_checks,
    )
    return primary_record


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fixed-engine",
        choices=("primary", "independent", "both"),
        default="both",
        help="select the expensive fixed-bank replay; default runs both routes",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    labels = ("primary", "independent") if args.fixed_engine == "both" else (
        args.fixed_engine,
    )
    results = {}
    controls = {}
    for label in labels:
        results[label], controls[label] = fixed_census(label)
    if args.fixed_engine == "both":
        require(results["primary"] == results["independent"],
                "primary/independent fixed records disagree")

    fixed_primary = fixed_module("primary", "extension_filter")
    fixed_independent = fixed_module("independent", "extension_filter")
    extension_record = divisor_extension(fixed_primary, fixed_independent)

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "assert node",
    )

    print("THM-3478 FACTORIAL DUAL SEVEN-EXIT BLOCK CLOSURE EXACT AUDIT")
    print("fixed_dependencies=%s:%s;%s:%s" % (
        PRIMARY_NAME, PRIMARY_SHA256, INDEPENDENT_NAME, INDEPENDENT_SHA256,
    ))
    print("divisor_dependencies=%s:%s;%s:%s" % (
        DIVISOR_PRIMARY_NAME, DIVISOR_PRIMARY_SHA256,
        DIVISOR_INDEPENDENT_NAME, DIVISOR_INDEPENDENT_SHA256,
    ))
    print("fixed_mode=%s implementations_completed=%s" % (args.fixed_engine, labels))
    print("fixed_universe=2501<=d<=2600 all seven-exit residuals after ordered filters")
    print("exit_order=(d prime,d-1 prime_power,d-2 prime,d-3 prime,d-4 prime,d-5 prime,d-6 prime)")
    print("fixed_counts=%s residual_count=%d residuals=%s" % (
        EXPECTED_COUNTS, len(EXPECTED_RESIDUALS), EXPECTED_RESIDUALS,
    ))
    print("fixed_bank=%s chunks=%s workers=%d" % (BANK, CHUNKS, MAX_WORKERS))
    print("fixed_schema=tuple((d,tuple((prime,sorted_post_intersection),...)),...)")
    print("fixed_serialization=json.dumps(payload,separators=(',',':'),sort_keys=True).encode('ascii')")
    print("fixed_chunk_semantic_sha256=%s" % (EXPECTED_CHUNK_DIGESTS,))
    print("fixed_global_semantic_sha256=%s" % EXPECTED_GLOBAL_DIGEST)
    print("fixed_killers=%s" % (EXPECTED_KILLERS,))
    print("fixed_killer_histogram=%s survivors=()" % (EXPECTED_HISTOGRAM,))
    print("fixed_positive_controls=%s" % (tuple(sorted(controls.items())),))
    print("fixed_status=FINITE-EXACT dual coefficient replay; not a structural all-height theorem")
    print("extension_universe=d=2601..2606 all six are seven-exit residuals")
    print("extension_schema=tuple((d,tuple((p,tuple((slope_num,slope_den,capacity),...),sorted_post),...),sorted_final),...)")
    print("extension_serialization=json.dumps(payload,separators=(',',':'),sort_keys=True).encode('ascii')")
    print("extension_dual_formula_records_agree=True semantic_sha256=%s" % (
        digest(extension_record),
    ))
    print("extension_killers=%s survivor_boundary=%s" % (
        EXPECTED_EXTENSION_KILLERS, EXPECTED_EXTENSION_SURVIVORS,
    ))
    print("extension_status=PROVED THM-3475 compiler + FINITE-EXACT six-row application")
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=2603 is nonzero in at least one slot")
    print("first_unclosed_by_this_block=r=2604 d=2606 divisor_degrees=(521,1042,1563,2084)")
    print("reference_runtime=primary_about_6300s;independent_about_4020s;each_4_workers;extension_under_1s_on_session_host")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
