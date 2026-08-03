#!/usr/bin/env python3
"""Dual exact extension audit for THM-3201 through d=2500."""

import ast
import hashlib
import importlib
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


PRIMARY_NAME = "factorial_seven_exit_first_flag_2400_thm3201.py"
PRIMARY_SHA256 = "1cf067eb3cc8da29e51ee8f560f4100fb66a4b9684a2da8b1f7ef14d380e4171"
INDEPENDENT_NAME = "factorial_seven_exit_first_flag_2400_independent_audit_thm3201.py"
INDEPENDENT_SHA256 = "afe3ec56c7efe99ff4a09c1af856b4e144312d48ec45e729d063a5fd83b0f0b0"
CHUNKS = ((2401, 2433), (2434, 2466), (2467, 2492), (2493, 2500))
EXPECTED_COUNTS = (90, 79, 68, 58, 49, 40, 35)
EXPECTED_RESIDUALS = (
    2406, 2407, 2408, 2409, 2410,
    2430, 2431, 2432, 2433, 2434, 2435, 2436,
    2454, 2455, 2456, 2457, 2458, 2466,
    2484, 2485, 2486, 2487, 2488, 2489, 2490, 2491, 2492,
    2493, 2494, 2495, 2496, 2497, 2498, 2499, 2500,
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def digest(value):
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def load_pinned(name, expected_hash):
    path = Path(__file__).resolve().with_name(name)
    require(path.is_file(), ("missing dependency", path))
    actual_hash = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual_hash == expected_hash, (path, actual_hash, expected_hash))
    return importlib.import_module(path.stem)


def flatten(results):
    records = tuple(record for result in results for record in result["records"])
    survivors = tuple(item for result in results for item in result["survivors"])
    return records, survivors


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256)
    independent = load_pinned(INDEPENDENT_NAME, INDEPENDENT_SHA256)

    base = primary.load_engine()
    companion = base.load_companion()
    independent_base = independent.load_engine()
    primary_counts, primary_residuals = primary.exit_counts(base, companion, 2401, 2500)
    independent_counts, independent_residuals = independent.exit_counts(
        independent_base, 2401, 2500
    )
    require(primary_counts == independent_counts == EXPECTED_COUNTS,
            (primary_counts, independent_counts))
    require(primary_residuals == independent_residuals == EXPECTED_RESIDUALS,
            (primary_residuals, independent_residuals))

    with ProcessPoolExecutor(max_workers=4) as pool:
        primary_results = list(pool.map(primary.scan_chunk, CHUNKS))
    with ProcessPoolExecutor(max_workers=4) as pool:
        independent_results = list(pool.map(independent.scan_chunk, CHUNKS))

    primary_records, primary_survivors = flatten(primary_results)
    independent_records, independent_survivors = flatten(independent_results)
    require(not primary_survivors and not independent_survivors,
            (primary_survivors, independent_survivors))
    require(primary_records == independent_records,
            ("semantic trace disagreement", digest(primary_records), digest(independent_records)))
    require(tuple(d for d, _ in primary_records) == EXPECTED_RESIDUALS,
            tuple(d for d, _ in primary_records))

    killers = tuple((d, trace[-1][0]) for d, trace in primary_records)
    killer_histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    require(killer_histogram == ((2, 1), (3, 5), (5, 9), (7, 12), (11, 6), (13, 2)),
            killer_histogram)

    next_primary = next(
        d for d in range(2501, 10000)
        if primary.seven_exit_residual(base, companion, d)
    )
    next_independent = next(
        d for d in range(2501, 10000)
        if independent.seven_exit_residual(independent_base, d)
    )
    require(next_primary == next_independent == 2501, (next_primary, next_independent))
    invoice = tuple(tuple(companion.factor(n)) for n in range(2501, 2494, -1))
    require(
        invoice
        == (
            ((41, 1), (61, 1)),
            ((2, 2), (5, 4)),
            ((3, 1), (7, 2), (17, 1)),
            ((2, 1), (1249, 1)),
            ((11, 1), (227, 1)),
            ((2, 6), (3, 1), (13, 1)),
            ((5, 1), (499, 1)),
        ),
        invoice,
    )

    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
            ("assert node", __file__))

    print("THM-3201 SEVEN-EXIT EXTENSION THROUGH d=2500 DUAL AUDIT")
    print("primary=%s sha256=%s" % (PRIMARY_NAME, PRIMARY_SHA256))
    print("independent=%s sha256=%s" % (INDEPENDENT_NAME, INDEPENDENT_SHA256))
    print("universe=2401<=d<=2500 census=%s" % (EXPECTED_COUNTS,))
    print("seven_exit_residuals=%d residual_digest=%s" %
          (len(EXPECTED_RESIDUALS), digest(EXPECTED_RESIDUALS)))
    print("balanced_chunks=%s" % (CHUNKS,))
    print("dual_semantic_trace_digest=%s" % digest(primary_records))
    print("closed=35 survivors=0")
    print("killer_histogram=%s max_killer_prime=%d" %
          (killer_histogram, max(prime for _, prime in killers)))
    print("first_next_seven_exit_residual=d=2501 r=2499")
    print("first_next_factorization=%s" % (invoice,))
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=2498 is nonzero in at least one slot")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
