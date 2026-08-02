#!/usr/bin/env python3
"""Exact primary companion for THM-3161.

It imports the proved THM-3152 arithmetic engine, scans the full finite
universe in deterministic parallel chunks, and prints theorem consequences
plus hashes of the exact narrowing traces.
"""

import hashlib
import importlib.util
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


REPORT_CHUNKS = (
    (1001, 1100),
    (1101, 1300),
    (1301, 1500),
    (1501, 1750),
    (1751, 2000),
)
WORK_CHUNKS = tuple(
    (start, min(start + 49, 2000)) for start in range(1001, 2001, 50)
)
MAX_WORKERS = 4
BANK = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
DEPENDENCY_NAME = "factorial_multiplace_newton_degree_barcode_thm3152.py"
DEPENDENCY_SHA256 = "f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0"


def companion_path():
    sibling = Path(__file__).resolve().with_name(DEPENDENCY_NAME)
    if sibling.exists():
        return sibling
    fallback = Path(
        "/tmp/math-wt-frontier-synthesis-20260802/04-computation/"
        "factorial_multiplace_newton_degree_barcode_thm3152.py"
    )
    if fallback.exists():
        return fallback
    raise RuntimeError("THM-3152 companion not found beside THM-3161")


def load_companion():
    source = companion_path()
    actual_hash = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    if actual_hash != DEPENDENCY_SHA256:
        raise RuntimeError((source, actual_hash, DEPENDENCY_SHA256))
    spec = importlib.util.spec_from_file_location("thm3152", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def three_exit_residual(c, d):
    return not c.is_prime(d) and not c.is_prime_power(d - 1) and not c.is_prime(d - 2)


def four_exit_residual(c, d):
    return three_exit_residual(c, d) and not c.is_prime(d - 3)


def trace_rows(c, rows):
    possible = set(range(1, min(row.degree() for row in rows) + 1))
    trace = []
    for prime in BANK:
        local, blocks = c.degree_barcode_many(rows, prime)
        narrowed = possible & local
        if narrowed != possible:
            trace.append(
                (
                    prime,
                    tuple(sorted(narrowed)),
                    tuple((str(s), cap, den) for s, cap, den in blocks),
                )
            )
        possible = narrowed
        if not possible:
            break
    return possible, tuple(trace)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def scan_chunk(chunk):
    start, end = chunk
    c = load_companion()
    three = []
    four = []
    current_records = []
    observer_survivors = []
    for d in range(start, end + 1):
        if not three_exit_residual(c, d):
            continue
        three.append(d)
        current = four_exit_residual(c, d)
        if current:
            four.append(d)
        # Every current row is load-bearing.  To locate the first observer
        # boundary, additionally replay every three-exit row through d=1384.
        if current or d <= 1384:
            p, q = c.pair(d)
            r = c.first_full_pseudoremainder(p, q, d)
            c.require(r.degree() == d - 3, (d, "remainder degree"))
            possible, trace = trace_rows(c, (p, q, r))
            if current:
                c.require(not possible, (d, possible, "current residual survivor"))
                current_records.append((d, trace))
            if d <= 1384 and possible:
                c.require(c.is_prime(d - 3), (d, possible, "survivor outside THM-3146"))
                observer_survivors.append((d, tuple(sorted(possible))))
    return {
        "chunk": chunk,
        "three": tuple(three),
        "four": tuple(four),
        "current_trace_digest": digest(current_records),
        "observer_survivors": tuple(observer_survivors),
    }


def boundary_1384():
    c = load_companion()
    d = 1384
    p, q = c.pair(d)
    r = c.first_full_pseudoremainder(p, q, d)
    possible, trace = trace_rows(c, (p, q, r))
    c.require(possible == {1}, (possible, "d=1384 fixed bank"))
    local173, blocks173 = c.degree_barcode_many((p, q, r), 173)
    c.require(not (possible & local173), (possible, local173, "p=173 must kill"))
    c.require(
        blocks173 == ((c.Fraction(1, 173), 1211, 173),),
        (blocks173, "p=173 block"),
    )
    return {
        "trace": tuple((prime, len(degrees), degrees) for prime, degrees, _ in trace),
        "fixed_final": (1,),
        "p173_positive": tuple(sorted(local173 - {0})),
        "p173_blocks": tuple((str(s), cap, den) for s, cap, den in blocks173),
    }


def main():
    c = load_companion()
    planted, _, planted_zero, _ = c.planted_control()
    c.require(1 in planted and 1 in planted_zero, "planted common factors")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        chunks = list(pool.map(scan_chunk, WORK_CHUNKS))
    three = tuple(d for result in chunks for d in result["three"])
    four = tuple(d for result in chunks for d in result["four"])
    survivors = tuple(
        record for result in chunks for record in result["observer_survivors"]
    )
    c.require(len(three) == 617, (len(three), "three-exit count"))
    c.require(len(four) == 511, (len(four), "four-exit count"))
    report_counts = tuple(
        sum(start <= d <= end for d in four) for start, end in REPORT_CHUNKS
    )
    c.require(report_counts == (45, 100, 101, 125, 140), report_counts)
    c.require(survivors == ((1384, (1,)),), (survivors, "first observer survivor"))
    boundary = boundary_1384()

    print("THM-3161 FOUR-EXIT FACTORIAL TAIL EXACT COMPANION")
    print("canonical_dependency=%s sha256=%s" % (DEPENDENCY_NAME, DEPENDENCY_SHA256))
    print("universe=1001<=d<=2000 inclusive")
    print("bank=%s" % (BANK,))
    print("three-exit residuals=%d" % len(three))
    print("current four-exit residuals=%d" % len(four))
    print("four-exit report-chunk counts=%s" % (report_counts,))
    print("current residual survivors=0")
    print("first three-exit observer survivor in scanned tail=%s" % (survivors[0],))
    print("all tail observer survivors, if any, have d-3 prime=True (by current closure)")
    print("work-chunk exact trace digests=%s" % (tuple(result["current_trace_digest"] for result in chunks),))
    print("d=1384 progressive trace=%s" % (boundary["trace"],))
    print("d=1384 fixed final=%s" % (boundary["fixed_final"],))
    print("d=1384 p=173 blocks=%s" % (boundary["p173_blocks"],))
    print("d=1384 p=173 positive degrees=%s" % (boundary["p173_positive"],))
    print("planted v+1 retained=%s planted v retained=%s" % (planted, planted_zero))
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=1998 is nonzero in at least one slot")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
