#!/usr/bin/env python3
"""Exact primary six-exit/first-Euclidean-flag audit for THM-3180."""

import hashlib
import importlib.util
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


BANK = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
WORK_CHUNKS = ((2001, 2050), (2051, 2100), (2101, 2150), (2151, 2200))
OLD_REPORT_CHUNKS = (
    (1001, 1100),
    (1101, 1300),
    (1301, 1500),
    (1501, 1750),
    (1751, 2000),
)
MAX_WORKERS = 4
DEPENDENCY_NAME = "factorial_multiplace_newton_degree_barcode_thm3152.py"
DEPENDENCY_SHA256 = "f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0"


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def companion_path():
    companion = Path(__file__).resolve().with_name(DEPENDENCY_NAME)
    require(companion.is_file(), ("canonical THM-3152 companion not found", companion))
    return companion


def load_companion():
    source = companion_path()
    actual_hash = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual_hash == DEPENDENCY_SHA256, (source, actual_hash, DEPENDENCY_SHA256))
    spec = importlib.util.spec_from_file_location("thm3152_six_exit_primary", source)
    require(spec is not None and spec.loader is not None, ("bad import spec", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exit_predicates(c):
    return (
        lambda d: c.is_prime(d),
        lambda d: c.is_prime_power(d - 1),
        lambda d: c.is_prime(d - 2),
        lambda d: c.is_prime(d - 3),
        lambda d: c.is_prime(d - 4),
        lambda d: c.is_prime(d - 5),
    )


def exit_counts(c, start, end):
    alive = list(range(start, end + 1))
    counts = []
    for exits in exit_predicates(c):
        alive = [d for d in alive if not exits(d)]
        counts.append(len(alive))
    return tuple(counts), tuple(alive)


def six_exit_residual(c, d):
    return all(not exits(d) for exits in exit_predicates(c))


def seven_exit_residual(c, d):
    return six_exit_residual(c, d) and not c.is_prime(d - 6)


def flag_trace(c, rows):
    possible = set(range(1, min(row.degree() for row in rows) + 1))
    trace = []
    for prime in BANK:
        local, _ = c.degree_barcode_many(rows, prime)
        possible &= local
        trace.append((prime, tuple(sorted(possible))))
        if not possible:
            break
    return possible, tuple(trace)


def scan_chunk(chunk):
    start, end = chunk
    c = load_companion()
    rows = []
    residuals = []
    for d in range(start, end + 1):
        if not six_exit_residual(c, d):
            continue
        residuals.append(d)
        p, q = c.pair(d)
        r = c.first_full_pseudoremainder(p, q, d)
        require(
            (p.degree(), q.degree(), r.degree()) == (d - 2, d - 1, d - 3),
            (d, p.degree(), q.degree(), r.degree()),
        )
        possible, trace = flag_trace(c, (p, q, r))
        require(not possible, (d, tuple(sorted(possible)), "six-exit survivor"))
        rows.append((d, trace))
    require(tuple(residuals) == tuple(d for d, _ in rows), (chunk, residuals, rows))
    return {
        "chunk": chunk,
        "residuals": tuple(residuals),
        "rows": tuple(rows),
        "semantic_trace_digest": digest(rows),
    }


def main():
    c = load_companion()

    planted, _, planted_zero, _ = c.planted_control()
    require(1 in planted, ("planted nonzero factor lost", planted))
    require(1 in planted_zero, ("planted coordinate factor lost", planted_zero))

    expected_censuses = {
        (3, 1000): (831, 642, 513, 390, 301, 211),
        (1001, 2000): (865, 725, 617, 511, 426, 341),
        (3, 2000): (1696, 1367, 1130, 901, 727, 552),
        (2001, 2200): (176, 149, 130, 111, 95, 79),
    }
    censuses = {}
    residual_lists = {}
    for bounds, expected in expected_censuses.items():
        counts, residuals = exit_counts(c, *bounds)
        require(counts == expected, (bounds, counts, expected))
        censuses[bounds] = counts
        residual_lists[bounds] = residuals

    old_report_progression = tuple(
        exit_counts(c, start, end)[0][2:] for start, end in OLD_REPORT_CHUNKS
    )
    require(
        old_report_progression
        == ((56, 45, 38, 30), (123, 100, 79, 58), (121, 101, 86, 72),
            (154, 125, 104, 83), (163, 140, 119, 98)),
        old_report_progression,
    )
    new_chunk_progression = tuple(
        exit_counts(c, start, end)[0][2:] for start, end in WORK_CHUNKS
    )
    require(
        new_chunk_progression
        == ((31, 25, 20, 15), (29, 24, 20, 16),
            (31, 26, 22, 18), (39, 36, 33, 30)),
        new_chunk_progression,
    )

    seven_exit_rows = tuple(
        d for d in residual_lists[(2001, 2200)] if seven_exit_residual(c, d)
    )
    require(len(seven_exit_rows) == 66, len(seven_exit_rows))
    first_unaudited = next(
        d for d in range(2201, 10000) if seven_exit_residual(c, d)
    )
    require(first_unaudited == 2201, first_unaudited)
    require(
        (2201, 2200, 2199, 2198, 2197, 2196, 2195)
        == (31 * 71, 2**3 * 5**2 * 11, 3 * 733, 2 * 7 * 157,
            13**3, 2**2 * 3**2 * 61, 5 * 439),
        "first unaudited factorization invoice",
    )

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        chunk_results = list(pool.map(scan_chunk, WORK_CHUNKS))

    residuals = tuple(d for result in chunk_results for d in result["residuals"])
    rows = tuple(row for result in chunk_results for row in result["rows"])
    require(residuals == residual_lists[(2001, 2200)], (residuals, residual_lists[(2001, 2200)]))
    require(len(rows) == 79, len(rows))
    require(all(not trace[-1][1] for _, trace in rows), "a row did not close")

    killers = tuple((d, trace[-1][0]) for d, trace in rows)
    killer_histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    trace_2009 = next(trace for d, trace in rows if d == 2009)
    require(trace_2009[-1][1] == (), (2009, trace_2009))

    print("THM-3180 SIX-EXIT FIRST-EUCLIDEAN-FLAG PRIMARY AUDIT")
    print("canonical_dependency=%s sha256=%s" % (DEPENDENCY_NAME, DEPENDENCY_SHA256))
    print("exit_order=(d prime,d-1 prime power,d-2 prime,d-3 prime,d-4 prime,d-5 prime)")
    print("census_3_1000=%s" % (censuses[(3, 1000)],))
    print("census_1001_2000=%s" % (censuses[(1001, 2000)],))
    print("census_3_2000=%s" % (censuses[(3, 2000)],))
    print("old_report_progression_after_exits_3_to_6=%s" % (old_report_progression,))
    print("census_2001_2200=%s" % (censuses[(2001, 2200)],))
    print("new_chunk_progression_after_exits_3_to_6=%s" % (new_chunk_progression,))
    print("new_six_exit_residual_count=%d residual_digest=%s" % (len(residuals), digest(residuals)))
    print("related_thm3176_seven_exit_residual_count=%d" % len(seven_exit_rows))
    print("bank=%s" % (BANK,))
    print("deterministic_work_chunks=%s max_workers=%d" % (WORK_CHUNKS, MAX_WORKERS))
    print("chunk_semantic_trace_digests=%s" % (tuple(result["semantic_trace_digest"] for result in chunk_results),))
    print("global_semantic_trace_digest=%s" % digest(rows))
    print("killer_histogram=%s max_killer_prime=%d" % (killer_histogram, max(prime for _, prime in killers)))
    print("hostile_control_d_2009_trace=%s" % (trace_2009,))
    print("positive_controls=planted(v+1) and planted(v), degree 1 retained")
    print("closed=79 survivors=0")
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=2198 is nonzero in at least one slot")
    print("first_unaudited=r=2199 d=2201 seven_exit_residual=True")
    print("first_unaudited_factorization=(31*71,2^3*5^2*11,3*733,2*7*157,13^3,2^2*3^2*61,5*439)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
