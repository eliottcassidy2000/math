#!/usr/bin/env python3
"""Independent determinant-lower-hull audit for THM-3201."""

import hashlib
import importlib.util
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


BANK = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
WORK_CHUNKS = ((2201, 2250), (2251, 2300), (2301, 2350), (2351, 2400))
MAX_WORKERS = 4
ENGINE_NAME = "factorial_six_exit_first_flag_2200_independent_audit_thm3180.py"
ENGINE_SHA256 = "5d00729ae39cdf1d085eca73d4612d4960a0788dd746f945ab2e8936ebbd466d"


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def engine_path():
    engine = Path(__file__).resolve().with_name(ENGINE_NAME)
    require(engine.is_file(), ("THM-3180 independent engine not found", engine))
    return engine


def load_engine():
    source = engine_path()
    actual_hash = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual_hash == ENGINE_SHA256, (source, actual_hash, ENGINE_SHA256))
    spec = importlib.util.spec_from_file_location("thm3180_independent_for_thm3201", source)
    require(spec is not None and spec.loader is not None, ("bad import spec", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def seven_exit_predicates(a):
    return a.exit_predicates() + (lambda d: a.prime(d - 6),)


def exit_counts(a, start, end):
    alive = list(range(start, end + 1))
    counts = []
    for exits in seven_exit_predicates(a):
        alive = [d for d in alive if not exits(d)]
        counts.append(len(alive))
    return tuple(counts), tuple(alive)


def seven_exit_residual(a, d):
    return all(not exits(d) for exits in seven_exit_predicates(a))


def scan_chunk(chunk):
    start, end = chunk
    a = load_engine()
    records = []
    survivors = []
    for d in range(start, end + 1):
        if not seven_exit_residual(a, d):
            continue
        p, q = a.moments(d)
        r = a.first_full_remainder(p, q, d)
        require(
            (p.degree(), q.degree(), r.degree()) == (d - 2, d - 1, d - 3),
            (d, p.degree(), q.degree(), r.degree()),
        )
        possible, trace = a.flag_trace((p, q, r))
        records.append((d, trace))
        if possible:
            survivors.append((d, tuple(sorted(possible))))
    return {
        "chunk": chunk,
        "records": tuple(records),
        "survivors": tuple(survivors),
        "semantic_trace_digest": digest(records),
    }


def main():
    a = load_engine()
    control_message = a.positive_controls()

    counts, residuals = exit_counts(a, 2201, 2400)
    require(counts == (170, 139, 115, 92, 75, 58, 50), counts)
    chunk_progressions = tuple(
        exit_counts(a, start, end)[0] for start, end in WORK_CHUNKS
    )
    require(
        chunk_progressions
        == ((43, 35, 29, 24, 20, 16, 13),
            (42, 34, 27, 20, 16, 12, 10),
            (44, 38, 34, 30, 26, 22, 20),
            (41, 32, 25, 18, 13, 8, 7)),
        chunk_progressions,
    )

    first_next = next(d for d in range(2401, 10000) if seven_exit_residual(a, d))
    require(first_next == 2406, first_next)
    require(
        (2406, 2405, 2404, 2403, 2402, 2401, 2400)
        == (2 * 3 * 401, 5 * 13 * 37, 2**2 * 601, 3**3 * 89,
            2 * 1201, 7**4, 2**5 * 3 * 5**2),
        "first next residual factorization invoice",
    )

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        chunk_results = list(pool.map(scan_chunk, WORK_CHUNKS))

    records = tuple(record for result in chunk_results for record in result["records"])
    scanned = tuple(d for d, _ in records)
    survivors = tuple(item for result in chunk_results for item in result["survivors"])
    require(scanned == residuals, (scanned, residuals))
    require(len(records) == 50, len(records))
    require(not survivors, survivors)

    closed_records = tuple((d, trace) for d, trace in records if not trace[-1][1])
    require(len(closed_records) == 50, len(closed_records))
    killers = tuple((d, trace[-1][0]) for d, trace in closed_records)
    killer_histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    trace_2201 = next(trace for d, trace in records if d == 2201)
    require(trace_2201 == ((2, (2048,)), (3, ())), trace_2201)

    print("THM-3201 SEVEN-EXIT FIRST-EUCLIDEAN-FLAG INDEPENDENT AUDIT")
    print("independent_engine=%s sha256=%s" % (ENGINE_NAME, ENGINE_SHA256))
    print("implementation=self-contained recurrence + trial division + determinant lower hull")
    print("exit_order=(d prime,d-1 prime power,d-2 prime,d-3 prime,d-4 prime,d-5 prime,d-6 prime)")
    print("universe=2201<=d<=2400 census=%s" % (counts,))
    print("chunk_progressions=%s" % (chunk_progressions,))
    print("seven_exit_residual_count=%d residual_digest=%s" % (len(residuals), digest(residuals)))
    print("bank=%s" % (BANK,))
    print("deterministic_work_chunks=%s max_workers=%d" % (WORK_CHUNKS, MAX_WORKERS))
    print("chunk_semantic_trace_digests=%s" % (tuple(result["semantic_trace_digest"] for result in chunk_results),))
    print("global_semantic_trace_digest=%s" % digest(records))
    print("closed=50 survivors=0")
    print("survivor_rows=()")
    print("killer_histogram=%s max_killer_prime=%d" % (killer_histogram, max(prime for _, prime in killers)))
    print("hostile_boundary_d_2201_trace=%s" % (trace_2201,))
    print("positive_controls=%s" % control_message)
    print("first_next_seven_exit_residual=d=2406 r=2404")
    print("first_next_factorization=(2*3*401,5*13*37,2^2*601,3^3*89,2*1201,7^4,2^5*3*5^2)")
    print("consequence=every exact {0,1,2} quadratic window 1<=r<=2403 is nonzero in at least one slot")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
