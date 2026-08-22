#!/usr/bin/env python3
"""Independent exact audit for THM-3161 using the THM-3152 referee engine."""

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
DEPENDENCY_NAME = "factorial_multiplace_newton_degree_barcode_independent_audit_thm3152.py"
DEPENDENCY_SHA256 = "427ba0e8c7e5b25efdf56b00200ab201b534f150c6702d3d0ba4a1ce956f32de"


def engine_path():
    sibling = Path(__file__).resolve().with_name(DEPENDENCY_NAME)
    if sibling.exists():
        return sibling
    fallback = Path(
        "/tmp/math-wt-frontier-synthesis-20260802/04-computation/"
        "factorial_multiplace_newton_degree_barcode_independent_audit_thm3152.py"
    )
    if fallback.exists():
        return fallback
    raise RuntimeError("independent THM-3152 engine not found")


def load_engine():
    source = engine_path()
    actual_hash = hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    if actual_hash != DEPENDENCY_SHA256:
        raise RuntimeError((source, actual_hash, DEPENDENCY_SHA256))
    spec = importlib.util.spec_from_file_location("thm3152_independent", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def residual(a, d):
    return (
        not a.prime(d)
        and a.factor_count(d - 1) != 1
        and not a.prime(d - 2)
        and not a.prime(d - 3)
    )


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def require(condition, data):
    if not condition:
        raise AssertionError(data)


def scan_chunk(chunk):
    start, end = chunk
    a = load_engine()
    ds = [d for d in range(start, end + 1) if residual(a, d)]
    traces = []
    for d in ds:
        p, q = a.moments(d)
        r = a.first_remainder(p, q, d)
        require(
            (p.degree(), q.degree(), r.degree()) == (d - 2, d - 1, d - 3),
            (d, p.degree(), q.degree(), r.degree()),
        )
        possible = set(range(1, r.degree() + 1))
        used = []
        for prime in BANK:
            narrowed = possible & a.allowed((p, q, r), prime)[0]
            if narrowed != possible:
                used.append((prime, tuple(sorted(narrowed))))
            possible = narrowed
            if not possible:
                break
        require(not possible, (d, tuple(sorted(possible))))
        traces.append((d, tuple(used)))
    return {
        "chunk": chunk,
        "ds": tuple(ds),
        "trace_digest": digest(traces),
    }


def boundary():
    a = load_engine()
    d = 1384
    p, q = a.moments(d)
    r = a.first_remainder(p, q, d)
    possible = set(range(1, r.degree() + 1))
    trace = []
    for prime in BANK:
        narrowed = possible & a.allowed((p, q, r), prime)[0]
        if narrowed != possible:
            trace.append((prime, len(narrowed), tuple(sorted(narrowed))))
        possible = narrowed
    require(possible == {1}, possible)
    local173, blocks173 = a.allowed((p, q, r), 173)
    require(not (possible & local173), (possible, local173))
    require(blocks173 == [((1, 173), 1211, 173)], blocks173)
    return trace, tuple(sorted(local173 - {0})), tuple(blocks173)


def main():
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        chunks = list(pool.map(scan_chunk, WORK_CHUNKS))
    ds = tuple(d for result in chunks for d in result["ds"])
    require(len(ds) == 511, len(ds))
    report_counts = tuple(
        sum(start <= d <= end for d in ds) for start, end in REPORT_CHUNKS
    )
    require(report_counts == (45, 100, 101, 125, 140), report_counts)
    trace, positive173, blocks173 = boundary()
    print("THM-3161 INDEPENDENT FOUR-EXIT TAIL AUDIT")
    print("independent_engine=%s sha256=%s" % (DEPENDENCY_NAME, DEPENDENCY_SHA256))
    print("universe=1001<=d<=2000 inclusive")
    print("bank=%s" % (BANK,))
    print("current four-exit residuals=%d" % len(ds))
    print("report-chunk counts=%s" % (report_counts,))
    print("closed=511 survivors=0")
    print("work-chunk independent trace digests=%s" % (tuple(result["trace_digest"] for result in chunks),))
    print("d=1384 fixed-bank trace=%s" % (trace,))
    print("d=1384 fixed-bank final=(1,)")
    print("d=1384 p=173 blocks=%s" % (blocks173,))
    print("d=1384 p=173 positive degrees=%s" % (positive173,))
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
