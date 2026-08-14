#!/usr/bin/env python3
"""Audit the combined 264--698 ledger via independently audited segments."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path


TARGET = F(186_636_088_362, 58_865_718_786_875)
EXPECTED_LEDGER = "932057abf10f674e4bb31f334c1ea94f39e4e17627c6a950bd5d727f8e595186"
EXPECTED_TASKS = "079e54a95dd55135dddcb9093bb98fba261939ce9cd9c9aeb7b38c8fe45866b7"
SEGMENTS = (
    (264, 454, "2865e79a9aed02b293329cc7ed219ca0f428a376175f9b4fbb73552b43127053", 397_502),
    (455, 678, "f12b9cc843fd3d4dd8e300a2461df00721778f8665ab096d322a31dcc53b908a", 734_566),
    (679, 698, "0915c0c3a6eaf9fa4577e4449900865992d8fa646184e4ba6600801e89a55421", 79_898),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def digest(data):
    return sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def main():
    parser = ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("segments", nargs=3, type=Path)
    args = parser.parse_args()
    raw = args.ledger.read_bytes()
    require(digest(raw) == EXPECTED_LEDGER, digest(raw))
    pieces = []
    for path, (lo, hi, wanted, count) in zip(args.segments, SEGMENTS):
        data = path.read_bytes()
        require(digest(data) == wanted, (path, digest(data)))
        require(len(data.splitlines()) == count, (path, len(data.splitlines())))
        pieces.append(data)
    require(b"".join(pieces) == raw, "combined ledger differs from audited segments")
    rows = tuple(tuple(map(int, line.split())) for line in raw.splitlines())
    expected = tuple(
        (p, q)
        for p in range(264, 699)
        for q in range(p + 1, 8 * p)
        if gcd(p, q) <= 3
    )
    require(tuple(row[:2] for row in rows) == expected, (len(rows), len(expected)))
    payload = "".join(f"{p} {q}\n" for p, q in expected).encode()
    require(sha256(payload).hexdigest() == EXPECTED_TASKS, "task digest")
    require(all(F(row[2], row[3]) > TARGET for row in rows), "threshold failure")
    weakest = min(rows, key=lambda row: F(row[2], row[3]))
    require(weakest == (698, 2559, 20682154, 1400220127, 168, 90, 12, 1), weakest)
    print("LRC14 DISCONNECTED RAW BRIDGE 264--698 INDEPENDENT ASSEMBLY AUDIT")
    print("rows", len(rows), "mass_comparisons", len(rows) * 2530, "failures", 0)
    print("segment_counts", tuple(x[3] for x in SEGMENTS))
    print("task_semantic_sha256", sha256(payload).hexdigest())
    print("ledger_sha256", digest(raw))
    print("weakest", weakest, "margin", F(weakest[2], weakest[3]) - TARGET)
    print("status=PASS; each segment has a separate all-argmin reference audit")


if __name__ == "__main__":
    main()
