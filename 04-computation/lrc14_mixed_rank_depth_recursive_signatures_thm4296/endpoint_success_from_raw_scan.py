#!/usr/bin/env python3
"""Derive a pair-success ledger from a completed raw carrier descent scan.

This parser deliberately uses only ``PAIR ... FAILURES n`` rows.  Endpoint
cutoffs, expected counts, and previously generated success ledgers play no
role.  It rejects partial scans and checks the scan's own completion totals.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PAIR_RE = re.compile(r"^PAIR\s+(\d+),(\d+)\s+.*\bFAILURES\s+(\d+)\b")
SUMMARY_RE = re.compile(
    r"^COMPLETED_ROWS\s+(\d+)\s+TOTAL_FAILURES\s+(\d+)\s+STOPPED\s+([01])\b"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse_scan(path: Path) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    rows: list[tuple[int, int, int]] = []
    summary: tuple[int, int, int] | None = None
    saw_verdict = False
    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            match = PAIR_RE.match(line)
            if match:
                q, r, failures = map(int, match.groups())
                require(0 < q < r, f"invalid pair in scan: {q},{r}")
                rows.append((q, r, failures))
                continue
            match = SUMMARY_RE.match(line)
            if match:
                require(summary is None, "duplicate completion summary")
                summary = tuple(map(int, match.groups()))
                continue
            if line.startswith("VERDICT PASS "):
                saw_verdict = True

    require(rows, "scan has no PAIR rows")
    require(summary is not None, "scan is partial: missing COMPLETED_ROWS summary")
    require(saw_verdict, "scan is partial or failed: missing PASS verdict")
    completed, total_failures, _stopped = summary
    require(completed == len(rows), "COMPLETED_ROWS disagrees with PAIR rows")
    require(total_failures == sum(item[2] for item in rows),
            "TOTAL_FAILURES disagrees with PAIR rows")

    pairs = [(q, r) for q, r, _ in rows]
    require(len(set(pairs)) == len(pairs), "duplicate PAIR row")
    expected_order = sorted(pairs, key=lambda pair: (-pair[1], pair[0]))
    require(pairs == expected_order, "PAIR rows are not descending-endpoint/ascending-q")

    success = sorted((q, r) for q, r, failures in rows if failures == 0)
    failure = sorted((q, r) for q, r, failures in rows if failures != 0)
    return success, failure


def write_pairs(path: Path, rows: list[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for q, r in rows:
            handle.write(f"{q},{r}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scan", type=Path)
    parser.add_argument("success_csv", type=Path)
    parser.add_argument("--failure-csv", type=Path)
    args = parser.parse_args()

    success, failure = parse_scan(args.scan)
    write_pairs(args.success_csv, success)
    if args.failure_csv is not None:
        write_pairs(args.failure_csv, failure)
    print(f"RAW_SCAN rows={len(success) + len(failure)} successes={len(success)} "
          f"failed_pairs={len(failure)}")
    print("FAILED_PAIR_SET " +
          (" ".join(f"({q},{r})" for q, r in failure) if failure else "EMPTY"))
    print(f"WROTE {args.success_csv}")


if __name__ == "__main__":
    main()
