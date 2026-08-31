#!/usr/bin/env python3
"""Maintained typed consumer for post-9017 raw descent rows."""

from __future__ import annotations

import csv
import pathlib
import re
import sys

MASK64 = (1 << 64) - 1
OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3


def fnv(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        for byte in word.to_bytes(8, "little"):
            state ^= byte
            state = (state * PRIME) & MASK64
    return state


def main() -> None:
    if len(sys.argv) != 3:
        raise RuntimeError("usage: consumer DESCENT_OUT FAILURE_CSV")
    text = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
    pair_pattern = re.compile(
        r"^PAIR (\d+),(\d+) .* FAILURES (\d+) FAILURE_FNV ([0-9a-f]+)$",
        re.MULTILINE,
    )
    pairs = [
        (int(match.group(1)), int(match.group(2)), int(match.group(3)),
         int(match.group(4), 16))
        for match in pair_pattern.finditer(text)
    ]
    if len(pairs) != 48 or len({(q, r) for q, r, _, _ in pairs}) != 48:
        raise RuntimeError("raw pair row count/distinctness changed")
    success = [(q, r) for q, r, failures, _ in pairs if failures == 0]
    failed = [(q, r, failures, digest)
              for q, r, failures, digest in pairs if failures != 0]
    if len(success) != 47 or failed != [(100, 628, 4, 0x2E7F9DDCCC403B9D)]:
        raise RuntimeError("raw success/failure partition changed")
    if "COMPLETED_ROWS 48 TOTAL_FAILURES 4 STOPPED 1 LEDGER_FNV be2cbab7fb58d91e" not in text:
        raise RuntimeError("descent summary changed")
    for endpoint in range(636, 628, -1):
        endpoint_rows = [row for row in pairs if row[1] == endpoint]
        if not endpoint_rows or any(row[2] != 0 for row in endpoint_rows):
            raise RuntimeError(f"endpoint {endpoint} not closed")

    with pathlib.Path(sys.argv[2]).open(newline="", encoding="utf-8") as handle:
        body_rows = list(csv.DictReader(handle))
    bodies = [(int(row["q"]), int(row["r"]), int(row["body_hex"], 16))
              for row in body_rows]
    expected_bodies = [
        (100, 628, 0x05346408),
        (100, 628, 0x15306408),
        (100, 628, 0x17581400),
        (100, 628, 0x27D01008),
    ]
    if bodies != expected_bodies or any(body.bit_count() != 9 for _, _, body in bodies):
        raise RuntimeError("raw failure bodies changed")
    if fnv([body for _, _, body in bodies]) != 0x2E7F9DDCCC403B9D:
        raise RuntimeError("failure body FNV changed")

    success_words = [word for q, r in success for word in (q, r)]
    failure_words = [word for q, r, count, digest in failed
                     for word in (q, r, count, digest)]
    print("R629_POST9017_DESCENT_CONSUMER_V1")
    print(f"RAW_PAIR_ROWS {len(pairs)} SUCCESS_ROWS {len(success)} "
          f"FAILURE_ROWS {len(failed)}")
    print(f"SUCCESS_PAIR_FNV {fnv(success_words):016x} "
          f"FAILURE_ROW_FNV {fnv(failure_words):016x}")
    print("CLOSED_ENDPOINTS 636..629")
    print("FIRST_FAILED_PAIR 100,628 FAILURES 4 FNV 2e7f9ddccc403b9d")
    print("FAILURE_BODIES " + " ".join(f"{body:08x}" for _, _, body in bodies))
    print("COMPLETED_ROWS 48 LEDGER_FNV be2cbab7fb58d91e")
    print("VERDICT PASS EXACT_RAW_SUCCESS_FAILURE_ROWS")


if __name__ == "__main__":
    main()
