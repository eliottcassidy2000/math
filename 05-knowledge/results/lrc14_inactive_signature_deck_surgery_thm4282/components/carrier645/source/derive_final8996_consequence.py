#!/usr/bin/env python3
"""Independent set consequence for the detached final8996 carrier audit."""

from __future__ import annotations

import csv
import hashlib
import pathlib
import sys

OFFSET = 14695981039346656037
PRIME = 1099511628211
MASK64 = (1 << 64) - 1


def fnv_pairs(rows: list[tuple[int, int]]) -> int:
    state = OFFSET
    for row in rows:
        for value in row:
            for byte in value.to_bytes(8, "little"):
                state ^= byte
                state = (state * PRIME) & MASK64
    return state


def read_pairs(path: pathlib.Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text().splitlines():
        q_text, r_text = line.split(",")
        row = (int(q_text), int(r_text))
        if rows and not rows[-1] < row:
            raise RuntimeError("combined residual is not strictly ordered")
        rows.append(row)
    return rows


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: SCRIPT COMBINED_RESIDUAL ROWS_TSV OUTPUT_CSV")
    combined_path, rows_path, output_path = map(pathlib.Path, sys.argv[1:])
    combined = read_pairs(combined_path)
    if len(combined) != 23637 or fnv_pairs(combined) != 0xE8B363D2B3D9BA6A:
        raise RuntimeError("combined residual identity changed")
    expected_target = [row for row in combined if 645 <= row[1] <= 663]
    if len(expected_target) != 90 or fnv_pairs(expected_target) != 0x942995BEE7469430:
        raise RuntimeError("derived endpoint slice changed")

    with rows_path.open(newline="") as source:
        audited_rows = list(csv.DictReader(source, delimiter="\t"))
    observed_target = [(int(row["q"]), int(row["r"])) for row in audited_rows]
    if observed_target != expected_target:
        raise RuntimeError("audit row ledger differs from derived endpoint slice")
    if any(int(row["stage5_cut45_fail"]) != 0 for row in audited_rows):
        raise RuntimeError("final carrier has a body-cover failure")

    removed = set(observed_target)
    final = [row for row in combined if row not in removed]
    if len(final) != 23547 or max(r for _, r in final) != 644:
        raise RuntimeError("final residual count/top endpoint changed")
    top = [row for row in final if row[1] == 644]
    if len(top) != 7:
        raise RuntimeError("final residual top fibre changed")
    output_bytes = "".join(f"{q},{r}\n" for q, r in final).encode()
    output_path.write_bytes(output_bytes)
    print("THM4282_FINAL8996_DETACHED_SET_CONSEQUENCE_V1")
    print(f"COMBINED {len(combined)} FNV {fnv_pairs(combined):016x}")
    print(f"REMOVED {len(expected_target)} FNV {fnv_pairs(expected_target):016x}")
    print(
        f"FINAL {len(final)} FNV {fnv_pairs(final):016x} "
        f"SHA256 {hashlib.sha256(output_bytes).hexdigest()} MAX_R 644"
    )
    print("TOP644 COUNT 7 ROWS " + " ".join(f"({q},{r})" for q, r in top))
    print("ACCOUNTING 24223 - 586 - 90 = 23547")
    print("VERDICT PASS FINAL8996_POST_COMBINED_CONSEQUENCE")


if __name__ == "__main__":
    main()
