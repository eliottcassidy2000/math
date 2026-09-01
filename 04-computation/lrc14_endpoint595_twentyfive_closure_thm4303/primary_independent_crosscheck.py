#!/usr/bin/env python3
"""Compare all common row fields in THM-4303's two raw implementations."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

Pair = tuple[int, int]
PATTERN = re.compile(
    r"^PAIR (\d+),(\d+) ACTIVE (\d+) ACTIVE_FNV ([0-9a-f]+) "
    r"JOINT (\d+) NONJOINT (\d+) EXPOSED (\d+) EXPOSED_FNV ([0-9a-f]+) "
    r"FAILURES (\d+) FAILURE_FNV ([0-9a-f]+)$"
)
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def read_primary(path: Path) -> dict[Pair, tuple[int | str, ...]]:
    answer: dict[Pair, tuple[int | str, ...]] = {}
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            pair = int(row["q"]), int(row["r"])
            answer[pair] = (
                int(row["active"]), f"{int(row['active_fnv'], 16):016x}",
                int(row["active_joint"]), int(row["active_nonjoint"]),
                int(row["exposed"]), f"{int(row['exposed_fnv'], 16):016x}",
                int(row["failures"]), f"{int(row['failure_fnv'], 16):016x}",
            )
    return answer


def read_independent(path: Path) -> dict[Pair, tuple[int | str, ...]]:
    answer: dict[Pair, tuple[int | str, ...]] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("PAIR "):
            continue
        match = PATTERN.fullmatch(line)
        require(match is not None, "malformed independent PAIR line")
        values = match.groups()
        pair = int(values[0]), int(values[1])
        require(pair not in answer, "duplicate independent pair")
        answer[pair] = (
            int(values[2]), f"{int(values[3], 16):016x}",
            int(values[4]), int(values[5]), int(values[6]),
            f"{int(values[7], 16):016x}", int(values[8]),
            f"{int(values[9], 16):016x}",
        )
    return answer


def pair_fnv(rows: set[Pair]) -> str:
    state = FNV_OFFSET
    for pair in sorted(rows):
        for value in pair:
            for byte in value.to_bytes(8, "little"):
                state = ((state ^ byte) * FNV_PRIME) & MASK64
    return f"{state:016x}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primary_csv", type=Path)
    parser.add_argument("independent_output", type=Path)
    args = parser.parse_args()

    primary = read_primary(args.primary_csv)
    independent = read_independent(args.independent_output)
    require(primary == independent, "primary/independent row values differ")
    require(len(primary) == 28 and pair_fnv(set(primary)) == "47981ce64825ef2a",
            "cross-checked row universe changed")
    failures = {pair: values[6] for pair, values in primary.items()
                if values[6] != 0}
    require(failures == {(96, 595): 116, (100, 595): 13, (210, 595): 16},
            "cross-checked failure split changed")

    print("THM4303_PRIMARY_INDEPENDENT_CROSSCHECK_V1")
    print("ROWS 28 FNV 47981ce64825ef2a COMMON_FIELDS 8")
    print("FIELDS active active_fnv joint nonjoint exposed exposed_fnv "
          "failures failure_fnv")
    print("FAILURE_PAIRS (96,595):116 (100,595):13 (210,595):16")
    print("SCOPE FINITE_EXACT_OUTPUT_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS EXACT_FIELD_MATCH_ALL_28_ROWS")


if __name__ == "__main__":
    main()
