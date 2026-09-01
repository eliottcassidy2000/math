#!/usr/bin/env python3
"""Compare clean-room ledgers with a separately implemented primary run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

PRIMARY_PAIR = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) BOUND_EXACT 1 MIN_TICKS (-?\d+) "
    r"MIN_BODY ([0-9a-f]{8}) NODES \d+ PRUNES \d+$"
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("independent_screen", type=Path)
    parser.add_argument("primary_screen", type=Path)
    parser.add_argument("independent_exact", type=Path)
    parser.add_argument("primary_exact_out", type=Path)
    args = parser.parse_args()

    with args.independent_screen.open(newline="", encoding="utf-8") as handle:
        independent = {
            (int(row["q"]), int(row["r"])): row for row in csv.DictReader(handle)
        }
    with args.primary_screen.open(newline="", encoding="utf-8") as handle:
        primary = {
            (int(row["q"]), int(row["r"])): row for row in csv.DictReader(handle)
        }
    assert len(independent) == len(primary) == 181_194
    assert independent.keys() == primary.keys()
    for pair, left in independent.items():
        right = primary[pair]
        assert left["grid"] == right["grid"]
        assert left["total_mass"] == right["rank2_total"]
        assert left["coarse_mass"] == right["degree_bound_mass"]
        assert left["coarse_ticks"] == right["degree_bound_ticks"]
        assert left["coarse_positive"] == right["positive"]
        assert left["top9_mask_hex"] == right["top9_hex"]

    with args.independent_exact.open(newline="", encoding="utf-8") as handle:
        exact = {
            (int(row["q"]), int(row["r"])): row for row in csv.DictReader(handle)
        }
    primary_exact: dict[tuple[int, int], tuple[str, str, str]] = {}
    for line in args.primary_exact_out.read_text(encoding="utf-8").splitlines():
        match = PRIMARY_PAIR.fullmatch(line.strip())
        if match:
            primary_exact[(int(match.group(1)), int(match.group(2)))] = (
                match.group(3), match.group(4), match.group(5)
            )
    assert len(exact) == len(primary_exact) == 107
    assert exact.keys() == primary_exact.keys()
    for pair, row in exact.items():
        grid, ticks, body = primary_exact[pair]
        assert row["grid"] == grid
        assert row["minimum_ticks"] == ticks
        assert row["minimum_body_hex"] == body

    print("LRC14_RANK2_PRIMARY_CLEANROOM_CROSSCHECK_V1")
    print("SCREEN_PAIRS 181194 FIELDWISE_MATCHES 181194")
    print("EXACT_PAIRS 107 MINIMUM_MATCHES 107")
    print(f"PRIMARY_SCREEN_SHA256 {sha(args.primary_screen)}")
    print(f"PRIMARY_EXACT_SHA256 {sha(args.primary_exact_out)}")
    print("INDEPENDENT_GRAPH_PATH EVENT_STATE_MERGE")
    print("PRIMARY_GRAPH_PATH SORTED_WALL_MIDPOINT_CLASSIFICATION")
    print("VERDICT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
