#!/usr/bin/env python3
"""Exact filter certificate for raw and grid-normalized rank-two minima."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

PAIR_RE = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) COARSE_TICKS (-?\d+) "
    r"MIN_MASS (\d+) MIN_TICKS (-?\d+) MIN_BODY ([0-9a-f]{8}) "
    r"DIRECT_MASS (\d+) RATIO_CROSS_VS_50_70 (-?\d+) "
    r"NODES (\d+) PRUNES (\d+)$"
)


def selected_rows(path: Path) -> dict[tuple[int, int], dict[str, int | str]]:
    rows: dict[tuple[int, int], dict[str, int | str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = PAIR_RE.fullmatch(line.strip())
        if not match:
            continue
        pair = (int(match.group(1)), int(match.group(2)))
        rows[pair] = {
            "grid": int(match.group(3)),
            "coarse_ticks": int(match.group(4)),
            "minimum_mass": int(match.group(5)),
            "minimum_ticks": int(match.group(6)),
            "body": match.group(7),
            "direct_mass": int(match.group(8)),
            "ratio_cross": int(match.group(9)),
        }
    assert len(rows) == 4
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("broad_screen", type=Path)
    parser.add_argument("hostile_exact", type=Path)
    parser.add_argument("selected_exact", type=Path)
    args = parser.parse_args()

    with args.broad_screen.open(newline="", encoding="utf-8") as handle:
        screen = list(csv.DictReader(handle))
    with args.hostile_exact.open(newline="", encoding="utf-8") as handle:
        hostile = list(csv.DictReader(handle))
    selected = selected_rows(args.selected_exact)

    assert len(screen) == 181_194
    assert len(hostile) == 107
    bad_pairs = {
        (int(row["q"]), int(row["r"]))
        for row in screen
        if row["coarse_positive"] == "0"
    }
    assert len(bad_pairs) == 107
    assert bad_pairs == {(int(row["q"]), int(row["r"])) for row in hostile}
    assert all(int(row["minimum_ticks"]) > 0 for row in hostile)

    candidate = min(
        hostile,
        key=lambda row: int(row["minimum_ticks"]) / int(row["grid"]),
    )
    # Recheck the floating-free order and the exact advertised witness.
    candidate_pair = (int(candidate["q"]), int(candidate["r"]))
    candidate_ticks = int(candidate["minimum_ticks"])
    candidate_grid = int(candidate["grid"])
    candidate_body = candidate["minimum_body_hex"]
    assert (candidate_pair, candidate_ticks, candidate_grid, candidate_body) == (
        (50, 70), 245_428_469_244, 91_205_797_082_400, "031c7400"
    )
    for row in hostile:
        assert (
            int(row["minimum_ticks"]) * candidate_grid
            >= candidate_ticks * int(row["grid"])
        )

    positive = [row for row in screen if row["coarse_positive"] == "1"]
    assert len(positive) == 181_087
    raw_filter = [
        (int(row["q"]), int(row["r"]))
        for row in positive
        if int(row["coarse_ticks"]) <= candidate_ticks
    ]
    ratio_filter = [
        (int(row["q"]), int(row["r"]))
        for row in positive
        if int(row["coarse_ticks"]) * candidate_grid
        <= candidate_ticks * int(row["grid"])
    ]
    assert raw_filter == [(100, 110)], raw_filter
    assert ratio_filter == [(50, 212), (50, 274), (100, 110)], ratio_filter

    for pair in raw_filter:
        assert int(selected[pair]["minimum_ticks"]) > candidate_ticks
    for pair in ratio_filter:
        row = selected[pair]
        assert (
            int(row["minimum_ticks"]) * candidate_grid
            > candidate_ticks * int(row["grid"])
        )

    divisor = math.gcd(candidate_ticks, candidate_grid)
    print("LRC14_RANK2_GLOBAL_FILTER_CERTIFICATE_V1")
    print(
        "COARSE_POSITIVE 181087 COARSE_NONPOSITIVE 107 "
        "HOSTILE_EXACT_POSITIVE 107"
    )
    print("RAW_COARSE_AT_OR_BELOW_CANDIDATE 1 PAIRS 100,110")
    print(
        "NORMALIZED_COARSE_AT_OR_BELOW_CANDIDATE 3 PAIRS "
        "50,212;50,274;100,110"
    )
    print(
        "GLOBAL_RANK2_MINIMUM_TICKS 245428469244 PAIR 50,70 "
        "BODY 031c7400 GRID 91205797082400"
    )
    print(
        "GLOBAL_NORMALIZED_RANK2_MINIMUM_RATIO "
        f"{candidate_ticks // divisor}/{candidate_grid // divisor} "
        "PAIR 50,70 BODY 031c7400"
    )
    print(
        "LOGIC ALL_COARSE_NONPOSITIVE_ROWS_EXACT_PLUS_ALL_COARSE_POSITIVE_"
        "ROWS_NOT_EXCLUDED_BY_RAW_OR_RATIO_LOWER_BOUND_EXACT"
    )
    print("SCOPE FINITE_EXACT_THM4231_REMAINDER181194_RANK2_FUNCTIONAL_ONLY")
    print("VERDICT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
