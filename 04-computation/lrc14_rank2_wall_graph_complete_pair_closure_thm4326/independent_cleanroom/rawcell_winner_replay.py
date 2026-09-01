#!/usr/bin/env python3
"""Independent midpoint/wall replay of every exact winning body.

Unlike the C++ event-state sweep, this verifier reconstructs a fresh sorted set
of all pool/q/r walls and classifies each open cell from its midpoint residue.
It sums only cells whose failure set has rank at most two and misses the body.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
BASE_GRID = 18_241_159_416_480
PAIR_RE = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) COARSE_TICKS (-?\d+) "
    r"MIN_MASS (\d+) MIN_TICKS (-?\d+) MIN_BODY ([0-9a-f]{8}) "
    r"DIRECT_MASS (\d+) RATIO_CROSS_VS_50_70 (-?\d+) "
    r"NODES (\d+) PRUNES (\d+)$"
)


def safe_midpoint(speed: int, grid: int, left: int, right: int) -> bool:
    residue = (speed * (left + right)) % (2 * grid)
    return 7 * residue >= grid and 7 * residue <= 13 * grid


def all_walls(q: int, r: int, grid: int) -> list[int]:
    walls = {0, grid}
    for speed in (*POOL, q, r):
        divisor = 14 * speed
        assert grid % divisor == 0
        unit = grid // divisor
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    return sorted(walls)


def replay(q: int, r: int, body: int) -> tuple[int, int, int]:
    grid = math.lcm(BASE_GRID, 14 * q, 14 * r)
    walls = all_walls(q, r, grid)
    mass = 0
    accepted_cells = 0
    for left, right in zip(walls, walls[1:]):
        if not safe_midpoint(q, grid, left, right):
            continue
        if not safe_midpoint(r, grid, left, right):
            continue
        failures = 0
        rank = 0
        for bit, speed in enumerate(POOL):
            if not safe_midpoint(speed, grid, left, right):
                failures |= 1 << bit
                rank += 1
                if rank > 2:
                    break
        if rank <= 2 and failures & body == 0:
            mass += right - left
            accepted_cells += 1
    return grid, mass, accepted_cells


def load_expected(exact_csv: Path, selected_out: Path) -> dict[tuple[int, int], tuple[int, int, int]]:
    expected: dict[tuple[int, int], tuple[int, int, int]] = {}
    with exact_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pair = (int(row["q"]), int(row["r"]))
            expected[pair] = (
                int(row["minimum_body_hex"], 16),
                int(row["minimum_mass"]),
                int(row["minimum_ticks"]),
            )
    for line in selected_out.read_text(encoding="utf-8").splitlines():
        match = PAIR_RE.fullmatch(line.strip())
        if not match:
            continue
        pair = (int(match.group(1)), int(match.group(2)))
        value = (int(match.group(7), 16), int(match.group(5)), int(match.group(6)))
        if pair in expected:
            assert expected[pair] == value
        else:
            expected[pair] = value
    assert len(expected) == 110, len(expected)
    return expected


def fnv_add(value: int, line: str) -> int:
    for byte in (line + "\n").encode("ascii"):
        value ^= byte
        value = (value * 1_099_511_628_211) & ((1 << 64) - 1)
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("exact_csv", type=Path)
    parser.add_argument("selected_out", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    expected = load_expected(args.exact_csv, args.selected_out)
    rows: list[str] = []
    for (q, r), (body, expected_mass, expected_ticks) in sorted(expected.items()):
        assert body.bit_count() == 9
        grid, mass, cells = replay(q, r, body)
        ticks = 63 * mass - 4 * grid
        matched = mass == expected_mass and ticks == expected_ticks
        if not matched:
            raise AssertionError((q, r, mass, ticks, expected_mass, expected_ticks))
        rows.append(
            f"{q},{r},{grid},{body:08x},{cells},{mass},{ticks},"
            f"{expected_mass},{expected_ticks},1"
        )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "q,r,grid,body_hex,accepted_cells,replay_mass,replay_ticks,"
        "expected_mass,expected_ticks,match"
    )
    args.output_csv.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="ascii")
    ledger_fnv = 14_695_981_039_346_656_037
    for row in rows:
        ledger_fnv = fnv_add(ledger_fnv, row)
    print("LRC14_RANK2_RAWCELL_WINNER_REPLAY_V1")
    print(f"PAIRS {len(rows)} MATCHES {len(rows)} DATA_FNV {ledger_fnv:016x}")
    print("PATH SORTED_WALLS_PLUS_DIRECT_MIDPOINT_RESIDUES")
    print("SCOPE FINITE_EXACT_WINNING_BODIES_NO_OPTIMALITY_CLAIM_FROM_THIS_REPLAY")
    print("VERDICT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
