#!/usr/bin/env python3
"""Exact direct-activity 2-for-1 audit for the frozen q50 95-cover."""

from __future__ import annotations

import csv
import itertools
import math
import sys
from collections import defaultdict
from pathlib import Path

POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def safe_midpoint(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return 7 * residue >= grid and 7 * residue <= 13 * grid


def build_geometry() -> tuple[int, tuple[tuple[int, int], ...]]:
    grid = 1
    for speed in (*POOL, 50, 589):
        grid = math.lcm(grid, 14 * speed)
    walls = {0, grid}
    for speed in (*POOL, 50, 589):
        unit = grid // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    classes: defaultdict[int, int] = defaultdict(int)
    for left, right in zip(ordered, ordered[1:]):
        if not (safe_midpoint(50, grid, left, right) and
                safe_midpoint(589, grid, left, right)):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe_midpoint(speed, grid, left, right):
                failure |= 1 << bit
        if failure.bit_count() <= 9:
            classes[failure] += right - left
    require(grid == 2827379709554400 and len(classes) == 2383,
            "q50 geometry changed")
    return grid, tuple(sorted(classes.items()))


def read_hex_column(path: Path, column: str, predicate=lambda row: True) -> list[int]:
    with path.open(newline="", encoding="ascii") as handle:
        return [int(row[column], 16) for row in csv.DictReader(handle) if predicate(row)]


def main() -> None:
    packet = Path(sys.argv[1])
    bodies = read_hex_column(
        packet / "input_endpoint589_failures.csv", "body_hex", lambda row: row["q"] == "50"
    )
    residual = [body for body in bodies if body & 0x0220932C]
    cover = read_hex_column(packet / "input_cover95.csv", "mask_hex")
    incidence = [[mask & body == 0 for body in residual] for mask in cover]
    multiplicity = [sum(column) for column in zip(*incidence)]
    require(len(residual) == 19134 and len(cover) == 95 and min(multiplicity) >= 1,
            "cover input changed")
    grid, classes = build_geometry()

    def active(mask: int) -> bool:
        mass = sum(width for failure, width in classes if failure & ~mask == 0)
        return 63 * mass - 4 * grid >= 0

    complement_histogram: dict[int, int] = {}
    candidates_tested = 0
    active_candidates = 0
    exchanges: list[tuple[int, int, int, int]] = []
    universe = (1 << 30) - 1
    for left, right in itertools.combinations(range(len(cover)), 2):
        required_union = 0
        required_count = 0
        for index, body in enumerate(residual):
            removed = int(incidence[left][index]) + int(incidence[right][index])
            if removed and multiplicity[index] == removed:
                required_union |= body
                required_count += 1
        available = universe ^ required_union
        width = available.bit_count()
        complement_histogram[width] = complement_histogram.get(width, 0) + 1
        if width < 8:
            continue
        labels = [bit for bit in range(30) if available >> bit & 1]
        found = 0
        for rank in (8, 9):
            for chosen in itertools.combinations(labels, rank):
                mask = sum(1 << bit for bit in chosen)
                candidates_tested += 1
                if active(mask):
                    active_candidates += 1
                    found = mask
                    break
            if found:
                break
        if found:
            exchanges.append((left, right, required_count, found))

    print("Q50_COVER95_TWO_FOR_ONE_DIRECT_AUDIT_V1")
    print(f"RESIDUAL {len(residual)} COVER {len(cover)} PAIRS {len(cover)*(len(cover)-1)//2}")
    print("AVAILABLE_WIDTH_HIST " + ";".join(f"{key}:{value}" for key, value in sorted(complement_histogram.items())))
    print(f"ENUMERATED {candidates_tested} ACTIVE_FOUND {active_candidates} EXCHANGES {len(exchanges)}")
    for left, right, required, mask in exchanges:
        print(f"EXCHANGE {left+1} {right+1} REQUIRED {required} MASK {mask:08x}")
    print("VERDICT " + ("IMPROVABLE" if exchanges else "TWO_FOR_ONE_LOCAL_MINIMUM"))


if __name__ == "__main__":
    main()
