#!/usr/bin/env python3
"""Independent raw-cell wall audit for endpoint-587 carrier misses."""

from __future__ import annotations

import argparse
import collections
import csv
import math
from pathlib import Path

POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
ENDPOINT = 587
FNV_BASIS = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
WORD_MASK = (1 << 64) - 1


class Fnv:
    def __init__(self) -> None:
        self.value = FNV_BASIS

    def add(self, word: int) -> None:
        word &= WORD_MASK
        for shift in range(0, 64, 8):
            self.value ^= (word >> shift) & 0xFF
            self.value = (self.value * FNV_PRIME) & WORD_MASK

    def add_i128(self, value: int) -> None:
        bits = value % (1 << 128)
        self.add(bits)
        self.add(bits >> 64)


def safe(speed: int, grid: int, left: int, right: int) -> bool:
    twice_phase = speed * (left + right) % (2 * grid)
    return 7 * twice_phase >= grid and 7 * twice_phase <= 13 * grid


def geometry(q: int) -> tuple[int, int, list[tuple[int, int]], int]:
    grid = 1
    for speed in (*POOL, q, ENDPOINT):
        grid = math.lcm(grid, 14 * speed)
    walls = {0, grid}
    for speed in (*POOL, q, ENDPOINT):
        unit, remainder = divmod(grid, 14 * speed)
        if remainder:
            raise AssertionError("nonintegral wall unit")
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    cells: list[tuple[int, int]] = []
    pair_ticks = 0
    for left, right in zip(ordered, ordered[1:]):
        if not safe(q, grid, left, right) or not safe(ENDPOINT, grid, left, right):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe(speed, grid, left, right):
                failure |= 1 << bit
        width = right - left
        pair_ticks += width
        cells.append((failure, width))
    return grid, len(ordered) - 1, cells, pair_ticks


def read_failures(path: Path) -> tuple[dict[int, list[int]], int]:
    rows: dict[int, list[int]] = collections.defaultdict(list)
    global_hash = Fnv()
    previous: tuple[int, int] | None = None
    with path.open("r", encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["q", "r", "body_hex"]:
            raise AssertionError("failure header changed")
        for record in reader:
            q = int(record["q"])
            r = int(record["r"])
            body = int(record["body_hex"], 16)
            if r != ENDPOINT or body.bit_count() != 9:
                raise AssertionError("failure row escaped endpoint/rank")
            if previous is not None and previous >= (q, body):
                raise AssertionError("failure row order changed")
            previous = (q, body)
            rows[q].append(body)
            global_hash.add(q)
            global_hash.add(r)
            global_hash.add(body)
    if not rows:
        raise AssertionError("empty failure ledger")
    return dict(rows), global_hash.value


def aggregate(cells: list[tuple[int, int]]) -> dict[int, int]:
    answer: dict[int, int] = collections.defaultdict(int)
    for failure, width in cells:
        answer[failure] += width
    return dict(answer)


def class_hash(classes: dict[int, int]) -> int:
    ledger = Fnv()
    for failure in sorted(classes):
        ledger.add(failure)
        ledger.add(classes[failure])
    return ledger.value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--detail", type=Path, required=True)
    args = parser.parse_args()

    failures, global_failure_hash = read_failures(args.failures)
    total_bodies = sum(map(len, failures.values()))
    total_low_positive = 0
    total_full_positive = 0
    total_equality = 0
    summary_hash = Fnv()
    output_lines: list[str] = []

    with args.detail.open("w", encoding="ascii", newline="") as detail:
        detail.write(
            "q,r,ordinal,body_hex,truncated_mass,full_mass,"
            "truncated_scaled_ticks,full_scaled_ticks\n"
        )
        for q in sorted(failures):
            bodies = failures[q]
            grid, all_cell_count, cells, pair_ticks = geometry(q)
            all_classes = aggregate(cells)
            low_cells = [cell for cell in cells if cell[0].bit_count() <= 9]
            low_classes = aggregate(low_cells)
            row_body_hash = Fnv()
            row_detail_hash = Fnv()
            values: list[tuple[int, int, int, int, int]] = []
            for ordinal, body in enumerate(bodies):
                row_body_hash.add(body)
                low_mass = sum(width for failure, width in low_cells
                               if failure & body == 0)
                full_mass = sum(width for failure, width in cells
                                if failure & body == 0)
                if low_mass > full_mass:
                    raise AssertionError("truncation exceeded full mass")
                low_ticks = 63 * low_mass - 4 * grid
                full_ticks = 63 * full_mass - 4 * grid
                detail.write(
                    f"{q},{ENDPOINT},{ordinal},{body:08x},{low_mass},"
                    f"{full_mass},{low_ticks},{full_ticks}\n"
                )
                row_detail_hash.add(q)
                row_detail_hash.add(ENDPOINT)
                row_detail_hash.add(ordinal)
                row_detail_hash.add(body)
                row_detail_hash.add(low_mass)
                row_detail_hash.add(full_mass)
                row_detail_hash.add_i128(low_ticks)
                row_detail_hash.add_i128(full_ticks)
                values.append((body, low_ticks, full_ticks, low_mass, full_mass))

            low_min_body, low_min, _, _, _ = min(
                values, key=lambda item: (item[1], item[0]))
            low_max_body, low_max, _, _, _ = max(
                values, key=lambda item: (item[1], -item[0]))
            full_min_body, _, full_min, _, _ = min(
                values, key=lambda item: (item[2], item[0]))
            full_max_body, _, full_max, _, _ = max(
                values, key=lambda item: (item[2], -item[0]))
            low_positive = sum(low_ticks > 0 for _, low_ticks, _, _, _ in values)
            full_positive = sum(full_ticks > 0 for _, _, full_ticks, _, _ in values)
            equality = sum(low_mass == full_mass
                           for _, _, _, low_mass, full_mass in values)
            total_low_positive += low_positive
            total_full_positive += full_positive
            total_equality += equality

            all_hash = class_hash(all_classes)
            low_hash = class_hash(low_classes)
            summary_hash.add(q)
            summary_hash.add(len(bodies))
            summary_hash.add(grid)
            summary_hash.add(all_cell_count)
            summary_hash.add(len(cells))
            summary_hash.add(pair_ticks)
            summary_hash.add(len(all_classes))
            summary_hash.add(len(low_classes))
            summary_hash.add(all_hash)
            summary_hash.add(low_hash)
            summary_hash.add_i128(low_min)
            summary_hash.add(low_min_body)
            summary_hash.add_i128(full_min)
            summary_hash.add(full_min_body)
            summary_hash.add(low_positive)
            summary_hash.add(full_positive)
            summary_hash.add(equality)
            summary_hash.add(row_detail_hash.value)

            output_lines.append(
                f"ROW {q},{ENDPOINT} BODIES {len(bodies)} BODY_FNV "
                f"{row_body_hash.value:x} GRID {grid} CELLS {all_cell_count} "
                f"PAIR_SAFE_CELLS {len(cells)} PAIR_TICKS {pair_ticks} "
                f"ALL_CLASSES {len(all_classes)} LOW_CLASSES {len(low_classes)} "
                f"ALL_CLASS_FNV {all_hash:x} LOW_CLASS_FNV {low_hash:x} "
                f"LOW_POSITIVE {low_positive} FULL_POSITIVE {full_positive} "
                f"EQUALITY {equality} LOW_MIN_TICKS {low_min} "
                f"LOW_MIN_BODY {low_min_body:08x} FULL_MIN_TICKS {full_min} "
                f"FULL_MIN_BODY {full_min_body:08x} LOW_MAX_TICKS {low_max} "
                f"LOW_MAX_BODY {low_max_body:08x} FULL_MAX_TICKS {full_max} "
                f"FULL_MAX_BODY {full_max_body:08x} DETAIL_FNV "
                f"{row_detail_hash.value:x}"
            )

    print("LRC14_ENDPOINT587_CLEANROOM_RAWCELL_WALL_AUDIT_V1")
    print(
        f"FAILURE_ROWS {len(failures)} BODIES {total_bodies} "
        f"FAILURE_FNV {global_failure_hash:x}"
    )
    for line in output_lines:
        print(line)
    print(
        f"SUMMARY BODIES {total_bodies} LOW_POSITIVE {total_low_positive} "
        f"FULL_POSITIVE {total_full_positive} EQUALITY {total_equality} "
        f"SUMMARY_FNV {summary_hash.value:x}"
    )
    print("REPRESENTATION RAW_PAIR_SAFE_WALL_CELLS_PYTHON_INTEGER_ARITHMETIC")
    print("INEQUALITY TRUNCATED_LOW_RANK_CELL_MASS_LE_LITERAL_SAFE_MASS")
    print(
        "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT587_FAILURE_BODIES_"
        "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14"
    )
    print("VERDICT " + ("PASS" if total_low_positive == total_bodies else "HOSTILE_FAIL"))


if __name__ == "__main__":
    main()
