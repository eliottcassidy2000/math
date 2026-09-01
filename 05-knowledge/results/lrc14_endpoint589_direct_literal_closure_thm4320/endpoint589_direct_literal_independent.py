#!/usr/bin/env python3
"""Clean-room exact wall replay for endpoint-589 carrier-failure bodies."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
ENDPOINT = 589
OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(state: int, word: int) -> int:
    word &= (1 << 128) - 1
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def fnv_i128(state: int, word: int) -> int:
    bits = word % (1 << 128)
    state = fnv_add(state, bits & MASK64)
    return fnv_add(state, bits >> 64)


def safe_midpoint(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return 7 * residue >= grid and 7 * residue <= 13 * grid


def build_geometry(q: int) -> tuple[int, int, int, dict[int, int], int]:
    grid = 1
    for speed in (*POOL, q, ENDPOINT):
        grid = math.lcm(grid, 14 * speed)
    walls = {0, grid}
    for speed in (*POOL, q, ENDPOINT):
        unit = grid // (14 * speed)
        require(unit * 14 * speed == grid, "nonintegral wall unit")
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    classes: defaultdict[int, int] = defaultdict(int)
    pair_ticks = 0
    for left, right in zip(ordered, ordered[1:]):
        if not (safe_midpoint(q, grid, left, right) and
                safe_midpoint(ENDPOINT, grid, left, right)):
            continue
        pair_ticks += right - left
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe_midpoint(speed, grid, left, right):
                failure |= 1 << bit
        if failure.bit_count() <= 9:
            classes[failure] += right - left
    class_fnv = OFFSET
    for failure, width in sorted(classes.items()):
        class_fnv = fnv_add(class_fnv, failure)
        class_fnv = fnv_add(class_fnv, width)
    return grid, len(ordered) - 1, pair_ticks, dict(classes), class_fnv


def read_failures(path: Path) -> dict[int, list[int]]:
    result = {50: [], 96: []}
    distinct = {50: set(), 96: set()}
    body_fnv = {50: OFFSET, 96: OFFSET}
    global_fnv = OFFSET
    with path.open(newline="", encoding="ascii") as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in rows:
            q, r = int(row["q"]), int(row["r"])
            body = int(row["body_hex"], 16)
            require(q in result and r == ENDPOINT and body.bit_count() == 9,
                    "malformed failure row")
            require(body not in distinct[q], "duplicate failure body")
            distinct[q].add(body)
            result[q].append(body)
            body_fnv[q] = fnv_add(body_fnv[q], body)
            for word in (q, r, body):
                global_fnv = fnv_add(global_fnv, word)
    require(len(result[50]) == 20025 and len(result[96]) == 11,
            "failure counts changed")
    require(body_fnv == {50: 0xFF421454F02D9099,
                         96: 0x6F70A2D28FF6A957},
            "row failure FNV changed")
    require(global_fnv == 0x2228C33B4EF3B01D,
            "global failure FNV changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("failures", type=Path)
    parser.add_argument("detail_csv", type=Path)
    args = parser.parse_args()
    failures = read_failures(args.failures)
    geometries = {q: build_geometry(q) for q in failures}
    summaries: dict[int, tuple[int, int, int, int, int]] = {}
    with args.detail_csv.open("w", newline="", encoding="ascii") as handle:
        output = csv.writer(handle, lineterminator="\n")
        output.writerow(("q", "r", "ordinal", "body_hex",
                         "truncated_mass", "scaled_ticks"))
        for q in (50, 96):
            grid, _, _, classes, _ = geometries[q]
            ledger = OFFSET
            values: list[tuple[int, int]] = []
            for ordinal, body in enumerate(failures[q]):
                mass = sum(width for failure, width in classes.items()
                           if failure & body == 0)
                ticks = 63 * mass - 4 * grid
                output.writerow((q, ENDPOINT, ordinal, f"{body:08x}",
                                 mass, ticks))
                for word in (q, ENDPOINT, ordinal, body, mass):
                    ledger = fnv_add(ledger, word)
                ledger = fnv_i128(ledger, ticks)
                values.append((ticks, body))
            minimum = min(values, key=lambda item: (item[0], item[1]))
            maximum = max(values, key=lambda item: (item[0], -item[1]))
            positive = sum(ticks > 0 for ticks, _ in values)
            summaries[q] = (minimum[0], minimum[1], maximum[0],
                            maximum[1], ledger)
            require(positive == len(values), f"nonpositive literal at q={q}")

    require(geometries[50][0] == 2827379709554400 and
            geometries[50][1] == 8389 and
            geometries[50][2] == 2077256602813392 and
            len(geometries[50][3]) == 2383 and
            geometries[50][4] == 0x88D3EB2D7A477232,
            "q50 geometry changed")
    require(geometries[96][0] == 1130951883821760 and
            geometries[96][1] == 8501 and
            geometries[96][2] == 830901383902590 and
            len(geometries[96][3]) == 2352 and
            geometries[96][4] == 0xCB74E72091A68363,
            "q96 geometry changed")
    require(summaries[50] == (14566818763788984, 0x013C6401,
                              26685137010259728, 0x21884443,
                              0x9A0CD88B508499A2) and
            summaries[96] == (7172391058639758, 0x0D0C6401,
                              8558758749944214, 0x35126400,
                              0x313A0CDBA0E5AC5C),
            "minimum literal surplus changed")

    print("LRC14_ENDPOINT589_DIRECT_LITERAL_INDEPENDENT_V1")
    for q in (50, 96):
        grid, cells, pair_ticks, classes, class_fnv = geometries[q]
        minimum, min_body, maximum, max_body, ledger = summaries[q]
        print(f"ROW {q},{ENDPOINT} BODIES {len(failures[q])} GRID {grid} "
              f"CELLS {cells} PAIR_TICKS {pair_ticks} LOW_CLASSES "
              f"{len(classes)} CLASS_FNV {class_fnv:016x} POSITIVE "
              f"{len(failures[q])} MIN_TICKS {minimum} MIN_BODY "
              f"{min_body:08x} MAX_TICKS {maximum} MAX_BODY {max_body:08x} "
              f"DETAIL_FNV {ledger:016x}")
    print("TOTAL_BODIES 20036 ALL_STRICTLY_POSITIVE 1")
    print("INEQUALITY TRUNCATED_LOW_RANK_CELL_MASS_LE_LITERAL_SAFE_MASS")
    print("SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT589_FAILURE_BODIES_"
          "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
