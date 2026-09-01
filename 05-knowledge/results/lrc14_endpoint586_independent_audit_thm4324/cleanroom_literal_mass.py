#!/usr/bin/env python3
"""Clean-room raw-cell literal-mass audit for endpoint-586 failures."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


POOL = (8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
        80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
        168, 170, 176, 190, 193, 240, 252, 264, 286, 290)
OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
MASK128 = (1 << 128) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def add_word(state: int, word: int) -> int:
    word &= MASK64
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 255)) * PRIME) & MASK64
    return state


def fnv_words(words) -> int:
    state = OFFSET
    for word in words:
        state = add_word(state, word)
    return state


def read_failures(path: Path) -> list[int]:
    lines = path.read_text(encoding="ascii").splitlines()
    require(lines and lines[0] == "q,r,body_hex", "bad failure header")
    bodies: list[int] = []
    global_words: list[int] = []
    for line in lines[1:]:
        fields = line.split(",")
        require(len(fields) == 3, "malformed failure row")
        q, r, body = int(fields[0]), int(fields[1]), int(fields[2], 16)
        require(q == 50 and r == 586 and body.bit_count() == 9,
                "failure escaped q50 endpoint586 rank9 universe")
        bodies.append(body)
        global_words.extend((q, r, body))
    require(len(bodies) == 4090 and bodies == sorted(set(bodies)),
            "failure count/order changed")
    require(fnv_words(bodies) == 0x14CE094F4AB4BA94 and
            fnv_words(global_words) == 0xFFB884B2B17E6EF4,
            "failure FNV changed")
    return bodies


def safe(speed: int, grid: int, left: int, right: int) -> bool:
    phase2 = speed * (left + right) % (2 * grid)
    return 7 * phase2 >= grid and 7 * phase2 <= 13 * grid


def make_cells() -> tuple[int, int, list[tuple[int, int, bool]], int,
                          list[tuple[int, int]], list[tuple[int, int]]]:
    grid = 1
    for speed in (*POOL, 50, 586):
        grid = math.lcm(grid, 14 * speed)
    walls = {0, grid}
    for speed in (*POOL, 50, 586):
        unit = grid // (14 * speed)
        require(unit * 14 * speed == grid, "nonintegral wall unit")
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    cuts = sorted(walls)
    raw: list[tuple[int, int, bool]] = []
    pair_mass = 0
    classes: dict[int, int] = {}
    for left, right in zip(cuts, cuts[1:]):
        if not safe(50, grid, left, right) or not safe(586, grid, left, right):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe(speed, grid, left, right):
                failure |= 1 << bit
        width = right - left
        raw.append((failure, width, failure.bit_count() <= 9))
        pair_mass += width
        classes[failure] = classes.get(failure, 0) + width
    all_classes = sorted(classes.items())
    low_classes = [item for item in all_classes if item[0].bit_count() <= 9]
    return grid, len(cuts) - 1, raw, pair_mass, all_classes, low_classes


def class_fnv(classes: list[tuple[int, int]]) -> int:
    return fnv_words(word for row in classes for word in row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--detail", type=Path, required=True)
    args = parser.parse_args()
    bodies = read_failures(args.failures)
    grid, open_cells, raw, pair_mass, all_classes, low_classes = make_cells()
    require(grid == 26723298545143200 and open_cells == 8381 and
            len(raw) == 5802, "wall geometry changed")

    output = ["q,r,ordinal,body_hex,truncated_mass,full_mass,"
              "truncated_scaled_ticks,full_scaled_ticks\n"]
    detail_state = OFFSET
    low_positive = full_positive = equality = 0
    low_min = full_min = low_max = full_max = None
    low_min_body = full_min_body = low_max_body = full_max_body = 0
    for ordinal, body in enumerate(bodies):
        low_mass = full_mass = 0
        for failed, width, low in raw:
            if failed & body:
                continue
            full_mass += width
            if low:
                low_mass += width
        require(0 <= low_mass <= full_mass, "truncation exceeded full mass")
        low_ticks = 63 * low_mass - 4 * grid
        full_ticks = 63 * full_mass - 4 * grid
        output.append(f"50,586,{ordinal},{body:08x},{low_mass},{full_mass},"
                      f"{low_ticks},{full_ticks}\n")
        words = (50, 586, ordinal, body, low_mass, full_mass,
                 low_ticks & MASK64, (low_ticks & MASK128) >> 64,
                 full_ticks & MASK64, (full_ticks & MASK128) >> 64)
        for word in words:
            detail_state = add_word(detail_state, word)
        low_positive += low_ticks > 0
        full_positive += full_ticks > 0
        equality += low_mass == full_mass
        if low_min is None or (low_ticks, body) < (low_min, low_min_body):
            low_min, low_min_body = low_ticks, body
        if full_min is None or (full_ticks, body) < (full_min, full_min_body):
            full_min, full_min_body = full_ticks, body
        if low_max is None or low_ticks > low_max or (
                low_ticks == low_max and body < low_max_body):
            low_max, low_max_body = low_ticks, body
        if full_max is None or full_ticks > full_max or (
                full_ticks == full_max and body < full_max_body):
            full_max, full_max_body = full_ticks, body
    args.detail.write_bytes("".join(output).encode("ascii"))
    detail_sha = hashlib.sha256(args.detail.read_bytes()).hexdigest()
    print("LRC14_ENDPOINT586_CLEANROOM_RAWCELL_LITERAL_MASS_V1")
    print("FAILURE_ROWS 1 BODIES 4090 FAILURE_FNV ffb884b2b17e6ef4")
    print(f"ROW 50,586 BODIES 4090 BODY_FNV 14ce094f4ab4ba94 GRID {grid} "
          f"CELLS {open_cells} PAIR_SAFE_CELLS {len(raw)} PAIR_TICKS {pair_mass} "
          f"ALL_CLASSES {len(all_classes)} LOW_CLASSES {len(low_classes)} "
          f"ALL_CLASS_FNV {class_fnv(all_classes):016x} "
          f"LOW_CLASS_FNV {class_fnv(low_classes):016x} "
          f"LOW_POSITIVE {low_positive} FULL_POSITIVE {full_positive} "
          f"EQUALITY {equality} LOW_MIN_TICKS {low_min} "
          f"LOW_MIN_BODY {low_min_body:08x} FULL_MIN_TICKS {full_min} "
          f"FULL_MIN_BODY {full_min_body:08x} LOW_MAX_TICKS {low_max} "
          f"LOW_MAX_BODY {low_max_body:08x} FULL_MAX_TICKS {full_max} "
          f"FULL_MAX_BODY {full_max_body:08x} DETAIL_FNV {detail_state:016x}")
    print(f"DETAIL_SHA256 {detail_sha}")
    print(f"SUMMARY BODIES 4090 LOW_POSITIVE {low_positive} "
          f"FULL_POSITIVE {full_positive} EQUALITY {equality}")
    print("REPRESENTATION RAW_PAIR_SAFE_WALL_CELLS_NO_CLASS_AGGREGATION")
    print("INEQUALITY TRUNCATED_RANK_LE9_MASS_LE_FULL_LITERAL_MASS")
    print("SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT586_FAILURE_BODIES_"
          "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT " + ("PASS" if low_positive == len(bodies) else "HOSTILE_FAIL"))


if __name__ == "__main__":
    main()
