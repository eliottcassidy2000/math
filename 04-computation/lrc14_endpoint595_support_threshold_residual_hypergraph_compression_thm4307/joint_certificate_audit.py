#!/usr/bin/env python3
"""Independent integer audit of the D=1667 dual and frozen ten-mask cover."""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

from solve_response_cover import BIT_POSITIONS, FROZEN_MIXED_ALL_DUAL

MASK64 = (1 << 64) - 1
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


class Fnv:
    def __init__(self):
        self.state = FNV_OFFSET

    def add(self, value: int):
        value &= MASK64
        for byte in range(8):
            self.state ^= value >> (8 * byte) & 0xFF
            self.state = self.state * FNV_PRIME & MASK64


def run(atlas_path: Path, additions_path: Path):
    weights = [0] * 145
    for index, value in FROZEN_MIXED_ALL_DUAL.items():
        weights[index] = value
    if sum(weights) != 15030:
        raise AssertionError("dual total changed")

    by_mask: dict[int, tuple[int, tuple[int, int, int, int]]] = {}
    histogram = Counter()
    load_ledger = Fnv()
    type_count = 0
    with atlas_path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            words = tuple(int(row[f"w{i}"], 16) for i in range(4))
            load = 0
            for logical, physical in enumerate(BIT_POSITIONS):
                word = physical // 64
                bit = physical % 64
                if words[word] >> bit & 1:
                    load += weights[logical]
            if load > 1667:
                raise AssertionError("dual overload")
            histogram[load] += 1
            for word in words:
                load_ledger.add(word)
            load_ledger.add(load)
            type_count += 1
            if row["least8"]:
                by_mask[int(row["least8"], 16)] = (8, words)
            if row["least9"]:
                by_mask[int(row["least9"], 16)] = (9, words)
    if type_count != 14619:
        raise AssertionError("response type count changed")

    additions = [int(line.strip(), 16) for line in additions_path.read_text().splitlines() if line.strip()]
    if len(additions) != 10 or len(set(additions)) != 10:
        raise AssertionError("addition ledger changed")
    union = [0, 0, 0, 0]
    addition_ledger = Fnv()
    cover_ledger = Fnv()
    for mask in additions:
        rank, words = by_mask[mask]
        if rank != 9:
            raise AssertionError("frozen addition is not rank nine")
        addition_ledger.add(mask)
        cover_ledger.add(mask)
        for index, word in enumerate(words):
            union[index] |= word
            cover_ledger.add(word)
    expected = [(1 << 64) - 1, (1 << 52) - 1, (1 << 13) - 1, (1 << 16) - 1]
    if union != expected:
        raise AssertionError("ten-mask response union incomplete")

    histogram_ledger = Fnv()
    for load, count in sorted(histogram.items()):
        histogram_ledger.add(load)
        histogram_ledger.add(count)
    print("LRC14_ENDPOINT595_JOINT_CERTIFICATE_AUDIT_V1")
    print("DUAL DENOMINATOR 1667 TOTAL 15030 NINE_TIMES_DENOMINATOR 15003 INTEGER_LOWER 10")
    print(
        f"DUAL_TYPES {type_count} MAX_LOAD {max(histogram)} "
        f"SATURATED_TYPES {histogram[1667]} DISTINCT_LOADS {len(histogram)} "
        f"LOAD_LEDGER_FNV {load_ledger.state:x} HISTOGRAM_FNV {histogram_ledger.state:x}"
    )
    print(
        f"COVER ADDITIONS 10 RANK9 10 MASK_FNV {addition_ledger.state:x} "
        f"MASK_RESPONSE_FNV {cover_ledger.state:x} UNION "
        + ":".join(f"{word:016x}" for word in union)
        + " EXACT_MINIMUM 10"
    )
    print("DUAL_SUPPORT " + " ".join(f"{i}:{v}" for i, v in sorted(FROZEN_MIXED_ALL_DUAL.items())))
    print("SCOPE COMPLETE_14619_TYPE_MIXED_RESPONSE_QUOTIENT_EXACT_INTEGER_AUDIT")
    print("VERDICT PASS")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: joint_certificate_audit.py ATLAS.csv ADDITIONS10.txt")
    run(Path(sys.argv[1]), Path(sys.argv[2]))
