#!/usr/bin/env python3
"""Primary exact owner-word census for THM-4075."""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd


CLOCKS = tuple(range(15, 41))
MAXIMA = tuple(range(25, 31))
FULL_DIVISOR_MASK = (1 << 13) - 1
PRIOR_CANDIDATES = 59_880

EXPECTED_COUNTS = {
    25: (22_953, 22_953, 0),
    26: (112_241, 112_207, 34),
    27: (149_354, 149_354, 0),
    28: (457_165, 457_062, 103),
    29: (252_939, 252_939, 0),
    30: (769_704, 769_452, 252),
}

EXPECTED_HISTOGRAMS = {
    25: {16: 1107, 17: 4422, 19: 9411, 20: 767, 21: 1439, 22: 644,
         23: 3411, 24: 209, 26: 677, 27: 660, 28: 34, 29: 8, 30: 2,
         31: 108, 32: 45, 34: 9},
    26: {15: 2500, 16: 18531, 17: 26273, 18: 4547, 19: 26034,
         20: 6752, 21: 7219, 22: 4885, 23: 9980, 24: 441, 25: 3641,
         27: 759, 28: 511, 29: 72, 31: 51, 32: 19, 34: 23, 35: 3},
    27: {15: 2619, 16: 15834, 17: 43515, 18: 7833, 19: 36641,
         20: 4471, 21: 12747, 22: 3647, 23: 13944, 24: 1011, 25: 5137,
         26: 772, 28: 437, 29: 49, 30: 132, 31: 326, 32: 68, 33: 96,
         34: 57, 35: 10, 36: 2, 37: 6},
    28: {15: 18992, 16: 59048, 17: 122819, 18: 21469, 19: 102944,
         20: 11791, 21: 39747, 22: 14139, 23: 45834, 24: 1003,
         25: 13782, 26: 2005, 27: 1877, 29: 188, 30: 42, 31: 1010,
         32: 214, 33: 108, 34: 55, 35: 93, 37: 5},
    29: {15: 3204, 16: 29559, 17: 84513, 18: 4998, 19: 74216,
         20: 3373, 21: 16496, 22: 5899, 23: 20637, 24: 251, 25: 7298,
         26: 930, 27: 842, 28: 252, 31: 121, 32: 187, 33: 23, 34: 98,
         35: 33, 37: 9},
    30: {16: 119093, 17: 217092, 18: 36930, 19: 206344, 20: 39420,
         21: 72092, 22: 13344, 23: 45222, 24: 805, 25: 13989,
         26: 1704, 27: 2459, 28: 300, 29: 197, 31: 115, 32: 61,
         33: 412, 34: 35, 35: 70, 37: 16, 38: 4},
}

EXPECTED_HARDEST = (
    (1, 3, 11, 23, 24, 25, 26, 27, 28, 29, 30),
    (1, 3, 22, 23, 24, 25, 26, 27, 28, 29, 30),
    (2, 3, 11, 23, 24, 25, 26, 27, 28, 29, 30),
    (2, 3, 22, 23, 24, 25, 26, 27, 28, 29, 30),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def abs_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def build_tables(maximum: int) -> tuple[dict[int, list[int]], dict[int, list[tuple[int, int]]]]:
    safe_masks: dict[int, list[int]] = {}
    tail_tables: dict[int, list[tuple[int, int]]] = {}
    for clock in CLOCKS:
        masks = [0] * (maximum + 1)
        for speed in range(1, maximum + 1):
            mask = 0
            for residue in range(clock):
                if 14 * abs_residue(speed * residue, clock) >= clock:
                    mask |= 1 << residue
            masks[speed] = mask
        safe_masks[clock] = masks

        tails = []
        for tail in range(1, 2 * clock, 2):
            eligible = 0
            owner = 0
            for residue in range(clock):
                if 7 * abs_residue(tail * residue, clock) < clock:
                    eligible |= 1 << residue
                    nearest = (2 * tail * residue + clock) // (2 * clock)
                    if nearest % 2:
                        owner |= 1 << residue
            tails.append((eligible, owner))
        tail_tables[clock] = tails
    return safe_masks, tail_tables


def core_packet(core: tuple[int, ...], clock: int, safe_masks: dict[int, list[int]]) -> int:
    packet = (1 << clock) - 1
    for speed in core:
        packet &= safe_masks[clock][speed]
    return packet


def owner_word_certificate(packet: int, tails: list[tuple[int, int]]) -> bool:
    if packet == 0:
        return False
    words = {
        owner & packet
        for eligible, owner in tails
        if packet & ~eligible == 0
    }
    return all((packet ^ word) not in words for word in words)


def divisor_masks(maximum: int) -> list[int]:
    masks = [0] * (maximum + 1)
    for speed in range(1, maximum + 1):
        for divisor in range(2, 15):
            if speed % divisor == 0:
                masks[speed] |= 1 << (divisor - 2)
    return masks


def census() -> tuple[list[dict[str, object]], list[tuple[int, ...]]]:
    safe_masks, tail_tables = build_tables(max(MAXIMA))
    div_masks = divisor_masks(max(MAXIMA))
    caches: dict[int, dict[int, bool]] = {clock: {} for clock in CLOCKS}
    rows: list[dict[str, object]] = []
    hardest: list[tuple[int, ...]] = []

    for maximum in MAXIMA:
        searched = 0
        candidates = 0
        primitive = 0
        nonprimitive = 0
        histogram: dict[int, int] = {}
        uncertified = 0

        for base in combinations(range(1, maximum), 10):
            searched += 1
            core = base + (maximum,)
            cover = div_masks[maximum]
            for speed in base:
                cover |= div_masks[speed]
            if cover != FULL_DIVISOR_MASK:
                continue

            candidates += 1
            if gcd(*core) == 1:
                primitive += 1
            else:
                nonprimitive += 1

            first_clock = 0
            for clock in CLOCKS:
                packet = core_packet(core, clock, safe_masks)
                cached = caches[clock].get(packet)
                if cached is None:
                    cached = owner_word_certificate(packet, tail_tables[clock])
                    caches[clock][packet] = cached
                if cached:
                    first_clock = clock
                    histogram[clock] = histogram.get(clock, 0) + 1
                    break
            if first_clock == 0:
                uncertified += 1
            elif first_clock == 38:
                hardest.append(core)

        expected_candidates, expected_primitive, expected_nonprimitive = EXPECTED_COUNTS[maximum]
        require(searched == comb(maximum - 1, 10), f"searched count at {maximum}")
        require(candidates == expected_candidates, f"candidate count at {maximum}")
        require(primitive == expected_primitive, f"primitive count at {maximum}")
        require(nonprimitive == expected_nonprimitive, f"nonprimitive count at {maximum}")
        require(uncertified == 0, f"uncertified core at {maximum}")
        require(histogram == EXPECTED_HISTOGRAMS[maximum], f"histogram at {maximum}")
        rows.append({
            "maximum": maximum,
            "searched": searched,
            "candidates": candidates,
            "primitive": primitive,
            "nonprimitive": nonprimitive,
            "histogram": histogram,
        })

    require(tuple(hardest) == EXPECTED_HARDEST, "clock-38 boundary cores")
    return rows, hardest


def format_histogram(histogram: dict[int, int]) -> str:
    return ",".join(f"{clock}:{histogram[clock]}" for clock in sorted(histogram))


def main() -> None:
    rows, hardest = census()
    new_candidates = sum(int(row["candidates"]) for row in rows)
    semantic = {
        "prior_candidates": PRIOR_CANDIDATES,
        "rows": rows,
        "hardest": hardest,
        "new_candidates": new_candidates,
        "combined_candidates": PRIOR_CANDIDATES + new_candidates,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    require(new_candidates == 1_764_356, "new candidate total")
    require(PRIOR_CANDIDATES + new_candidates == 1_824_236, "combined candidate total")

    print("LRC14 DIVISOR-COMPLETE CLOSURE THROUGH 30 -- PRIMARY OWNER WORDS")
    print(f"inherited THM-2066 max<=24 candidates={PRIOR_CANDIDATES}")
    for row in rows:
        print(
            f"M={row['maximum']} searched={row['searched']} "
            f"divisor_complete={row['candidates']} primitive={row['primitive']} "
            f"nonprimitive={row['nonprimitive']} certified={row['candidates']} "
            f"hist={format_histogram(row['histogram'])}"
        )
    print(f"new_candidates={new_candidates} combined_candidates={PRIOR_CANDIDATES + new_candidates}")
    print("clock_38_boundary=" + ";".join(",".join(map(str, core)) for core in hardest))
    print(f"semantic_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
