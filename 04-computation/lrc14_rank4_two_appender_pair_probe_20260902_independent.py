#!/usr/bin/env python3
"""Independent exact wall/candidate audit for the rank-four pair probe.

This implementation imports no code from the C++ census.  It rebuilds the
literal wall arrangement with an unfurled two-sided midpoint test and checks
the full rank ledger, the rank-four-minimizing candidates reported by the
exhaustive census, their retained mass, and their complete safe mass.  It
does not independently enumerate all C(30,8) bodies; optimizer independence
is supplied by the direct-cell branch-and-bound in the primary program.
"""

from __future__ import annotations

import csv
import hashlib
import math
from fractions import Fraction
from pathlib import Path


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
REMAINDER = Path(
    "05-knowledge/results/"
    "lrc14_rank2_wall_graph_complete_pair_closure_thm4326/"
    "thm4231_remainder181194.csv"
)
REMAINDER_SHA256 = (
    "9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1"
)

EXPECTED = {
    (50, 70): {
        "grid": 91_205_797_082_400,
        "pair_safe": 67_008_340_713_600,
        "rank_mass": (
            2_736_428_937_434,
            3_624_615_539_184,
            5_940_679_612_750,
            9_782_779_215_280,
            12_121_120_132_736,
        ),
        "rank3": 22_084_503_304_648,
        "body": 0x003CC402,
        "retained": 11_108_019_386_090,
        "full": 14_349_106_369_428,
        "ticks": 261_308_990_696_490,
        "component_bound": 2_011,
        "appender_cutoff": 6_021,
    },
    (509, 640): {
        "grid": 74_278_001_143_906_560,
        "pair_safe": 54_571_555_450_299_888,
        "rank_mass": (
            3_425_342_246_785_080,
            4_228_127_642_044_974,
            5_663_341_764_242_158,
            7_737_111_571_740_786,
            8_694_506_884_621_480,
        ),
        "rank3": 21_053_923_224_812_998,
        "body": 0x10386401,
        "retained": 11_595_224_392_491_506,
        "full": 14_171_259_234_541_482,
        "ticks": 419_267_167_784_466_066,
        "component_bound": 3_040,
        "appender_cutoff": 5_295,
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def midpoint_safe(speed: int, grid: int, left: int, right: int) -> bool:
    # If z=2D*{speed*(left+right)/(2D)}, safety is exactly
    # D <= 7z <= 13D.  This deliberately avoids the C++ fold-to-half-circle
    # implementation.
    residue = speed * (left + right) % (2 * grid)
    return grid <= 7 * residue <= 13 * grid


def audit_pair(q: int, r: int) -> str:
    expected = EXPECTED[(q, r)]
    grid = 1
    for speed in POOL + (q, r):
        grid = math.lcm(grid, 14 * speed)
    require(grid == expected["grid"], f"{q},{r}: grid changed")

    walls = {0, grid}
    for speed in POOL + (q, r):
        unit = grid // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)

    rank_mass = [0] * 31
    rank_cells = [0] * 31
    pair_safe = 0
    retained_candidate = 0
    full_candidate = 0
    body = expected["body"]
    for left, right in zip(ordered, ordered[1:]):
        if not midpoint_safe(q, grid, left, right):
            continue
        if not midpoint_safe(r, grid, left, right):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not midpoint_safe(speed, grid, left, right):
                failure |= 1 << bit
        width = right - left
        rank = failure.bit_count()
        pair_safe += width
        rank_mass[rank] += width
        rank_cells[rank] += 1
        if failure & body == 0:
            full_candidate += width
            if rank <= 4:
                retained_candidate += width

    require(pair_safe == expected["pair_safe"], f"{q},{r}: pair-safe mass")
    require(
        tuple(rank_mass[:5]) == expected["rank_mass"],
        f"{q},{r}: rank-zero-through-four ledger",
    )
    require(sum(rank_mass[:4]) == expected["rank3"], f"{q},{r}: rank3 prefix")
    require(sum(rank_mass) == pair_safe, f"{q},{r}: full rank partition")
    require(
        retained_candidate == expected["retained"],
        f"{q},{r}: candidate retained mass",
    )
    require(full_candidate == expected["full"], f"{q},{r}: candidate full mass")
    ticks = 81 * retained_candidate - 7 * grid
    require(ticks == expected["ticks"] and ticks > 0, f"{q},{r}: target ticks")
    minimum_ratio = Fraction(retained_candidate, grid)
    gamma = Fraction(6, 7) * minimum_ratio - Fraction(4, 63)
    require(gamma > 0, f"{q},{r}: one-appender limiting surplus")
    component_bound = expected["component_bound"]
    threshold = Fraction(6 * component_bound, 49) / gamma
    cutoff = (threshold.numerator + threshold.denominator - 1) // threshold.denominator
    require(cutoff == expected["appender_cutoff"], f"{q},{r}: appender cutoff")
    lower_at_cutoff = (
        Fraction(6, 7) * minimum_ratio
        - Fraction(6 * component_bound, 49 * cutoff)
    )
    require(lower_at_cutoff > Fraction(4, 63), f"{q},{r}: cutoff strictness")
    labels = ",".join(str(POOL[i]) for i in range(30) if body >> i & 1)
    return (
        f"PAIR {q},{r} GRID {grid} RANK_CELLS_0_TO_4 "
        + ",".join(str(value) for value in rank_cells[:5])
        + f" CANDIDATE {body:08x} LABELS {labels} "
        + f"RETAINED {retained_candidate} FULL {full_candidate} "
        + f"TICKS {ticks} STATUS MATCH"
        + f"\nAPPENDER {q},{r} COMPONENT_BOUND {component_bound} GAMMA {gamma} "
        + f"THRESHOLD {threshold} CUTOFF {cutoff} "
        + f"LOWER_AT_CUTOFF {lower_at_cutoff} "
        + f"SURPLUS {lower_at_cutoff-Fraction(4, 63)} STATUS STRICT_PASS"
    )


def main() -> None:
    raw = REMAINDER.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == REMAINDER_SHA256, "THM-4231 remainder hash changed")
    wanted = {("50", "70"), ("509", "640")}
    found: list[tuple[str, str]] = []
    with REMAINDER.open(newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and tuple(row[:2]) in wanted:
                found.append((row[0], row[1]))
    require(found == [("50", "70"), ("509", "640")], "residual membership")

    print("LRC14_RANK4_TWO_APPENDER_PAIR_INDEPENDENT_AUDIT_V1")
    print(f"REMAINDER_SHA256 {digest} MEMBERSHIP 50,70;509,640 STATUS MATCH")
    for pair in EXPECTED:
        print(audit_pair(*pair))
    print("METHOD CLEANROOM_PYTHON_LITERAL_WALLS_TWO_SIDED_MIDPOINT_DIRECT_CANDIDATE_SUM")
    print("SCOPE WALL_LEDGER_AND_REPORTED_CANDIDATES_NOT_INDEPENDENT_GLOBAL_ENUMERATION")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
