#!/usr/bin/env python3
"""Exact content-tower and pigeonhole sidecar for THM-4160.

Two different censuses must agree: a literal loop over all 11,100,375
labelled bodies, and a grouped min/max reconstruction using binomial counts.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)
NEWCOMERS = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)
REPAIR_SETS = {
    5: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
    66: (42, 80, 85, 88, 95, 145, 168, 170, 176, 193, 240, 252, 264, 286, 290),
    182: (85, 95, 145, 168, 176, 193, 240, 252, 290),
    298: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
    336: (85, 88, 95, 145, 168, 176, 193, 252, 264),
    340: (85, 88, 95, 145, 168, 240, 252, 264),
    380: (8, 16, 85, 88, 95, 132, 145, 168, 170, 176, 193, 240, 252, 264, 286, 290),
    386: (8, 16, 42, 85, 88, 95, 132, 145, 168, 170, 176, 190, 193, 240, 252, 264, 286, 290),
    528: (8, 85, 88, 95, 145, 168, 170, 193, 240, 252, 286, 290),
    572: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
}


def require(predicate: bool, label: object) -> None:
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def classify(defect: int) -> str:
    if defect <= 0:
        return "stable"
    if defect <= 13:
        return "transient"
    return "beyond"


def empty_census() -> dict[str, object]:
    return {
        "count": 0,
        "defects": Counter(),
        "maximums": Counter(),
        "stable_maximums": Counter(),
        "transient_maximums": Counter(),
        "beyond_maximums": Counter(),
    }


def record(census: dict[str, object], minimum: int, maximum: int, count: int = 1) -> None:
    defect = 16 * maximum - 156 * minimum
    category = classify(defect)
    census["count"] += count
    census["defects"][defect] += count
    census["maximums"][maximum] += count
    census[f"{category}_maximums"][maximum] += count


def grouped_fixed_family(fixed: tuple[int, ...], choice_size: int) -> dict[str, object]:
    """Binomial min/max reconstruction, with no literal subset loop."""
    answer = empty_census()
    fixed_minimum = min(fixed)
    fixed_maximum = max(fixed)
    minimum_candidates = tuple(value for value in OPTIONAL if value < fixed_minimum) + (fixed_minimum,)
    maximum_candidates = (fixed_maximum,) + tuple(
        value for value in OPTIONAL if value > fixed_maximum
    )
    for minimum in minimum_candidates:
        for maximum in maximum_candidates:
            forced = int(minimum != fixed_minimum) + int(maximum != fixed_maximum)
            interior = sum(minimum < value < maximum for value in OPTIONAL)
            count = choose(interior, choice_size - forced)
            if count:
                record(answer, minimum, maximum, count)
    return answer


def add_census(target: dict[str, object], source: dict[str, object]) -> None:
    target["count"] += source["count"]
    for key in (
        "defects", "maximums", "stable_maximums", "transient_maximums",
        "beyond_maximums",
    ):
        target[key].update(source[key])


def count_split(census: dict[str, object]) -> tuple[int, int, int]:
    defects = census["defects"]
    return (
        sum(count for defect, count in defects.items() if defect <= 0),
        sum(count for defect, count in defects.items() if 1 <= defect <= 13),
        sum(count for defect, count in defects.items() if defect >= 14),
    )


def asymptotic_constants(census: dict[str, object]) -> tuple[Fraction, Fraction, Fraction]:
    def constant(histogram: Counter[int]) -> Fraction:
        return sum(
            (Fraction(count, maximum) for maximum, count in histogram.items()),
            Fraction(0),
        )

    return (
        constant(census["maximums"]),
        constant(census["stable_maximums"]),
        constant(census["beyond_maximums"]),
    )


def ordered(counter: Counter[int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(counter.items()))


def main() -> None:
    require(len(OPTIONAL) == 27, "optional universe")
    require(set(NEWCOMERS).isdisjoint(POOL), "newcomers outside pool")

    literal_old = empty_census()
    old_thm4148_pass = 0
    for choice in combinations(OPTIONAL, 8):
        body = ANCHORS + choice
        minimum = min(body)
        maximum = max(body)
        record(literal_old, minimum, maximum)
        old_thm4148_pass += int(
            27 * (13 * minimum - maximum) >= 4 * minimum * maximum
        )

    literal_new = empty_census()
    literal_new_by_q = {newcomer: empty_census() for newcomer in NEWCOMERS}
    minimum_omitted_repairs = {newcomer: len(OPTIONAL) for newcomer in NEWCOMERS}
    new_thm4148_pass = 0
    for choice in combinations(OPTIONAL, 7):
        choice_set = frozenset(choice)
        for newcomer in NEWCOMERS:
            body = ANCHORS + (newcomer,) + choice
            minimum = min(body)
            maximum = max(body)
            record(literal_new, minimum, maximum)
            record(literal_new_by_q[newcomer], minimum, maximum)
            new_thm4148_pass += int(
                27 * (13 * minimum - maximum) >= 4 * minimum * maximum
            )
            omitted = len(set(REPAIR_SETS[newcomer]) - choice_set)
            require(omitted >= 1, ("pigeonhole witness", newcomer, choice))
            minimum_omitted_repairs[newcomer] = min(
                minimum_omitted_repairs[newcomer], omitted
            )

    require(literal_old["count"] == comb(27, 8) == 2_220_075, "old count")
    require(literal_new["count"] == 10 * comb(27, 7) == 8_880_300, "new count")
    require(old_thm4148_pass == new_thm4148_pass == 0, "THM-4148 hostile gate")
    require(minimum_omitted_repairs == {
        newcomer: len(REPAIR_SETS[newcomer]) - 7 for newcomer in NEWCOMERS
    }, "literal pigeonhole minima")

    grouped_old = grouped_fixed_family(ANCHORS, 8)
    grouped_new = empty_census()
    grouped_new_by_q = {}
    for newcomer in NEWCOMERS:
        grouped = grouped_fixed_family(ANCHORS + (newcomer,), 7)
        grouped_new_by_q[newcomer] = grouped
        add_census(grouped_new, grouped)

    for label, literal, grouped in (
        ("old", literal_old, grouped_old),
        ("new", literal_new, grouped_new),
    ):
        require(literal["count"] == grouped["count"], (label, "count"))
        for key in (
            "defects", "maximums", "stable_maximums", "transient_maximums",
            "beyond_maximums",
        ):
            require(literal[key] == grouped[key], (label, key))
    for newcomer in NEWCOMERS:
        for key in (
            "count", "defects", "maximums", "stable_maximums",
            "transient_maximums", "beyond_maximums",
        ):
            require(literal_new_by_q[newcomer][key] == grouped_new_by_q[newcomer][key],
                    ("newcomer grouped agreement", newcomer, key))

    old_split = count_split(literal_old)
    new_split = count_split(literal_new)
    require(old_split == (344_366, 0, 1_875_709), "old content split")
    require(new_split == (1_052_735, 0, 7_827_565), "new content split")
    total_split = tuple(old_split[index] + new_split[index] for index in range(3))
    require(total_split == (1_397_101, 0, 9_703_274), "total content split")

    expected_new_by_q = {
        5: (0, 0, 888030),
        66: (187639, 0, 700391),
        182: (182920, 0, 705110),
        298: (116280, 0, 771750),
        336: (116280, 0, 771750),
        340: (116280, 0, 771750),
        380: (116280, 0, 771750),
        386: (116280, 0, 771750),
        528: (50388, 0, 837642),
        572: (50388, 0, 837642),
    }
    require({
        newcomer: count_split(literal_new_by_q[newcomer])
        for newcomer in NEWCOMERS
    } == expected_new_by_q, "newcomer content splits")

    combined = empty_census()
    add_census(combined, literal_old)
    add_census(combined, literal_new)
    gaps = {}
    for label, census in (("old", literal_old), ("new", literal_new), ("all", combined)):
        defects = census["defects"]
        gaps[label] = (
            max(defect for defect in defects if defect <= 0),
            min(defect for defect in defects if defect > 0),
        )
    require(gaps == {"old": (-20, 192), "new": (-20, 88), "all": (-20, 88)},
            "content defect gaps")

    old_constants = asymptotic_constants(literal_old)
    new_constants = asymptotic_constants(literal_new)
    all_constants = asymptotic_constants(combined)
    require(all_constants == tuple(
        old_constants[index] + new_constants[index] for index in range(3)
    ), "combined asymptotic constants")
    require(all(
        values[0] == values[1] + values[2]
        for values in (old_constants, new_constants, all_constants)
    ), "stable/beyond asymptotic partitions")

    semantic = {
        "old_split": old_split,
        "new_split": new_split,
        "total_split": total_split,
        "new_by_q": expected_new_by_q,
        "gaps": gaps,
        "omitted": minimum_omitted_repairs,
        "old_maximums": ordered(literal_old["maximums"]),
        "new_maximums": ordered(literal_new["maximums"]),
        "all_maximums": ordered(combined["maximums"]),
        "old_stable_maximums": ordered(literal_old["stable_maximums"]),
        "new_stable_maximums": ordered(literal_new["stable_maximums"]),
        "all_stable_maximums": ordered(combined["stable_maximums"]),
        "old_beyond_maximums": ordered(literal_old["beyond_maximums"]),
        "new_beyond_maximums": ordered(literal_new["beyond_maximums"]),
        "all_beyond_maximums": ordered(combined["beyond_maximums"]),
        "constants": tuple(tuple(str(value) for value in values)
                           for values in (old_constants, new_constants, all_constants)),
    }

    print("LRC14_ANCHORED_DELETION_COVER_THM4160_CONTENT_TOWER_20260826")
    print("literal_universe=old:C(27,8);new:10*C(27,7);"
          "grouped_reconstruction=forced_min_max_plus_binomial_interior")
    print(f"literal_counts=old:{literal_old['count']},new:{literal_new['count']},"
          f"all:{combined['count']};grouped_counts_match=True")
    print(f"old_split={old_split};new_split={new_split};total_split={total_split}")
    print(f"new_splits_by_q={expected_new_by_q}")
    print(f"defect_gaps={gaps};transient_D_1_through_13=0")
    print(f"minimum_omitted_repairs={minimum_omitted_repairs};"
          "literal_pigeonhole_bodies_checked=8880300")
    for label, census, constants in (
        ("old", literal_old, old_constants),
        ("new", literal_new, new_constants),
        ("all", combined, all_constants),
    ):
        print(f"{label}_maximum_histogram={ordered(census['maximums'])}")
        print(f"{label}_stable_maximum_histogram={ordered(census['stable_maximums'])}")
        print(f"{label}_beyond_maximum_histogram={ordered(census['beyond_maximums'])}")
        print(f"{label}_asymptotic_constants=all:{constants[0]},"
              f"stable:{constants[1]},beyond:{constants[2]}")
    print("bounded_height_law=N(X)=C*X+O(11100375),"
          "with C given by each all/stable/beyond asymptotic constant")
    print("content_consequence=for every positive common dilation,"
          "stable count 1397101 and beyond-both-gates count 9703274")
    print("semantic_sha256=" + sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest())


if __name__ == "__main__":
    main()
