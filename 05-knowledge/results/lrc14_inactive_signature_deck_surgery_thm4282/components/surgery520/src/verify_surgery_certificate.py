#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from collections import defaultdict
from pathlib import Path


DELETED = [57, 107, 222, 275, 345]
WITNESS = [0x0090492C, 0x018C2114, 0x09A2A040,
           0x108C1112, 0x1C81A100, 0x38120016]
DUAL_DEN = 11
DUAL = {
    0: 2, 7: 1, 8: 3, 9: 3, 10: 6, 11: 1, 12: 2, 13: 1,
    14: 2, 16: 3, 17: 3, 18: 3, 20: 1, 24: 1, 25: 2,
    29: 2, 30: 4, 32: 2, 34: 1, 36: 2, 37: 1, 38: 3,
    42: 1, 43: 2, 49: 4, 50: 1,
}


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_masks(masks: list[int]) -> int:
    state = 0xCBF29CE484222325
    for mask in masks:
        for byte in range(8):
            state ^= (mask >> (8 * byte)) & 0xFF
            state = (state * 0x100000001B3) & ((1 << 64) - 1)
    return state


def read_deck(path: Path) -> list[int]:
    masks = [int(line, 16) for line in path.read_text().splitlines() if line]
    need(len(masks) == len(set(masks)), f"duplicate mask in {path}")
    need(all(mask.bit_count() == 8 and mask < (1 << 30) for mask in masks),
         f"bad mask in {path}")
    return masks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-deck", type=Path, required=True)
    parser.add_argument("--rebuilt-deck", type=Path, required=True)
    parser.add_argument("--lost", type=Path, required=True)
    parser.add_argument("--classes", type=Path, required=True)
    parser.add_argument("--active-masks", type=Path, required=True)
    parser.add_argument("--primary-common", type=Path, required=True)
    parser.add_argument("--literal-common", type=Path, required=True)
    args = parser.parse_args()

    original = read_deck(args.original_deck)
    rebuilt = read_deck(args.rebuilt_deck)
    need(len(original) == 421 and fnv_masks(original) == 0x20D63DD42FE8150E,
         "original deck identity changed")
    retained = [mask for index, mask in enumerate(original) if index not in DELETED]
    need(rebuilt == retained + WITNESS and len(rebuilt) == 422 and
         fnv_masks(rebuilt) == 0x813801C9BD1676BA,
         "rebuilt deck transformation/identity changed")

    lost_rows = args.lost.read_text().splitlines()
    need(lost_rows[0] == "ordinal,body_hex,deleted_response_hex" and
         len(lost_rows) == 54, "lost-body CSV shape changed")
    bodies: list[int] = []
    for expected, line in enumerate(lost_rows[1:]):
        ordinal, body_hex, response_hex = line.split(",")
        body = int(body_hex, 16)
        response = int(response_hex, 16)
        need(int(ordinal) == expected and body.bit_count() == 9 and
             0 < response < 32, "bad lost-body row")
        actual = 0
        for local, index in enumerate(DELETED):
            if body & original[index] == 0:
                actual |= 1 << local
        need(actual == response and
             all(body & original[index] for index in range(421)
                 if index not in DELETED),
             "lost-body response not confined to deletion")
        bodies.append(body)
    need(bodies == sorted(set(bodies)), "lost bodies not ordered/distinct")

    class_lines = args.classes.read_text().splitlines()
    need(class_lines[0] ==
         "class\tpattern_hex\tleast_mask_hex\tmultiplicity\tcover\tmaximal",
         "class header changed")
    classes: dict[int, tuple[int, int, int, int, int]] = {}
    for line in class_lines[1:]:
        fields = line.split("\t")
        need(len(fields) == 6, "bad response-class row")
        index = int(fields[0])
        pattern = int(fields[1], 16)
        least = int(fields[2], 16)
        multiplicity = int(fields[3])
        cover = int(fields[4])
        maximal = int(fields[5])
        need(index == len(classes) and cover == pattern.bit_count() and
             multiplicity > 0 and maximal in (0, 1), "bad class fields")
        classes[index] = (pattern, least, multiplicity, cover, maximal)
    need(len(classes) == 7124 and sum(row[4] for row in classes.values()) == 810,
         "response-class counts changed")

    counts: dict[int, int] = defaultdict(int)
    least_by_class: dict[int, int] = {}
    witness_patterns: dict[int, int] = {}
    active_lines = args.active_masks.read_text().splitlines()
    need(active_lines[0] == "colex_rank\tmask_hex\tclass\tpattern_hex" and
         len(active_lines) == 2879148, "active-mask TSV shape changed")
    previous_rank = -1
    for line in active_lines[1:]:
        rank_s, mask_s, class_s, pattern_s = line.split("\t")
        rank = int(rank_s)
        mask = int(mask_s, 16)
        class_index = int(class_s)
        pattern = int(pattern_s, 16)
        need(rank > previous_rank and mask.bit_count() == 8,
             "active-mask order/rank changed")
        previous_rank = rank
        need(class_index in classes and classes[class_index][0] == pattern,
             "active mask mapped to wrong class")
        counts[class_index] += 1
        least_by_class[class_index] = min(least_by_class.get(class_index, mask),
                                          mask)
        if mask in WITNESS:
            witness_patterns[mask] = pattern
    for class_index, (_, least, multiplicity, _, _) in classes.items():
        need(counts[class_index] == multiplicity and
             least_by_class[class_index] == least,
             "class multiplicity/least representative mismatch")

    dual_sum = sum(DUAL.values())
    maximum_class_weight = 0
    for pattern, _, _, _, _ in classes.values():
        weight = sum(value for obligation, value in DUAL.items()
                     if (pattern >> obligation) & 1)
        maximum_class_weight = max(maximum_class_weight, weight)
        need(weight <= DUAL_DEN, "dual class inequality failed")
    need(dual_sum == 57 and maximum_class_weight == 11,
         "dual objective/boundary changed")

    need(set(witness_patterns) == set(WITNESS), "witness mask missing")
    witness_union = 0
    for mask in WITNESS:
        witness_union |= witness_patterns[mask]
    need(witness_union == (1 << 53) - 1,
         "six-mask witness does not cover all obligations")

    need(args.primary_common.read_bytes() == args.literal_common.read_bytes(),
         "primary/literal common CSV disagreement")
    need(len(args.primary_common.read_text().splitlines()) == 145122,
         "rebuilt common count changed")

    print("THM4281_SIGNATURE_SURGERY_DETACHED_CERTIFICATE_V1")
    print("ORIGINAL_DECK", len(original), f"FNV {fnv_masks(original):016x}")
    print("REBUILT_DECK", len(rebuilt), f"FNV {fnv_masks(rebuilt):016x}")
    print("LOST_OBLIGATIONS", len(bodies), "SHA256", sha(args.lost))
    print("ACTIVE_MASKS", len(active_lines) - 1, "SHA256", sha(args.active_masks))
    print("RESPONSE_CLASSES", len(classes), "MAXIMAL", 810,
          "SHA256", sha(args.classes))
    print("WITNESS", len(WITNESS),
          "UNION", f"{witness_union:014x}")
    print("DUAL NUMERATOR_SUM", dual_sum, "DENOMINATOR", DUAL_DEN,
          "MAX_CLASS_NUMERATOR", maximum_class_weight, "CEILING", 6)
    print("COMMON", 145122, "SHA256", sha(args.primary_common),
          "PRIMARY_LITERAL_IDENTICAL 1")
    print("VERDICT PASS DETACHED_EXACT_LEDGER_DUAL_AND_WITNESS_AUDIT")


if __name__ == "__main__":
    main()
