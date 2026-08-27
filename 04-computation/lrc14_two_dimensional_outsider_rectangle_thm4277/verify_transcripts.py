#!/usr/bin/env python3
"""Cross-check the deterministic THM-4277 primary and literal transcripts."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lines(path: Path) -> list[str]:
    raw = path.read_bytes()
    require(b"\r" not in raw, f"{path}: CR byte found")
    require(raw.endswith(b"\n"), f"{path}: missing final LF")
    decoded = raw.decode("ascii").splitlines()
    require(not any(line.startswith("SECONDS ") for line in decoded),
            f"{path}: nondeterministic timing line was not removed")
    return decoded


def unique_prefix(rows: list[str], prefix: str) -> str:
    matches = [row for row in rows if row.startswith(prefix)]
    require(len(matches) == 1, f"prefix {prefix!r} occurs {len(matches)} times")
    return matches[0]


WEAKEST = re.compile(
    r"^WEAKEST_NORMALIZED_GAP_PAIR (?P<pair>[0-9]+,[0-9]+) "
    r"REPAIR (?P<mask>0x[0-9a-f]+) \{(?P<labels>[^}]*)\} "
    r"MASS (?P<mass>[0-9]+/[0-9]+) GAP (?P<gap>[0-9]+/[0-9]+) "
    r"RAW_MARGIN (?P<raw>[0-9]+)$"
)


def parse_weakest(row: str) -> dict[str, str]:
    match = WEAKEST.fullmatch(row)
    require(match is not None, "malformed weakest-cell line")
    return match.groupdict()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary", type=Path)
    parser.add_argument("literal", type=Path)
    args = parser.parse_args()
    primary = lines(args.primary)
    literal = lines(args.literal)

    require(len(primary) == 8, "primary transcript line count changed")
    require(len(literal) == 9, "literal transcript line count changed")
    require(primary[0] == "LRC14_2D_RECTANGLE_PRIMARY_V1",
            "primary header changed")
    require(literal[0] == "LRC14_2D_RECTANGLE_LITERAL_V1",
            "literal header changed")
    require(unique_prefix(literal, "LARGEST_LITERAL_GRID ") ==
            "LARGEST_LITERAL_GRID 5889213041088817440 AT 499,647",
            "largest literal grid changed")

    expected_common = {
        "RECT ":
            "RECT Q 450..499 R 600..650 PAIRS 2550 PAIR_FNV 32d6ca3501ec5d84",
        "CANDIDATES_8192_FNV ":
            "CANDIDATES_8192_FNV 60148ca1fc61dbcb "
            "CANDIDATES_16384_FNV adf20f0ef1cadc1f",
        "MATRIX_CELLS ":
            "MATRIX_CELLS 41779200 NONNEGATIVE 36657425 EQUALITIES 0 "
            "SIGN_BYTE_FNV 9d3e995e23a7695a",
        "HOSTILE_COMMON ":
            "HOSTILE_COMMON 2572 FNV 44a01e1ab114723e BODY_FAILURES 7 "
            "FIRST_FAILURE 0x51c7008 {16,85,88,95,143,145,168,193,252} "
            "CHECKS 508099830 MAX_CHECKS 2572",
        "THEOREM_COMMON ":
            "THEOREM_COMMON 5257 FNV 60f329212844f8ac BODY_FAILURES 0 "
            "CHECKS 508103822 MAX_CHECKS 4129 WORST_BODY 0x7187008 "
            "{16,85,88,95,145,168,193,240,252} WITNESS "
            "{63,80,120,143,190,264,286,290}",
        "VERDICT ": "VERDICT PASS",
    }
    for prefix, expected in expected_common.items():
        p_row = unique_prefix(primary, prefix)
        l_row = unique_prefix(literal, prefix)
        require(p_row == l_row, f"cross-engine mismatch at {prefix!r}")
        require(p_row == expected, f"frozen ledger changed at {prefix!r}")

    p_weak = parse_weakest(unique_prefix(primary, "WEAKEST_NORMALIZED_GAP_PAIR "))
    l_weak = parse_weakest(unique_prefix(literal, "WEAKEST_NORMALIZED_GAP_PAIR "))
    invariant_keys = ("pair", "mask", "labels", "mass", "gap")
    require(all(p_weak[key] == l_weak[key] for key in invariant_keys),
            "cross-engine weakest invariant mismatch")
    expected_weak = {
        "pair": "462,626",
        "mask": "0x6042229",
        "labels": "8,16,30,63,88,143,240,252",
        "mass": "16477591782853/259521949879920",
        "gap": "7663493/259521949879920",
    }
    require(all(p_weak[key] == value for key, value in expected_weak.items()),
            "corrected weakest invariant changed")
    require(int(p_weak["raw"]) > 0 and int(l_weak["raw"]) > 0,
            "weakest raw margin is not positive")

    print("LRC14_2D_RECTANGLE_TRANSCRIPT_CROSSCHECK_V1")
    print("REPORTED_PAIR_CANDIDATE_MATRIX_DECK_BODY_LEDGER_LINES BYTE_EQUAL PASS")
    print("WEAKEST_INVARIANTS " + " ".join(
        f"{key.upper()} {p_weak[key]}" for key in invariant_keys))
    print(f"RAW_MARGINS PRIMARY {p_weak['raw']} LITERAL {l_weak['raw']}")
    print("VERDICT PASS PRIMARY_LITERAL_EXACT_AGREEMENT")


if __name__ == "__main__":
    main()
