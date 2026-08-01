#!/usr/bin/env python3
"""Exact projected k=2 scalar atlas on 1580..1599 and z1=1599 closure.

The all-3003-body atlas has 26 rows on seven occupied heights.  Its top
height z1=1599 consists of four HIGH-TAIL rows.  For each row the exact
one-high ray maximum, even after adding the three largest distinct low
labels, remains strictly below the scalar floor.  Thus the entire top shell
is scalar-empty uniformly over all later labels and the projected k=2 cap
drops from 1599 to 1598.

This is a finite-exact input for reserved THM-2995, not a promotion of that
stub and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
AUDIT_NAME = "lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980"
AUDIT_SOURCE = ROOT / f"04-computation/{AUDIT_NAME}.py"
AUDIT_OUTPUT = ROOT / f"05-knowledge/results/{AUDIT_NAME}.out"
ENGINE_NAME = "lrc14_j7_k2_z1608_z1604_terminal_descent_thm2980"
ENGINE_SOURCE = ROOT / f"04-computation/{ENGINE_NAME}.py"
ENGINE_OUTPUT = ROOT / f"05-knowledge/results/{ENGINE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1599_forced_high_closure_thm2995.out"

EXPECTED_AUDIT_SOURCE_SHA256 = "56a074db49eb30a5579522d34bf21fd434680ae2212843d49bf959931ae45bcd"
EXPECTED_AUDIT_OUTPUT_SHA256 = "433145c039ee083ae1784ca2e89f12e6179520367f92e0a18b2f5931dd769fa8"
EXPECTED_ENGINE_SOURCE_SHA256 = "23298a37076b7b291fa2ba5dbbe849ed98b1dff33e94be8f2e06cc5e47a83c8a"
EXPECTED_ENGINE_OUTPUT_SHA256 = "65d6d371e2ddd846167f9452b6e41c7f26535a4f4b2c70d20b39610194e750ea"
EXPECTED_ATLAS_PROFILE_SHA256 = "0f02d797f9a3bc39871d4a4ccf72c42525ca3ca4109773469b8fe95b780780d3"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "aab25b18f0bb9e2cb33eaa378176e720e1375c19891be9c97c67371f72ae2636"
EXPECTED_CLOSURE_PROFILE_SHA256 = "8945d0a0602c1782d5b57b2a281a92b3b8a179f4ac2260431acbb1ab45896d3e"
EXPECTED_SEMANTIC_SHA256 = "edf05b0f5fb03e8daa00ff7d30a72db6670d045dc035485b531bcdfd5e5cfa66"

START = 1580
END = 1599
EXPECTED_CANDIDATE_ROWS = 60_060
EXPECTED_HEIGHTS = ((1581, 3), (1586, 11), (1588, 1), (1590, 3), (1594, 2), (1595, 2), (1599, 4))
TOP_CASES = tuple((1599, body) for body in (
    (1, 2, 10, 12, 13, 14),
    (1, 3, 10, 12, 13, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
))
EXPECTED_GAPS = {
    TOP_CASES[0]: (F(110655592, 142563170295), F(3340381379, 14541443370090)),
    TOP_CASES[1]: (F(41353591, 57025268118), F(1398881597, 7128158514750)),
    TOP_CASES[2]: (F(9709149, 12672281804), F(6614593, 38016845412)),
    TOP_CASES[3]: (F(318065119, 450897468840), F(2183240303, 75299877296280)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(AUDIT_SOURCE) == EXPECTED_AUDIT_SOURCE_SHA256, "integrated audit source changed")
require(file_sha256(AUDIT_OUTPUT) == EXPECTED_AUDIT_OUTPUT_SHA256, "integrated audit output changed")
require(file_sha256(ENGINE_SOURCE) == EXPECTED_ENGINE_SOURCE_SHA256, "forced-high engine source changed")
require(file_sha256(ENGINE_OUTPUT) == EXPECTED_ENGINE_OUTPUT_SHA256, "forced-high engine output changed")
require("terminal_cap=1599" in AUDIT_OUTPUT.read_text(encoding="utf-8"), "inherited cap changed")

ASPEC = spec_from_file_location("k2_z1599_integrated_audit", AUDIT_SOURCE)
require(ASPEC is not None and ASPEC.loader is not None, "cannot load integrated atlas audit")
AUDIT = module_from_spec(ASPEC)
ASPEC.loader.exec_module(AUDIT)
ESPEC = spec_from_file_location("k2_z1599_forced_high_engine", ENGINE_SOURCE)
require(ESPEC is not None and ESPEC.loader is not None, "cannot load forced-high engine")
ENGINE_PARENT = module_from_spec(ESPEC)
ESPEC.loader.exec_module(ENGINE_PARENT)
ENGINE = ENGINE_PARENT.ENGINE


def configure_atlas():
    AUDIT.ATLAS.B.START = START
    AUDIT.ATLAS.B.END = END
    AUDIT.ATLAS.B.LABELS = np.arange(START, AUDIT.ATLAS.B.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return AUDIT.ATLAS.B.profile(body)


def closure_profile(case):
    return ENGINE.high_profile(case)


def render(atlas_profiles, closure_profiles):
    require(len(atlas_profiles) == comb(14, 6) == 3003, "six-body universe changed")
    candidate_rows = sum(row[7] for row in atlas_profiles)
    survivors = tuple(
        sorted(
            (survivor for row in atlas_profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    keys = tuple((row[1], row[0]) for row in survivors)
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0])
        for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    top_keys = tuple(key for key in keys if key[0] == max(first for first, _body in keys))
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate ledger changed")
    require(len(survivors) == 26 and heights == EXPECTED_HEIGHTS, "atlas survivor profile changed")
    require((len(high_keys), len(keys) - len(high_keys)) == (17, 9), "atlas route partition changed")
    require(top_keys == TOP_CASES and all(case in high_keys for case in TOP_CASES), "z1599 top shell changed")
    require(not any(1596 <= first <= 1598 for first, _body in keys), "descent gap occupied")
    atlas_profile_hash = sha256(repr(tuple(atlas_profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    require(atlas_profile_hash == EXPECTED_ATLAS_PROFILE_SHA256, "atlas profile digest changed")
    require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "atlas survivor digest changed")

    require(tuple((row[1], row[2]) for row in closure_profiles) == TOP_CASES, "closure universe changed")
    require(all(row[0] == "FORCED-HIGH" for row in closure_profiles), "closure route changed")
    for row in closure_profiles:
        case = (row[1], row[2])
        require((row[10], row[11]) == EXPECTED_GAPS[case], (case, "forced-high gaps changed"))
        require(row[10] > 0 and row[11] > 0, (case, "scalar branch survived"))
    closure_hash = sha256(repr(tuple(closure_profiles)).encode()).hexdigest()
    if EXPECTED_CLOSURE_PROFILE_SHA256 is not None:
        require(closure_hash == EXPECTED_CLOSURE_PROFILE_SHA256, "closure profile digest changed")
    minimum_gap = min(row[11] for row in closure_profiles)
    minimum_case = next((row[1], row[2]) for row in closure_profiles if row[11] == minimum_gap)

    semantic_payload = (
        EXPECTED_AUDIT_SOURCE_SHA256,
        EXPECTED_AUDIT_OUTPUT_SHA256,
        EXPECTED_ENGINE_SOURCE_SHA256,
        EXPECTED_ENGINE_OUTPUT_SHA256,
        START,
        END,
        candidate_rows,
        heights,
        high_keys,
        TOP_CASES,
        EXPECTED_GAPS,
        minimum_gap,
        1598,
        atlas_profile_hash,
        survivor_hash,
        closure_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1599 forced-high closure and 1580..1599 atlas",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"integrated_audit_source_sha256={file_sha256(AUDIT_SOURCE)}",
        f"integrated_audit_output_sha256={file_sha256(AUDIT_OUTPUT)}",
        f"forced_high_engine_source_sha256={file_sha256(ENGINE_SOURCE)}",
        f"forced_high_engine_output_sha256={file_sha256(ENGINE_OUTPUT)}",
        "inherited_cap=1599;scope=projected k=2 scalar necessary sector;all distinct later nonaligned labels;no finite label horizon",
        f"atlas_band={START}..{END};candidate_rows={candidate_rows};survivors={len(survivors)};height_counts={heights};high_tail_rows={len(high_keys)}",
        "descent_gap=1596..1598 empty;next_occupied_height=1595;z1599_rows=4",
        "closure_partition=forced_high_scalar_empty:4;empty:4;survivors:0",
    ]
    for survivor in survivors:
        body, first, delta1, suffix, upper, lower, gap, h, components, modulus, high_floor = survivor
        lines.append(
            f"ATLAS;z1={first};E={','.join(map(str, body))};"
            f"route={'HIGH-TAIL' if (first, body) in high_keys else 'EXACT'};"
            f"h={ftext(h)};r={components};L={modulus};high_floor={high_floor};"
            f"delta1={ftext(delta1)};lower={ftext(lower)};upper={ftext(upper)};"
            f"gap={ftext(gap)};suffix={suffix}"
        )
    for row in closure_profiles:
        lines.append(
            f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
            f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
            f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
            f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
            f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
            f"amplitude_sha256={row[16]};conclusion=SCALAR-EMPTY"
        )
    lines.extend((
        f"hostile_control=minimum exact one-high gap {ftext(minimum_gap)} at z1={minimum_case[0]},E={','.join(map(str,minimum_case[1]))};strictly positive",
        "consequence=projected k=2 first drift label z1<=1598",
        "scope_boundary=finite-exact THM-2995 input only;remaining 1580..1598 band, other sectors, and LRC(14) remain open",
        f"atlas_profile_sha256={atlas_profile_hash}",
        f"atlas_survivor_sha256={survivor_hash}",
        f"closure_profile_sha256={closure_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = "0"
    configure_atlas()
    bodies = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        atlas_profiles = tuple(atlas_profile(body) for body in bodies)
        closure_profiles = tuple(closure_profile(case) for case in TOP_CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=configure_atlas) as pool:
            atlas_profiles = tuple(pool.imap(atlas_profile, bodies, chunksize=4))
        with mp.get_context("spawn").Pool(args.workers) as pool:
            closure_profiles = tuple(pool.imap(closure_profile, TOP_CASES))
    payload = render(atlas_profiles, closure_profiles)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
