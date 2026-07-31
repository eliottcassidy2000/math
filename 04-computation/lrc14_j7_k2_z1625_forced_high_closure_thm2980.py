#!/usr/bin/env python3
"""Exact forced-high scalar closure of the singleton projected k=2 z1625 shell."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1631_z1630_status_descent_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1625_forced_high_closure_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "8200ba7c4bb40db669aa9b830cdbef0707b6349412e9f66ee26a51cc05018574"
EXPECTED_BASE_OUTPUT_SHA256 = "ebed197453c407e38fd3f07cafd09c4addf868866759234ee6d25afbad6486eb"
EXPECTED_PROFILE_SHA256 = "1f8244b6fca00ce97cd4fd99ce482d944e889c66f43ea634e05ba2f7b44ab72c"
EXPECTED_SEMANTIC_SHA256 = "863e9fd6fc6cd34e8658c8399f381558300c32e451ff433617e06eefc86dd07c"
CASE = (1625, (1, 2, 10, 12, 13, 14))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1630 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1630 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1629" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1630 handoff changed",
)
SPEC = spec_from_file_location("k2_z1625_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1630 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def render(row):
    require((row[1], row[2]) == CASE and row[0] == "FORCED-HIGH", "z1625 universe changed")
    require((row[5], row[6]) == (76440, 6625), "forced-high wall changed")
    require(row[10] == F(17490371, 22881893580) > 0, "two-high gap changed")
    require(row[11] == F(142815671, 572047339500) > 0, "one-high gap changed")
    require((row[12], row[13], row[14], row[15]) == (
        (2340, 1836, 1750), 490, F(239, 766948), 6708,
    ), "forced-high maximizer changed")
    profile_hash = sha256(repr(row).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASE,
        row[10], row[11], row[12:16], 1624, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    return "\n".join((
        "LRC14 projected k=2 exact z1625 forced-high closure",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1629;atlas_gap=1626..1629 empty;boundary_rows=z1625:1",
        (
            f"CASE;z1={row[1]};E={','.join(map(str,row[2]))};h={ftext(row[3])};r={row[4]};"
            f"L={row[5]};high_floor={row[6]};delta1={ftext(row[7])};lower={ftext(row[8])};"
            f"required={ftext(row[9])};two_high_gap={ftext(row[10])};"
            f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
            f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
            f"amplitude_sha256={row[16]};conclusion=SCALAR-EMPTY"
        ),
        "consequence=projected k=2 first drift label z1<=1624",
        "scope_boundary=finite-exact THM-2980 input only;full reserved theorem and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
        "",
    ))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    row = BASE.BASE.BASE.high_profile(CASE)
    payload = render(row)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
