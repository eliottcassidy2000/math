#!/usr/bin/env python3
"""Exact forced-high scalar closure of the unique projected k=2 row at z1=1768."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_z1776_ray_status_closure_thm2941.py"
PARENT_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1776_ray_status_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1768_high_wall_closure_thm2941.out"
EXPECTED_PARENT_SOURCE_SHA256 = "b8de0dd88fbc3ce0a8322d9a424fddb9bcdcd22af733f2a59e21d2a1a29cceab"
EXPECTED_PARENT_OUTPUT_SHA256 = "e9f942598572a4e8ba94ea904687b5da3c2140668b8502f3228c698b289c4bc7"
EXPECTED_PROFILE_SHA256 = "32b8ba9c891fb99495d3923a587ac76a270b64b61ae3c469ae7e2fc2469ed9e2"
EXPECTED_SEMANTIC_SHA256 = "8657e811fb1bbcae75e02c8e0990bec2a894f33a64836aab6fad3fbf6d8e01f0"

CASE = (1768, (1, 8, 10, 12, 13, 14))
HIGH_OBLIGATION = "one later label at or beyond the projected wall"
PROJECTED_WALL_RATIO = F(13, 150)
SUFFIX_SLOTS = 4
EXPECTED_GAP = F(-906233, 6582732520)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256, "parent source changed")
require(file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256, "parent output changed")
SPEC = spec_from_file_location("k2_z1768_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1776 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)
H = P.P.H
H.EXPECTED_GAPS[CASE] = EXPECTED_GAP
require(tuple(key for key in P.P.ATLAS_KEYS if key[0] == CASE[0]) == (CASE,), "z1768 row changed")
require(tuple(key for key in P.P.ATLAS_HIGH_KEYS if key[0] == CASE[0]) == (CASE,), "HIGH-TAIL dependency changed")
require(H.HIGH_WALL_RATIO == PROJECTED_WALL_RATIO, "projected-wall ratio changed")
require(H.SUFFIX_SLOTS == SUFFIX_SLOTS, "suffix arity changed")
require(max(first for first, _body in P.P.ATLAS_KEYS if first < CASE[0]) == 1758, "next height changed")


def render(row):
    require((row[0], row[1]) == CASE and row[-1] == EXPECTED_GAP < 0, "profile changed")
    profile_hash = sha256(repr(row).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, CASE, HIGH_OBLIGATION, PROJECTED_WALL_RATIO, SUFFIX_SLOTS, EXPECTED_GAP, profile_hash)
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    (first, body, h, components, L, high_floor, first_delta, lower, positive, negative, zero, ray_digest, rank4, omitted, best_high, branch, constrained, upper, gap) = row
    lines = [
        "LRC14 projected k=2 forced-high exact closure at z1=1768",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=unique z1768 HIGH-TAIL row;{HIGH_OBLIGATION};no finite label horizon",
        "forced_high_dependency=atlas HIGH-TAIL;projected wall=floor(13L/150)+1;four suffix slots",
        "ray_law=delta(r+mL)=A(r)/(r+mL);optimum=unrestricted top three plus best wall-eligible point",
        f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};high_floor={high_floor};delta1={ftext(first_delta)};lower={ftext(lower)};ray_signs=+{positive}/-{negative}/0:{zero};ray_sha256={ray_digest};unrestricted_top4={rank4};first_omitted={omitted};best_high={best_high};branch={branch};constrained={constrained};upper={ftext(upper)};gap={ftext(gap)};conclusion=SCALAR-EMPTY",
        "global_z1768_rows=1;empty=1;survivors=0",
        "atlas=next occupied height:1758",
        "consequence=projected k=2 first drift label z1<=1758",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    output = render(H.profile(CASE))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
