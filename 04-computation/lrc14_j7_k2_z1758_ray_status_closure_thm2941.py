#!/usr/bin/env python3
"""Uniform ray/status closure of the unique projected k=2 row at z1=1758."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_z1768_high_wall_closure_thm2941.py"
PARENT_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1768_high_wall_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1758_ray_status_closure_thm2941.out"
EXPECTED_PARENT_SOURCE_SHA256 = "96fed14134bc57267b74636ef3174499c8f4bdb1e66806e865df5a02d640b5b8"
EXPECTED_PARENT_OUTPUT_SHA256 = "e8a5c6fe160680d2191c05110ae5f6d15f77ecbbf612bad86e1dddb9950a2b63"
EXPECTED_PROFILE_SHA256 = "670a20fc93b692a70d08b6bf02a31c63beb4030684286cb9a3ca0026a43e1221"
EXPECTED_SEMANTIC_SHA256 = "a84561137936a327b5a7c3767854327983969ad6a82e694be41b017be97fd646"

CASE = (1758, (1, 4, 8, 10, 12, 14))
QUANTIFIER = "distinct later nonaligned labels"
EXPECTED_COUNTS = (2, 0, 2, 0)
EXPECTED_STAGE_SHA256 = (
    "e78b52de33839ff3bb4fc51d6b2b477886877dbdca09f99b178e9ce4659b4741",
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    "41b3cdfdd6f1e191a9ff5c7fb359b4b0a34e3915b9492f4dea258a81523ee615",
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
)
EXPECTED_RAY_SHA256 = "276909b5f679a0675eb6c79a12d9cf2b5024f1d3528cc0a7df7014a30f607674"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256, "parent source changed")
require(file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256, "parent output changed")
SPEC = spec_from_file_location("k2_z1758_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1768 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)
E = P.P.E
ATLAS_KEYS = P.P.P.ATLAS_KEYS
require(tuple(key for key in ATLAS_KEYS if key[0] == CASE[0]) == (CASE,), "z1758 row changed")
require(max(first for first, _body in ATLAS_KEYS if first < CASE[0]) == 1750, "next height changed")


def profile():
    first, body = CASE
    carrier = E.U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), E.U.suffix.A.RULER)
    lower = h * E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = E.ray_and_status(first, body, carrier, h, lower, L)
    (_amplitudes, ray_digest, divisor_count, trials, first_delta, first_d, scalar, crude, status, states, stages) = result
    counts = (len(scalar), len(crude), len(status), len(states))
    require(counts == EXPECTED_COUNTS, ("counts changed", counts))
    require(stages == EXPECTED_STAGE_SHA256, ("stage digests changed", stages))
    require(ray_digest == EXPECTED_RAY_SHA256, "ray digest changed")
    require(not states, "a common-status state survived")
    return (first, body, h, len(carrier), L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages)


def render(row):
    profile_hash = sha256(repr(row).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, ("profile digest changed", profile_hash))
    semantic_payload = (EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, CASE, QUANTIFIER, EXPECTED_COUNTS, EXPECTED_STAGE_SHA256, profile_hash)
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic_hash))
    (first, body, h, components, L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages) = row
    lines = [
        "LRC14 projected k=2 exact ray/status closure at z1=1758",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=unique z1758 atlas row;all {QUANTIFIER};no finite label horizon",
        "pipeline=residue rays;denominator multisets;all-divisor crude;exact common-K5 Hunter status",
        f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};ray_sha256={ray_digest};denominators={divisor_count};trials={trials};counts={counts};stage_sha256={stages};conclusion=STATUS-EMPTY",
        "global_z1758_rows=1;empty=1;survivors=0",
        "atlas=next occupied height:1750",
        "consequence=projected k=2 first drift label z1<=1750",
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
    output = render(profile())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
