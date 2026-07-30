#!/usr/bin/env python3
"""Uniform ray/status closure of the unique projected k=2 row at z1=1776."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_z1780_hybrid_closure_thm2941.py"
PARENT_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1780_hybrid_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1776_ray_status_closure_thm2941.out"
EXPECTED_PARENT_SOURCE_SHA256 = "be378cd295799b1069793c86557d09eb8f682bf3c8a890c8fedcc0ff42ee100b"
EXPECTED_PARENT_OUTPUT_SHA256 = "afa96047c84a544b77142e34dfd25ad3923bace4d54f3f5ce47b0910f99ee27b"
EXPECTED_PROFILE_SHA256 = "245ae162ab60c6fafb5276b02d500eb59fb961c30bb597392102eb5c6b9bef6f"
EXPECTED_SEMANTIC_SHA256 = "0fc768d027fffc1ea690e266ca44a38f47c6545384e9ba2fb7037a1c1f7d29a1"

CASE = (1776, (1, 4, 8, 10, 12, 14))
QUANTIFIER = "distinct later nonaligned labels"
EXPECTED_COUNTS = (77, 29, 48, 0)
EXPECTED_STAGE_SHA256 = (
    "12373ed0c53612af7b54c8e7a054448c2b6162cf593f25d124e4b68e43df8dc6",
    "bccb33564b4579942021106d5ce2076f89e9685022ef4789db2ce41a06342dc4",
    "5d550387aedab5fddeeba12d8afc66d142d37e0ec79667c26f9ad7d5d922ca09",
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
SPEC = spec_from_file_location("k2_z1776_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1780 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)
E = P.E
require(tuple(key for key in P.ATLAS_KEYS if key[0] == CASE[0]) == (CASE,), "z1776 row changed")
require(max(first for first, _body in P.ATLAS_KEYS if first < CASE[0]) == 1768, "next height changed")


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
    require(stages == EXPECTED_STAGE_SHA256, "stage digests changed")
    require(ray_digest == EXPECTED_RAY_SHA256, "ray digest changed")
    require(not states, "a common-status state survived")
    return (first, body, h, len(carrier), L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages)


def render(row):
    profile_hash = sha256(repr(row).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, CASE, QUANTIFIER, EXPECTED_COUNTS, EXPECTED_STAGE_SHA256, profile_hash)
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    (first, body, h, components, L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages) = row
    lines = [
        "LRC14 projected k=2 exact ray/status closure at z1=1776",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=unique z1776 atlas row;all {QUANTIFIER};no finite label horizon",
        "pipeline=residue rays;denominator multisets;all-divisor crude;exact common-K5 Hunter status",
        f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};ray_sha256={ray_digest};denominators={divisor_count};trials={trials};counts={counts};stage_sha256={stages};conclusion=STATUS-EMPTY",
        "global_z1776_rows=1;empty=1;survivors=0",
        "atlas=next occupied height:1768",
        "consequence=projected k=2 first drift label z1<=1768",
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
