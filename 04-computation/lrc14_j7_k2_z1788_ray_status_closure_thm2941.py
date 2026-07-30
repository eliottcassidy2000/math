#!/usr/bin/env python3
"""Uniform ray/status closure of the unique projected k=2 row at z1=1788.

After the z1=1790 closure, the exact all-body atlas has no row at 1789 and
has exactly one row at 1788, with body (1,4,8,10,12,14).  This verifier
enumerates every scalar-eligible denominator multiset for four distinct later
nonaligned labels.  All 39 states fail either an all-divisor capacity test or
the exact common-K5 Hunter-status feasibility problem, so no literal packet
or label horizon is needed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z1790_exact_descent_closure_thm2941.py"
)
PARENT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1788_ray_status_closure_thm2941.out"
)
EXPECTED_PARENT_SOURCE_SHA256 = (
    "c80c152300e60a2830d7cff4af1ac104d0f7ca1dfea978033dbaf03a11809d5b"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "b03b46a6c438773ef1c433c435828a5426e18c998999cbee543541904b85f20b"
)
EXPECTED_PROFILE_SHA256 = (
    "faf85ab98de47cf4a32c9fc019356739b81c95d2712ef10ddd30aab63f40dbe0"
)
EXPECTED_SEMANTIC_SHA256 = (
    "16bef50e606fab62abb7e3d9737be4ebc93f466ae5847c570dd0e975631556a1"
)

CASE = (1788, (1, 4, 8, 10, 12, 14))
QUANTIFIER = "distinct later nonaligned labels"
EXPECTED_COUNTS = (39, 19, 20, 0)
EXPECTED_STAGE_SHA256 = (
    "2715b8e5c1d0c60676548d5890bfed0cf9e74585b4595cfe8c523ab5af6c470f",
    "b7f971d2ccafdb7329e6d2941025c7edc80a480ce588a5cb31aed4ca133a2013",
    "400c7e0575eda1edc440d6c8e64f396a026ac3806ee5e7c8d536d98d76ca0ccb",
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
)
EXPECTED_RAY_SHA256 = (
    "276909b5f679a0675eb6c79a12d9cf2b5024f1d3528cc0a7df7014a30f607674"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256,
    "z1790 parent source changed",
)
require(
    file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256,
    "z1790 parent output changed",
)
SPEC = spec_from_file_location("k2_z1788_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1790 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)

require(
    tuple(key for key in P.ATLAS_KEYS if key[0] == CASE[0]) == (CASE,),
    "z1788 is not the unique atlas row",
)
require(not any(first == 1789 for first, _body in P.ATLAS_KEYS), "z1789 occupied")
require(
    max(first for first, _body in P.ATLAS_KEYS if first < CASE[0]) == 1784,
    "next occupied atlas height changed",
)


def profile():
    first, body = CASE
    carrier = P.E.U.suffix.A.carrier_for(body)
    require(P.E.P.A.carrier_for(body) == carrier, "carrier engines disagree")
    h = F(sum(right - left for left, right in carrier), P.E.U.suffix.A.RULER)
    lower = h * P.E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(P.E.P.A.RULER % L == 0, "body ruler left master ruler")
    (
        _amplitudes,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = P.E.ray_and_status(first, body, carrier, h, lower, L)
    counts = (len(scalar), len(crude_kills), len(status_kills), len(states))
    require(counts == EXPECTED_COUNTS, ("stage counts changed", counts))
    require(stage_digests == EXPECTED_STAGE_SHA256, "stage digests changed")
    require(ray_digest == EXPECTED_RAY_SHA256, "ray digest changed")
    require(not states, "a common-status state survived")
    return (
        first,
        body,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
    )


def render(row):
    profile_hash = sha256(repr(row).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_PARENT_SOURCE_SHA256,
        EXPECTED_PARENT_OUTPUT_SHA256,
        CASE,
        QUANTIFIER,
        EXPECTED_COUNTS,
        EXPECTED_STAGE_SHA256,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    (
        first,
        body,
        h,
        components,
        L,
        lower,
        first_delta,
        first_d,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
    ) = row
    lines = [
        "LRC14 projected k=2 exact ray/status closure at z1=1788",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=unique z1788 atlas row;all {QUANTIFIER};no finite label horizon",
        (
            "pipeline=residue rays;denominator multisets;all-divisor crude;"
            "exact common-K5 Hunter status"
        ),
        (
            f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};"
            f"r={components};L={L};lower={ftext(lower)};delta1={ftext(first_delta)};"
            f"first_d={first_d};ray_sha256={ray_digest};"
            f"denominators={divisor_count};trials={trials};counts={counts};"
            f"stage_sha256={stage_digests};conclusion=STATUS-EMPTY"
        ),
        "global_z1788_rows=1;empty=1;survivors=0",
        "atlas=empty:1789..1790 after parent closure;next occupied height:1784",
        "consequence=projected k=2 first drift label z1<=1784",
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
    row = profile()
    output = render(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
