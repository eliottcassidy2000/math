#!/usr/bin/env python3
"""Close the unique projected k=2 scalar-atlas row at z1=1660 exactly."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, lcm
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1668_composite_closure_thm2941"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
if not BASE_SOURCE.exists():
    BASE_SOURCE = HERE.parent / f"{BASE_NAME}.py"
if not BASE_OUTPUT.exists():
    BASE_OUTPUT = HERE.parent / f"{BASE_NAME}.out"
OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k2_z1660_status_descent_closure_thm2941.out"
    if HERE.parent.name == "04-computation"
    else HERE.with_suffix(".out")
)

EXPECTED_BASE_SOURCE_SHA256 = "4719dea2b29b0b1da645bd9902fbd883b602b6fad4f1d63300aedcfb30198764"
EXPECTED_BASE_OUTPUT_SHA256 = "3b3a5d577e41afd4b9a2c9748c74dfb23db70b6723f1101e66d5fae939adff7a"
EXPECTED_ATLAS_PROFILE_SHA256 = "3d7a1a61c5278876444296cd52c059c69a1d58bf6beeb0abd18e05c2b4641120"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "f36b8b517d97b5552de4ed791006a75184c786ac6c09186f5102968bbcf439fc"
EXPECTED_CLOSURE_PROFILE_SHA256 = "79ac8bab2f4b3b6412231814589864111fe91f3a1893565338b9b061b701ae7a"
EXPECTED_SEMANTIC_SHA256 = "3d0ee0388a8ea9ef5f49b923423a75410c7a51a9c6f48a4384bf49b4c3a689ee"

START = 1656
END = 1660
EXPECTED_CANDIDATE_ROWS = 15_015
CLOSURE_CASE = (1660, (1, 4, 8, 10, 12, 14))
EXPECTED_HEIGHTS = ((1656, 5), (1660, 1))
EXPECTED_KEYS = tuple(
    sorted(
        (
            (1656, (1, 2, 8, 10, 12, 14)),
            (1656, (1, 4, 8, 10, 12, 14)),
            (1656, (1, 6, 8, 10, 12, 14)),
            (1656, (1, 8, 10, 12, 13, 14)),
            (1656, (2, 4, 8, 10, 12, 14)),
            CLOSURE_CASE,
        )
    )
)
EXPECTED_HIGH_KEYS = ((1656, (1, 8, 10, 12, 13, 14)),)
EXPECTED_STAGE_DIGESTS = (
    "a8b36d9993c1b049a6f283e637aa05685e8bc176491718a9203c842dcce4385f",
    sha256(b"()").hexdigest(),
    "902bc4e219c6dbeb1bfae07eda99eb2a29985143303ea004b7c688ad738284d2",
    sha256(b"()").hexdigest(),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1668/z1665 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1668/z1665 output changed")
SPEC = spec_from_file_location("k2_z1660_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1668/z1665 base")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE


def configure_atlas():
    BASE.BASE.ATLAS.B.START = START
    BASE.BASE.ATLAS.B.END = END
    BASE.BASE.ATLAS.B.LABELS = np.arange(START, BASE.BASE.ATLAS.B.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return BASE.BASE.ATLAS.B.profile(body)


def closure_profile():
    first, body = CLOSURE_CASE
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, L)
    require((len(result[6]), len(result[7]), len(result[8]), len(result[9])) == (1, 0, 1, 0), "closure counts changed")
    require(result[10] == EXPECTED_STAGE_DIGESTS, "stage digests changed")
    certificate_count = BASE.BASE.exact_status_audit(CLOSURE_CASE, L, result[8], result[10][2])
    require(certificate_count == 1, "exact certificate count changed")
    return (
        first, body, h, len(carrier), L, lower, result[4], result[5], result[1],
        result[2], result[3], result[10], certificate_count,
    )


def render(atlas_profiles, closure):
    require(len(atlas_profiles) == comb(14, 6) == 3003, "atlas universe changed")
    candidate_rows = sum(row[7] for row in atlas_profiles)
    survivors = tuple(
        sorted((survivor for row in atlas_profiles for survivor in row[10]), key=lambda row: (row[1], row[0]))
    )
    keys = tuple((row[1], row[0]) for row in survivors)
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0]) for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate ledger changed")
    require(keys == EXPECTED_KEYS, ("atlas keys changed", keys))
    require(heights == EXPECTED_HEIGHTS, ("height ledger changed", heights))
    require(high_keys == EXPECTED_HIGH_KEYS, ("HIGH partition changed", high_keys))
    require(not any(1657 <= first <= 1659 for first, _body in keys), "descent gap occupied")
    atlas_hash = sha256(repr(tuple(atlas_profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    closure_hash = sha256(repr(closure).encode()).hexdigest()
    if EXPECTED_ATLAS_PROFILE_SHA256 is not None:
        require(atlas_hash == EXPECTED_ATLAS_PROFILE_SHA256, "atlas profile digest changed")
    if EXPECTED_ATLAS_SURVIVOR_SHA256 is not None:
        require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "atlas survivor digest changed")
    if EXPECTED_CLOSURE_PROFILE_SHA256 is not None:
        require(closure_hash == EXPECTED_CLOSURE_PROFILE_SHA256, "closure profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, START, END,
        candidate_rows, heights, high_keys, CLOSURE_CASE, EXPECTED_STAGE_DIGESTS,
        1656, atlas_hash, survivor_hash, closure_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact common-status descent at z1=1660",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1660;scope=unique z1660 scalar-atlas row;all distinct later nonaligned labels;no finite label horizon",
        f"atlas_band=1656..1660;candidate_rows={candidate_rows};height_counts={heights};empty_firsts=1657..1659;next_occupied_height=1656",
        "closure_counts=scalar:1;crude:0;status:1;status_survivors:0;literal_packets:0;survivors:0",
        "status_replay=one infeasibility certificate verified exactly over Q;solver-selected dual basis excluded from stage and semantic digests",
        (
            f"CASE;z1={closure[0]};E={','.join(map(str, closure[1]))};h={ftext(closure[2])};"
            f"r={closure[3]};L={closure[4]};lower={ftext(closure[5])};"
            f"delta1={ftext(closure[6])};first_d={closure[7]};ray_sha256={closure[8]};"
            f"divisors={closure[9]};trials={closure[10]};stage_sha256={closure[11]};"
            f"exact_status_certificates={closure[12]};conclusion=EMPTY"
        ),
        "consequence=projected k=2 first drift label z1<=1656",
        f"atlas_profile_sha256={atlas_hash}",
        f"atlas_survivor_sha256={survivor_hash}",
        f"closure_profile_sha256={closure_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
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
        atlas_profiles = [atlas_profile(body) for body in bodies]
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=configure_atlas) as pool:
            atlas_profiles = list(pool.imap(atlas_profile, bodies, chunksize=4))
    atlas_profiles.sort(key=lambda row: row[0])
    payload = render(tuple(atlas_profiles), closure_profile())
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
