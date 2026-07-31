#!/usr/bin/env python3
"""Exact projected k=2 singleton descent through z1=1631 and z1=1630."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1645_z1642_status_descent_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1631_z1630_status_descent_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "f1b4ee77b3749b52be5e67a3b733282b2f1ccb8c70197240af5a820c7a3786c8"
EXPECTED_BASE_OUTPUT_SHA256 = "e98282c6266f474378ddeb725b18f7d2542b82e58ca5b6a6512f8a7600a39148"
EXPECTED_PROFILE_SHA256 = "ac4dca9bc74e2a860c39a0be8cf402c131f47868d189fa7915727c0193aa5367"
EXPECTED_SEMANTIC_SHA256 = "455f5935c9e83011ba1e9eae34708e92ec4114ad9e919a198b0a572f19bc92ca"

BODY = (1, 4, 8, 10, 12, 14)
CASES = ((1631, BODY), (1630, BODY))
EMPTY = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1642 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1642 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1641" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1642 handoff changed",
)
SPEC = spec_from_file_location("k2_z1631_z1630_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1642 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def profile(case):
    return BASE.BASE.ordinary_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "singleton universe changed")
    require(tuple(row[12] for row in rows) == ((1, 0, 0, 1), (1, 0, 1, 0)), "status counts changed")
    require(tuple(row[14] for row in rows) == (0, 1), "certificate counts changed")
    require((rows[0][16], rows[0][17], rows[0][18], rows[0][21]) == (1, 1, F(10555, 21203), 1), "z1631 projected closure changed")
    require((rows[1][16], rows[1][17], rows[1][18]) == (0, 0, None), "z1630 residual changed")
    require(rows[0][13] == (
        "20821d5f0aca13b9c49ead7fa072d5cab9903e651cd3c32ac1083653732907dd",
        EMPTY, EMPTY,
        "fcfe0ac8275fafde0a745063210af4348e303af9647560962eaa8ea60b59d24a",
    ), "z1631 stage digests changed")
    require(rows[1][13] == (
        "c79449b4a9e2f0f63180d71a556d4a678981af46307727bb6692f46ca98b1cd9",
        EMPTY,
        "ceef3adde664e209d46dc760c8d6196b429bc501401428d366e7e66b58f1df3d",
        EMPTY,
    ), "z1630 stage digests changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASES,
        tuple(row[12] for row in rows), tuple(row[14] for row in rows),
        rows[0][18], 1629, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact z1631/z1630 singleton descent",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1641;atlas_gap=1632..1641 empty;boundary_rows=z1631:1,z1630:1",
        "closure_ledger=scalar:2;crude:0;status:1;status_survivors:1;exact_status_certificates:1;literal_packets:1;projected_kills:1",
    ]
    for row in rows:
        lines.append(
            f"CASE;z1={row[1]};E={','.join(map(str,row[2]))};h={ftext(row[3])};r={row[4]};"
            f"L={row[5]};lower={ftext(row[6])};delta1={ftext(row[7])};first_d={row[8]};"
            f"ray_sha256={row[9]};divisors={row[10]};trials={row[11]};counts={row[12]};"
            f"stage_sha256={row[13]};exact_status_certificates={row[14]};cells={row[15]};"
            f"packets={row[16]};kills={row[17]};minimum_margin={ftext(row[18])};"
            f"maximum_prefix_cells={row[19]};direct_mass={ftext(row[21])};"
            f"state_sha256={row[22]};conclusion=EMPTY"
        )
    lines.extend((
        "consequence=projected k=2 first drift label z1<=1629",
        "scope_boundary=finite-exact THM-2980 input only;full reserved theorem and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(2, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = "0"
    if args.workers == 1:
        rows = tuple(profile(case) for case in CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = tuple(pool.imap(profile, CASES))
    payload = render(rows)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
