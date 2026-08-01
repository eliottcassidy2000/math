#!/usr/bin/env python3
"""Exact projected k=2 singleton descent through z1=1618 and z1=1616.

The pinned 1600..1655 atlas has no scalar row at 1619 or 1617 and has one
row at each intervening boundary.  The z1618 row leaves three literal packets
after exact common status, all killed by projected cells.  The z1616 row is
already empty by exact common status.  Hence z1<=1615.

This remains a finite-exact boundary input for reserved THM-2980 only.
"""

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
BASE_NAME = "lrc14_j7_k2_z1620_mixed_shell_closure_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1618_z1616_status_descent_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "eef6a5c0e7e988783357460c8f4144f3ce3a1a4e0cd7a52b949c6575661ce1e8"
EXPECTED_BASE_OUTPUT_SHA256 = "bb4e5df57439de20b40b72481257a20e9ebacd278d42a3a25b947b8d8f39504b"
EXPECTED_PROFILE_SHA256 = "555fe793a7505a96bc6efa4fbb6798e405d3faa83f564c6f0608600d7af80793"
EXPECTED_SEMANTIC_SHA256 = "889fbb401744fc4ec809e2de5e283b708610053b31824986833ab1b8952e7e67"

CASES = (
    (1618, (2, 6, 8, 10, 12, 14)),
    (1616, (1, 4, 8, 10, 12, 14)),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1620 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1620 output changed")
BASE_TEXT = BASE_OUTPUT.read_text(encoding="utf-8")
require(
    "consequence=projected k=2 first drift label z1<=1619" in BASE_TEXT
    and "height_counts=((1604, 1), (1608, 2), (1612, 19), (1616, 1), (1618, 1), (1620, 8), (1625, 1), (1630, 1), (1631, 1), (1642, 1), (1645, 1), (1650, 12))"
    in (ROOT / "05-knowledge/results/lrc14_j7_k2_z1650_shell_closure_thm2980.out").read_text(encoding="utf-8"),
    "z1620/atlas handoff changed",
)
SPEC = spec_from_file_location("k2_z1618_z1616_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1620 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE


def profile(case):
    return ENGINE.ordinary_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "singleton universe changed")
    require(tuple(row[12] for row in rows) == ((3, 0, 0, 3), (5, 0, 5, 0)), "status counts changed")
    require(tuple(row[14] for row in rows) == (0, 5), "exact certificate counts changed")
    require((rows[0][16], rows[0][17]) == (3, 3), "z1618 projected packet changed")
    require(
        rows[0][18] == F(1296, 5551) and rows[0][19] == 5 and rows[0][21] == 1,
        "z1618 projected margin changed",
    )
    require((rows[1][16], rows[1][17], rows[1][18]) == (0, 0, None), "z1616 residual changed")
    require(rows[0][13] == (
        "1e6e51507d7ce4bc7f29db2238ef3366896cf3c6b739a22908e0f6a45205ce68",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        "1d826a3626f89c0213d036c0a4e92b701e40a27ee73dadd46e61a03b115d615d",
    ), "z1618 stage digests changed")
    require(rows[1][13] == (
        "7e71ae3d9ae39fe728ee893dfc1c882d3670f0ba0b56a54586d1e608b976546f",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        "4a054b713a2122541924100dfa12c76a482adac74db2948a52e8448c2eac9d96",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ), "z1616 stage digests changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASES,
        tuple(row[12] for row in rows), tuple(row[14] for row in rows),
        rows[0][18], 1615, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact z1618/z1616 singleton descent",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1619;atlas_gaps=1619,1617 empty;boundary_rows=z1618:1,z1616:1",
        "closure_ledger=scalar:8;crude:0;status:5;status_survivors:3;exact_status_certificates:5;literal_packets:3;projected_kills:3",
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
        "consequence=projected k=2 first drift label z1<=1615",
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
