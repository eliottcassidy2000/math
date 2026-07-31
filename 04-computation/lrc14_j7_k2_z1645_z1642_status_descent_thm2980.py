#!/usr/bin/env python3
"""Exact projected k=2 singleton descent through z1=1645 and z1=1642.

The pinned 1600..1655 atlas has no scalar rows on 1646..1649 or 1643..1644
and has one row at each intervening boundary.  The z1645 row leaves one
literal packet after exact common status, killed by projected cells.  The
z1642 row is already empty by exact common status.  Hence z1<=1641.

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
BASE_NAME = "lrc14_j7_k2_z1650_shell_closure_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1645_z1642_status_descent_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "d38929a73bcdd36824273b727ca10a19c842f5c70e35310d1c041e8d3d23377e"
EXPECTED_BASE_OUTPUT_SHA256 = "a1a3cfe48e4e8034760c61b7fdbae282a10e052e3d0374a837788fd08e0e0b26"
EXPECTED_PROFILE_SHA256 = "cc71cf7da4b44c632138022909c14465ca3e39ba86fc944bae5ad86fc3eb2f21"
EXPECTED_SEMANTIC_SHA256 = "4526c21a7417b69ff0bf51d4b6a5b6d42b532530c4e65069b92c763cd4b02425"

BODY = (1, 4, 8, 10, 12, 14)
CASES = ((1645, BODY), (1642, BODY))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1650 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1650 output changed")
BASE_TEXT = BASE_OUTPUT.read_text(encoding="utf-8")
require(
    "height_counts=((1604, 1), (1608, 2), (1612, 19), (1616, 1), (1618, 1), (1620, 8), (1625, 1), (1630, 1), (1631, 1), (1642, 1), (1645, 1), (1650, 12))" in BASE_TEXT
    and "consequence=projected k=2 first drift label z1<=1649" in BASE_TEXT,
    "z1650 atlas handoff changed",
)
SPEC = spec_from_file_location("k2_z1645_z1642_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1650 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def profile(case):
    return BASE.ordinary_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "singleton universe changed")
    require(tuple(row[12] for row in rows) == ((5, 0, 4, 1), (16, 0, 16, 0)), "status counts changed")
    require(tuple(row[14] for row in rows) == (4, 16), "exact certificate counts changed")
    require((rows[0][16], rows[0][17]) == (1, 1), "z1645 projected packet changed")
    require(rows[0][18] == F(46807, 260897) and rows[0][21] == 1, "z1645 projected margin changed")
    require((rows[1][16], rows[1][17], rows[1][18]) == (0, 0, None), "z1642 residual changed")
    require(rows[0][13] == (
        "e2d07d9695dd56efaadccd257e9298fe981a54ef811ca263a686c9b9934fbcfa",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        "f2dd91da2a847390bdf66c6979e24b9ff7a96c9576407caf58c74f819061fc2c",
        "38aca26c0fc6c35ef414f73eddf67369bf71d87e0e225db164c5fc8d948ae4c0",
    ), "z1645 stage digests changed")
    require(rows[1][13] == (
        "e129d941ce6cd0b32e0e7200a31bd3ccae5a6bcc873e3a0354832a4aa05d2cd2",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        "a22b7b98e832f97b43b6b087dfb1ee180aa46ac786a5843c282ac10a1e458c84",
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ), "z1642 stage digests changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASES,
        tuple(row[12] for row in rows), tuple(row[14] for row in rows),
        rows[0][18], 1641, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact z1645/z1642 singleton descent",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1649;atlas_gaps=1646..1649,1643..1644 empty;boundary_rows=z1645:1,z1642:1",
        "closure_ledger=scalar:21;crude:0;status:20;status_survivors:1;exact_status_certificates:20;literal_packets:1;projected_kills:1",
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
        "consequence=projected k=2 first drift label z1<=1641",
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
