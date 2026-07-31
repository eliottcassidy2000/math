#!/usr/bin/env python3
"""Exact projected k=2 status/projected descent at z1=1590.

The pinned atlas has three ordinary rows at z1=1590 and no rows at
1591..1593.  One row closes entirely by exact common-status certificates;
the other two leave 34 status survivors, which expand to 51 literal packets.
Every packet is projected-empty, so the cap drops to 1589.
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
BASE_NAME = "lrc14_j7_k2_z1595_z1594_mixed_descent_thm2995"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1590_status_projected_descent_thm2995.out"

EXPECTED_BASE_SOURCE_SHA256 = "53fdc5c3c146d4b66a2656d249ea2c8ddec6721342c69063b49c41976e2c021f"
EXPECTED_BASE_OUTPUT_SHA256 = "7afd64dcbb294fca43ff8d254fcf72389592f8c414ce4f7f49319f1d4353dc7a"
EXPECTED_PROFILE_SHA256 = "c7785733f36acb06e7ad7f596407414980ae13dc3f89f5a65f828804609044a4"
EXPECTED_SEMANTIC_SHA256 = "e7e7f3b4a9789e74e14593f4d800015f4688882e7ad00b675948491455793c25"

CASES = tuple((1590, body) for body in (
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
))
EXPECTED_ROWS = {
    CASES[0]: ((15, 0, 15, 0), 15, 0, None, 0),
    CASES[1]: ((26, 4, 13, 9), 13, 17, F(1026, 16471), 2),
    CASES[2]: ((29, 2, 2, 25), 2, 34, F(1296, 5551), 5),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1595/z1594 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1595/z1594 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1593" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1595/z1594 handoff changed",
)
SPEC = spec_from_file_location("k2_z1590_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1595/z1594 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE
require(all(case in BASE.ATLAS_KEYS for case in CASES), "z1590 atlas rows changed")
require(not any(1591 <= first <= 1593 for first, _body in BASE.ATLAS_KEYS), "z1591..1593 gap changed")


def profile(case):
    return ENGINE.ordinary_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "z1590 universe changed")
    require(all(row[0] == "ALL-LABEL" for row in rows), "z1590 route changed")
    for row in rows:
        case = (row[1], row[2])
        counts, certificates, packets, margin, prefix = EXPECTED_ROWS[case]
        require(row[12] == counts and row[14] == certificates, (case, "status ledger changed"))
        require(row[16] == row[17] == packets, (case, "packet ledger changed"))
        require((row[18], row[19]) == (margin, prefix), (case, "projected control changed"))

    totals = tuple(sum(row[12][index] for row in rows) for index in range(4))
    require(totals == (70, 6, 30, 34), "aggregate status ledger changed")
    require(sum(row[14] for row in rows) == 30, "certificate count changed")
    require(sum(row[16] for row in rows) == sum(row[17] for row in rows) == 51, "packet count changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256,
        EXPECTED_BASE_OUTPUT_SHA256,
        CASES,
        EXPECTED_ROWS,
        totals,
        30,
        51,
        1589,
        profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact status/projected descent at z1590",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1593;atlas_gap=1591..1593 empty;boundary_rows=z1590:3",
        "closure_partition=all_label_status_projected:3;empty:3;survivors:0",
        "closure_ledger=scalar:70;crude:6;status:30;status_survivors:34;exact_status_certificates:30;literal_packets:51;projected_kills:51",
    ]
    for row in rows:
        lines.append(
            f"CASE;route=ALL-LABEL;z1={row[1]};E={','.join(map(str,row[2]))};"
            f"h={ftext(row[3])};r={row[4]};L={row[5]};lower={ftext(row[6])};"
            f"delta1={ftext(row[7])};first_d={row[8]};ray_sha256={row[9]};"
            f"divisors={row[10]};trials={row[11]};counts={row[12]};"
            f"stage_sha256={row[13]};exact_status_certificates={row[14]};"
            f"cells={row[15]};packets={row[16]};kills={row[17]};"
            f"minimum_margin={ftext(row[18])};maximum_prefix_cells={row[19]};"
            f"direct_mass={ftext(row[21])};state_sha256={row[22]};conclusion=EMPTY"
        )
    lines.extend((
        "hostile_control=34 status survivors expand rather than being identified with packets;the exact packet count is 51 and the minimum projected margin is 1026/16471",
        "consequence=projected k=2 first drift label z1<=1589",
        "scope_boundary=finite-exact THM-2995 input only;remaining 1580..1589 band, other sectors, and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(3, mp.cpu_count() or 1))
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
