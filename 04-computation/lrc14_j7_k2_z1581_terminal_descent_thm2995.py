#!/usr/bin/env python3
"""Terminal projected k=2 descent at z1=1581 for the THM-2995 band.

The pinned atlas has three rows at z1=1581 and no rows at 1580 or 1582..1585.
The ordinary row leaves one literal packet, which is projected-empty.  Both
HIGH-TAIL rows are scalar-empty by exact one-high ray maxima.  Hence every
projected k=2 scalar row in 1580..1599 is empty and the cap is z1<=1579,
pending the integrated 60,060-row promotion audit.
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
BASE_NAME = "lrc14_j7_k2_z1588_z1586_mixed_section_descent_thm2995"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1581_terminal_descent_thm2995.out"

EXPECTED_BASE_SOURCE_SHA256 = "036f42af2e6d8ec40cc1cb233c2f2d0bdfd191777ba0e08048b4e62926c67dd4"
EXPECTED_BASE_OUTPUT_SHA256 = "22c96c5d09087d2dab941346eceb278cc3408923802b6e39b8c446f1445e235d"
EXPECTED_PROFILE_SHA256 = "9619c193ea40a1fa61489cf46dbe6eaa8e5cd54505c9cfccdf1867a9ef9e4c8d"
EXPECTED_SEMANTIC_SHA256 = "b56b5d091b0e6083061be59a2b4dc2b8f11307b7152da89e29660d4d00db5f26"

STANDARD_BODY = (1, 4, 8, 10, 12, 14)
CASES = (
    (1581, STANDARD_BODY),
    (1581, (1, 4, 10, 12, 13, 14)),
    (1581, (1, 8, 10, 12, 13, 14)),
)
EXPECTED_HIGH_GAPS = {
    CASES[1]: (F(543763, 635608155), F(3458937, 13135901870)),
    CASES[2]: (F(69113, 80863965), F(163705343, 918291186540)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1588/z1586 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1588/z1586 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1585" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1588/z1586 handoff changed",
)
SPEC = spec_from_file_location("k2_z1581_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1588/z1586 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE
require(all(case in BASE.ATLAS_KEYS for case in CASES), "z1581 atlas rows changed")
require(not any(first == 1580 or 1582 <= first <= 1585 for first, _body in BASE.ATLAS_KEYS), "terminal atlas gap changed")


def profile(case):
    return ENGINE.ordinary_profile(case) if case[1] == STANDARD_BODY else ENGINE.high_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "z1581 universe changed")
    ordinary = rows[0]
    high = rows[1:]
    require(ordinary[0] == "ALL-LABEL" and all(row[0] == "FORCED-HIGH" for row in high), "route split changed")
    require(ordinary[12] == (1, 0, 0, 1) and ordinary[14] == 0, "ordinary status ledger changed")
    require(ordinary[16] == ordinary[17] == 1, "ordinary packet ledger changed")
    require((ordinary[18], ordinary[19]) == (F(7482, 47957), 1), "ordinary projected control changed")
    for row in high:
        case = (row[1], row[2])
        require((row[10], row[11]) == EXPECTED_HIGH_GAPS[case], (case, "forced-high gaps changed"))
        require(row[10] > 0 and row[11] > 0, (case, "forced-high branch survived"))

    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256,
        EXPECTED_BASE_OUTPUT_SHA256,
        CASES,
        ordinary[12],
        ordinary[14],
        ordinary[16],
        ordinary[18],
        EXPECTED_HIGH_GAPS,
        1579,
        profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1581 terminal descent for the 1580..1599 band",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1585;atlas_gaps=1582..1585,1580 empty;boundary_rows=z1581:3",
        "closure_partition=all_label_projected:1;forced_high_scalar:2;empty:3;survivors:0",
        "ordinary_ledger=scalar:1;crude:0;status:0;status_survivors:1;exact_status_certificates:0;literal_packets:1;projected_kills:1",
    ]
    for row in rows:
        if row[0] == "ALL-LABEL":
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
        else:
            lines.append(
                f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
                f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
                f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
                f"amplitude_sha256={row[16]};conclusion=SCALAR-EMPTY"
            )
    lines.extend((
        "hostile_control=the ordinary row is not status-empty and requires its unique projected packet;both HIGH-TAIL rows have strictly positive exact one-high gaps",
        "band_consequence=all projected k=2 scalar rows with 1580<=z1<=1599 are empty",
        "consequence=projected k=2 first drift label z1<=1579",
        "scope_boundary=finite-exact intended THM-2995 band input complete;THM-2995 remains reserved pending integrated promotion audit;other sectors and LRC(14) remain open",
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
