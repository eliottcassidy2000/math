#!/usr/bin/env python3
"""Terminal finite-exact projected k=2 descent through z1608 and z1604.

The pinned atlas has two rows at z1608, one at z1604, and no rows at
1609..1611, 1605..1607, or 1600..1603.  Exact common status closes the two
ordinary rows, while the remaining high-tail row has a positive exact
one-high gap.  Thus every projected k=2 scalar row in 1600..1679 is empty
and the finite-exact cap is z1<=1599.

This completes the intended band input but does not itself promote reserved
THM-2980; an integrated theorem audit is still required.
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
BASE_NAME = "lrc14_j7_k2_z1612_shell_closure_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
ATLAS_NAME = "lrc14_j7_k2_z1650_shell_closure_thm2980"
ATLAS_SOURCE = ROOT / f"04-computation/{ATLAS_NAME}.py"
ATLAS_OUTPUT = ROOT / f"05-knowledge/results/{ATLAS_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1608_z1604_terminal_descent_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "8861a4d02bab7a49c8f19a5b5c118b914e57ce504cfafe5926b6a647574fedb3"
EXPECTED_BASE_OUTPUT_SHA256 = "77c8e4d359c0649143bdca3a6cfcfedb69625251f4eeaaa31599ce2d5d6fa99b"
EXPECTED_ATLAS_SOURCE_SHA256 = "d38929a73bcdd36824273b727ca10a19c842f5c70e35310d1c041e8d3d23377e"
EXPECTED_ATLAS_OUTPUT_SHA256 = "a1a3cfe48e4e8034760c61b7fdbae282a10e052e3d0374a837788fd08e0e0b26"
EXPECTED_PROFILE_SHA256 = "be809cc1d19c45e0f451fd77404666ce397e79a55c41dbf87fd451efe48026c9"
EXPECTED_SEMANTIC_SHA256 = "42238eae59f0df241ca68919b5987444f482f7a009eda6e813b23991a3553796"

STANDARD_BODY = (1, 4, 8, 10, 12, 14)
HIGH_BODY = (1, 4, 10, 12, 13, 14)
CASES = (
    (1608, STANDARD_BODY),
    (1608, HIGH_BODY),
    (1604, STANDARD_BODY),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1612 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1612 output changed")
require(file_sha256(ATLAS_SOURCE) == EXPECTED_ATLAS_SOURCE_SHA256, "1600..1655 atlas source changed")
require(file_sha256(ATLAS_OUTPUT) == EXPECTED_ATLAS_OUTPUT_SHA256, "1600..1655 atlas output changed")
require(
    "consequence=projected k=2 first drift label z1<=1611"
    in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1612 handoff changed",
)
ATLAS_TEXT = ATLAS_OUTPUT.read_text(encoding="utf-8")
require(
    "atlas_band=1600..1655;candidate_rows=168168;survivors=49" in ATLAS_TEXT
    and "height_counts=((1604, 1), (1608, 2), (1612, 19), (1616, 1), (1618, 1), (1620, 8), (1625, 1), (1630, 1), (1631, 1), (1642, 1), (1645, 1), (1650, 12))"
    in ATLAS_TEXT,
    "terminal atlas ledger changed",
)
SPEC = spec_from_file_location("k2_z1608_z1604_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1612 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE


def profile(case):
    return ENGINE.high_profile(case) if case[1] == HIGH_BODY else ENGINE.ordinary_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "terminal universe changed")
    ordinary = (rows[0], rows[2])
    high = rows[1]
    require(tuple(row[0] for row in ordinary) == ("ALL-LABEL", "ALL-LABEL"), "ordinary typing changed")
    require(high[0] == "FORCED-HIGH", "high typing changed")
    require(tuple(row[12] for row in ordinary) == ((1, 0, 1, 0), (135, 22, 113, 0)), "status counts changed")
    require(tuple(row[14] for row in ordinary) == (1, 113), "certificate counts changed")
    require(all((row[16], row[17], row[18]) == (0, 0, None) for row in ordinary), "ordinary residual changed")
    require(high[10] == F(264874483, 310625444220) > 0, "two-high gap changed")
    require(high[11] == F(80927993, 310625444220) > 0, "exact one-high gap changed")
    require((high[12], high[13], high[14], high[15]) == (
        (1612, 2340, 1736), 490, F(53, 176988), 6708,
    ), "high maximizer changed")

    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256,
        EXPECTED_ATLAS_SOURCE_SHA256, EXPECTED_ATLAS_OUTPUT_SHA256,
        CASES, tuple(row[12] for row in ordinary), tuple(row[14] for row in ordinary),
        high[10], high[11], 1599, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1608/z1604 terminal descent",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        f"atlas_source_sha256={file_sha256(ATLAS_SOURCE)}",
        f"atlas_output_sha256={file_sha256(ATLAS_OUTPUT)}",
        "inherited_cap=1611;atlas_gaps=1609..1611,1605..1607,1600..1603 empty;boundary_rows=z1608:2,z1604:1",
        "closure_partition=all_label_status:2;forced_high_scalar:1;empty:3;survivors:0",
        "ordinary_ledger=scalar:136;crude:22;status:114;status_survivors:0;exact_status_certificates:114;literal_packets:0;projected_kills:0",
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
        "band_consequence=all projected k=2 scalar rows with 1600<=z1<=1679 are empty",
        "consequence=projected k=2 first drift label z1<=1599",
        "scope_boundary=finite-exact intended THM-2980 band input complete;THM-2980 remains reserved pending integrated promotion audit;LRC(14) remains open",
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
