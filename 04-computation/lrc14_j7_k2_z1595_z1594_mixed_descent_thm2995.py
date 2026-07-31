#!/usr/bin/env python3
"""Exact projected k=2 mixed descent through z1=1595 and z1=1594.

The pinned 1580..1599 atlas has two rows at each height.  At each height one
ordinary row closes by exact common-status/projected certificates and one
HIGH-TAIL row is scalar-empty by its exact one-high ray maximum.  Together
with the empty atlas gap 1596..1598 this lowers the projected cap to 1593.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import re
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1599_forced_high_closure_thm2995"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1595_z1594_mixed_descent_thm2995.out"

EXPECTED_BASE_SOURCE_SHA256 = "97dc7a8bdf3103d1e3cf223fe862be217dec4fe128647bf710fb7d26fcf2bf84"
EXPECTED_BASE_OUTPUT_SHA256 = "435e186056ab2727b518e1f4f24f08e6c76708464c6880ae89f8c92d33168159"
EXPECTED_PROFILE_SHA256 = "e6e15becad5e29ab6bf7526daaeafa9b2ecd5d1cece3ed52e123613fbcdae207"
EXPECTED_SEMANTIC_SHA256 = "804b8363a34ecd2ef20bf238bdf7dd294b2e97ea273868aa587edf65cf66cd34"

STANDARD_BODY = (1, 4, 8, 10, 12, 14)
CASES = (
    (1595, STANDARD_BODY),
    (1595, (1, 8, 10, 12, 13, 14)),
    (1594, STANDARD_BODY),
    (1594, (1, 4, 10, 12, 13, 14)),
)
EXPECTED_ORDINARY_COUNTS = {
    CASES[0]: ((22, 3, 0, 19), 0, 23, F(16, 273)),
    CASES[2]: ((24, 4, 20, 0), 20, 0, None),
}
EXPECTED_HIGH_GAPS = {
    CASES[1]: (F(2515378673, 3508202257560), F(23779715761, 585869777012520)),
    CASES[3]: (F(3253576657, 3695051926020), F(355145689, 1231683975340)),
}
ROW_RE = re.compile(r"(?:^|;)z1=(\d+);E=([0-9]+(?:,[0-9]+){5})(?:;|$)")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1599 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1599 output changed")
BASE_TEXT = BASE_OUTPUT.read_text(encoding="utf-8")
require("consequence=projected k=2 first drift label z1<=1598" in BASE_TEXT, "z1599 handoff changed")
ATLAS_KEYS = []
for line in BASE_TEXT.splitlines():
    if line.startswith("ATLAS;"):
        match = ROW_RE.search(line)
        require(match is not None, "malformed inherited ATLAS row")
        ATLAS_KEYS.append((int(match.group(1)), tuple(map(int, match.group(2).split(",")))))
require(all(case in ATLAS_KEYS for case in CASES), "z1595/z1594 atlas rows changed")
require(not any(1596 <= first <= 1598 for first, _body in ATLAS_KEYS), "inherited descent gap changed")

SPEC = spec_from_file_location("k2_z1595_z1594_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1599 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE


def profile(case):
    return ENGINE.ordinary_profile(case) if case[1] == STANDARD_BODY else ENGINE.high_profile(case)


def render(rows):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "closure universe changed")
    ordinary = tuple(row for row in rows if row[0] == "ALL-LABEL")
    high = tuple(row for row in rows if row[0] == "FORCED-HIGH")
    require((len(ordinary), len(high)) == (2, 2), "route partition changed")
    for row in ordinary:
        case = (row[1], row[2])
        counts, certificates, packets, margin = EXPECTED_ORDINARY_COUNTS[case]
        require(row[12] == counts and row[14] == certificates, (case, "status ledger changed"))
        require(row[16] == row[17] == packets and row[18] == margin, (case, "projected ledger changed"))
    for row in high:
        case = (row[1], row[2])
        require((row[10], row[11]) == EXPECTED_HIGH_GAPS[case], (case, "forced-high gaps changed"))
        require(row[10] > 0 and row[11] > 0, (case, "forced-high branch survived"))

    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_BASE_SOURCE_SHA256,
        EXPECTED_BASE_OUTPUT_SHA256,
        CASES,
        EXPECTED_ORDINARY_COUNTS,
        EXPECTED_HIGH_GAPS,
        1593,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact mixed descent through z1595 and z1594",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1598;atlas_gap=1596..1598 empty;boundary_rows=z1595:2,z1594:2",
        "closure_partition=ordinary_all_label:2;forced_high_scalar_empty:2;empty:4;survivors:0",
        "ordinary_ledger=scalar:46;crude:7;status:20;status_survivors:19;exact_status_certificates:20;literal_packets:23;projected_kills:23",
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
        "hostile_controls=the z1595 ordinary row has 19 status survivors requiring 23 exact projected packets;both HIGH-TAIL rows retain strictly positive exact one-high gaps",
        "consequence=projected k=2 first drift label z1<=1593",
        "scope_boundary=finite-exact THM-2995 input only;remaining 1580..1593 band, other sectors, and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
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
