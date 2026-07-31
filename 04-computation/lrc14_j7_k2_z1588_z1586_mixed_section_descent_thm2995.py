#!/usr/bin/env python3
"""Exact projected k=2 mixed descent through z1=1588 and z1=1586.

The pinned atlas has one ordinary row at 1588 and eleven rows at 1586.
Three ordinary rows close by exact status/projected certificates.  Five
HIGH-TAIL rows are scalar-empty.  The remaining four one-high rows close
pointwise by the sharp translated complete-section capacity gate, with no
residue-pair choice.  The cap therefore drops to 1585.
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
BASE_NAME = "lrc14_j7_k2_z1590_status_projected_descent_thm2995"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1588_z1586_mixed_section_descent_thm2995.out"

EXPECTED_BASE_SOURCE_SHA256 = "04e00a86c6d5b3a7554bb837dfa7866cf1ca4149c7313958d4f8bd76a54a20f6"
EXPECTED_BASE_OUTPUT_SHA256 = "d74cd340c7004b5e7d7c6f13d4365afb76be7fa6959cc7c01bdece2035d3c6bd"
EXPECTED_PROFILE_SHA256 = "33a26ff1b960b28d89a132ab0815e94b812696e53399b6f4f895437489f239bc"
EXPECTED_SECTION_SHA256 = "384eff9bfaac0a8f1d1e16c1add095c4960263ef2edc38f04b7ae6a49a82be90"
EXPECTED_SEMANTIC_SHA256 = "49495a8c600f655aa109da196c20daac22b67d8187d75ecd1ba00d5b269ae25b"

STANDARD_BODY = (1, 4, 8, 10, 12, 14)
SECOND_ORDINARY_BODY = (2, 4, 8, 10, 12, 14)
CASES = (
    (1588, STANDARD_BODY),
    (1586, (1, 2, 10, 12, 13, 14)),
    (1586, (1, 3, 10, 12, 13, 14)),
    (1586, STANDARD_BODY),
    (1586, (1, 4, 10, 12, 13, 14)),
    (1586, (1, 6, 8, 10, 13, 14)),
    (1586, (1, 6, 10, 12, 13, 14)),
    (1586, (1, 8, 10, 12, 13, 14)),
    (1586, (1, 10, 11, 12, 13, 14)),
    (1586, (2, 3, 10, 12, 13, 14)),
    (1586, SECOND_ORDINARY_BODY),
    (1586, (2, 8, 10, 12, 13, 14)),
)
ORDINARY_CASES = (CASES[0], CASES[3], CASES[10])
HARD_CASES = (CASES[1], CASES[4], CASES[7], CASES[8])
EXPECTED_ORDINARY = {
    CASES[0]: ((1, 0, 1, 0), 1, 0, None, 0),
    CASES[3]: ((25, 4, 21, 0), 21, 0, None, 0),
    CASES[10]: ((5, 0, 2, 3), 2, 3, F(1296, 5551), 1),
}
EXPECTED_SECTION_COUNTS = {
    HARD_CASES[0]: (1, 1, 0),
    HARD_CASES[1]: (2, 2, 0),
    HARD_CASES[2]: (175, 175, 0),
    HARD_CASES[3]: (1007, 941, 66),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1590 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1590 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1589" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1590 handoff changed",
)
SPEC = spec_from_file_location("k2_z1588_z1586_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1590 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE
SECTION_AUDIT = BASE.BASE.BASE.ENGINE_PARENT.BASE.high_section_audit
ATLAS_KEYS = BASE.BASE.ATLAS_KEYS
require(all(case in ATLAS_KEYS for case in CASES), "z1588/z1586 atlas rows changed")
require(not any(first in {1587, 1589} for first, _body in ATLAS_KEYS), "z1587/z1589 gap changed")


def closure_profile(case):
    return ENGINE.ordinary_profile(case) if case in ORDINARY_CASES else ENGINE.high_profile(case)


def section_profile(case):
    return SECTION_AUDIT(case)


def render(rows, audits):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "closure universe changed")
    ordinary = tuple(row for row in rows if row[0] == "ALL-LABEL")
    high = tuple(row for row in rows if row[0] == "FORCED-HIGH")
    scalar_high = tuple(row for row in high if row[11] > 0)
    hard_high = tuple(row for row in high if row[11] <= 0)
    require((len(ordinary), len(scalar_high), len(hard_high)) == (3, 5, 4), "route split changed")
    require(tuple((row[1], row[2]) for row in hard_high) == HARD_CASES, "hard one-high universe changed")
    require(tuple(audit[0] for audit in audits) == HARD_CASES, "section audit universe changed")

    for row in ordinary:
        case = (row[1], row[2])
        counts, certificates, packets, margin, prefix = EXPECTED_ORDINARY[case]
        require(row[12] == counts and row[14] == certificates, (case, "status ledger changed"))
        require(row[16] == row[17] == packets, (case, "packet ledger changed"))
        require((row[18], row[19]) == (margin, prefix), (case, "projected control changed"))
    require(all(row[10] > 0 for row in high), "two-high hostile scalar bound changed")
    for audit in audits:
        case = audit[0]
        classes, count_gate, residue_gate = EXPECTED_SECTION_COUNTS[case]
        require((audit[3], audit[7], audit[8]) == (classes, count_gate, residue_gate), (case, "section ledger changed"))
        require(not audit[9] and audit[7] + audit[8] == audit[3], (case, "section failure"))

    ordinary_totals = tuple(sum(row[12][index] for row in ordinary) for index in range(4))
    require(ordinary_totals == (31, 4, 24, 3), "ordinary aggregate changed")
    require(sum(row[14] for row in ordinary) == 24, "certificate total changed")
    require(sum(row[16] for row in ordinary) == sum(row[17] for row in ordinary) == 3, "packet total changed")
    require(sum(audit[3] for audit in audits) == 1185, "section class total changed")
    require(sum(audit[7] for audit in audits) == 1119 and sum(audit[8] for audit in audits) == 66, "section gate total changed")

    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    section_hash = sha256(repr(tuple(audits)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    if EXPECTED_SECTION_SHA256 is not None:
        require(section_hash == EXPECTED_SECTION_SHA256, "section digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256,
        EXPECTED_BASE_OUTPUT_SHA256,
        CASES,
        ORDINARY_CASES,
        HARD_CASES,
        EXPECTED_ORDINARY,
        EXPECTED_SECTION_COUNTS,
        tuple((row[1], row[2], row[11]) for row in scalar_high),
        tuple((audit[0], audit[3], audit[7], audit[8], audit[14], audit[17]) for audit in audits),
        1585,
        profile_hash,
        section_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact mixed translated-section descent through z1588 and z1586",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1589;atlas_gaps=1589,1587 empty;boundary_rows=z1588:1,z1586:11",
        "closure_partition=all_label:3;forced_high_scalar:5;forced_high_translated_section:4;empty:12;survivors:0",
        "ordinary_ledger=scalar:31;crude:4;status:24;status_survivors:3;exact_status_certificates:24;literal_packets:3;projected_kills:3",
        "translated_section_ledger=classes:1185;count_gate:1119;residue_gate:66;failures:0",
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
            conclusion = "SCALAR-EMPTY" if row[11] > 0 else "SECTION-CLOSURE-BELOW"
            lines.append(
                f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
                f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
                f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
                f"amplitude_sha256={row[16]};conclusion={conclusion}"
            )
    for audit in audits:
        lines.append(
            f"HIGH-SECTION;z1={audit[0][0]};E={','.join(map(str,audit[0][1]))};"
            f"L={audit[1]};complete_cells={audit[2]};classes={audit[3]};"
            f"low_triples={audit[4]};safe_cell_range={audit[5]};"
            f"minimum_excess={ftext(audit[6])};count_gate={audit[7]};"
            f"residue_gate={audit[8]};failures={len(audit[9])};"
            f"active_denominators={audit[10]};minimum_count_surplus={audit[12]};"
            f"minimum_residue_surplus={audit[13]};section_sha256={audit[14]};"
            f"amplitude_sha256={audit[15]};zero_supremum_classes={audit[16]};"
            f"scalar_vector_control_sha256={audit[17]};conclusion=EMPTY"
        )
    lines.extend((
        "hostile_controls=four exact one-high gaps are nonpositive;all 1185 classes close pointwise by the sharp translated ceil(d/7) complete-section gate, including 66 residue rather than raw-count gates",
        "consequence=projected k=2 first drift label z1<=1585",
        "scope_boundary=finite-exact THM-2995 input only;remaining 1580..1585 band, other sectors, and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"section_sha256={section_hash}",
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
        rows = tuple(closure_profile(case) for case in CASES)
        audits = tuple(section_profile(case) for case in HARD_CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = tuple(pool.imap(closure_profile, CASES))
            audits = tuple(pool.imap(section_profile, HARD_CASES))
    payload = render(rows, audits)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
