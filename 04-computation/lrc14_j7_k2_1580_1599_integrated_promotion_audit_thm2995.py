#!/usr/bin/env python3
"""Integrated exact promotion audit for projected k=2 band 1580..1599.

This referee pins the five component source/output pairs, verifies their
dependency handoffs and 26 closure CASE records, checks the four translated
complete-section witnesses, and independently reruns the all-3003-body scalar
atlas as one monolithic 20-height band.  The atlas survivor set must equal the
component closure set row for row.

The resulting cap z1<=1579 concerns only the projected k=2 scalar necessary
sector.  This computation does not itself promote reserved THM-2995 and does
not prove LRC(14).
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import re
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
ATLAS_ENGINE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1810_1835_thm2941.py"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_1580_1599_integrated_promotion_audit_thm2995.out"

EXPECTED_ATLAS_ENGINE_SHA256 = "a09b13e994ad6ab35e3324bd336e773b435b07859e6a3c924b84ec77f3e2aced"
EXPECTED_ATLAS_PROFILE_SHA256 = "0f02d797f9a3bc39871d4a4ccf72c42525ca3ca4109773469b8fe95b780780d3"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "aab25b18f0bb9e2cb33eaa378176e720e1375c19891be9c97c67371f72ae2636"
EXPECTED_KEY_SHA256 = "3aa39633058336c148fb3fc730e44fe88fbe817827faaa4506b613be4add804c"
EXPECTED_CASE_LINE_SHA256 = "3cbdc8581a0f68a4edefaec32884dace147b4e246a61dc9ebb51a25c902102c1"
EXPECTED_MANIFEST_SHA256 = "f2816cc161b27f7b2df0c3cb1c99104bede7624e10c9dc19de85be78cc1c9811"
EXPECTED_SEMANTIC_SHA256 = "4ba1d451dd931bb86458d7559ad1c44da9ca7622dd5d753818efe1566251aff9"

START = 1580
END = 1599
EXPECTED_CANDIDATE_ROWS = 60_060
EXPECTED_HEIGHTS = ((1581, 3), (1586, 11), (1588, 1), (1590, 3), (1594, 2), (1595, 2), (1599, 4))
EXPECTED_COMPONENT_CASE_COUNTS = (4, 4, 3, 12, 3)
EXPECTED_CONSEQUENCE_CAPS = (1598, 1593, 1589, 1585, 1579)
EXPECTED_CONCLUSIONS = (("EMPTY", 9), ("SCALAR-EMPTY", 13), ("SECTION-CLOSURE-BELOW", 4))

COMPONENTS = (
    (
        "lrc14_j7_k2_z1599_forced_high_closure_thm2995",
        "97dc7a8bdf3103d1e3cf223fe862be217dec4fe128647bf710fb7d26fcf2bf84",
        "435e186056ab2727b518e1f4f24f08e6c76708464c6880ae89f8c92d33168159",
    ),
    (
        "lrc14_j7_k2_z1595_z1594_mixed_descent_thm2995",
        "53fdc5c3c146d4b66a2656d249ea2c8ddec6721342c69063b49c41976e2c021f",
        "7afd64dcbb294fca43ff8d254fcf72389592f8c414ce4f7f49319f1d4353dc7a",
    ),
    (
        "lrc14_j7_k2_z1590_status_projected_descent_thm2995",
        "04e00a86c6d5b3a7554bb837dfa7866cf1ca4149c7313958d4f8bd76a54a20f6",
        "d74cd340c7004b5e7d7c6f13d4365afb76be7fa6959cc7c01bdece2035d3c6bd",
    ),
    (
        "lrc14_j7_k2_z1588_z1586_mixed_section_descent_thm2995",
        "036f42af2e6d8ec40cc1cb233c2f2d0bdfd191777ba0e08048b4e62926c67dd4",
        "22c96c5d09087d2dab941346eceb278cc3408923802b6e39b8c446f1445e235d",
    ),
    (
        "lrc14_j7_k2_z1581_terminal_descent_thm2995",
        "4949b90b09163c651b11040ded94b1f9c3975005637a1f73f1ea6c6c461b6ccd",
        "01ea98702d51f9eed7c706592560031f6158167e2fbb9f54621edfc9985a32a9",
    ),
)

ROW_RE = re.compile(r"(?:^|;)z1=(\d+);E=([0-9]+(?:,[0-9]+){5})(?:;|$)")
CONCLUSION_RE = re.compile(r"(?:^|;)conclusion=([A-Z-]+)(?:;|$)")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(file_sha256(ATLAS_ENGINE) == EXPECTED_ATLAS_ENGINE_SHA256, "atlas engine changed")
ASPEC = spec_from_file_location("k2_1580_1599_integrated_atlas", ATLAS_ENGINE)
require(ASPEC is not None and ASPEC.loader is not None, "cannot load scalar-atlas engine")
ATLAS = module_from_spec(ASPEC)
ASPEC.loader.exec_module(ATLAS)


def configure_atlas():
    ATLAS.B.START = START
    ATLAS.B.END = END
    ATLAS.B.LABELS = np.arange(START, ATLAS.B.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return ATLAS.B.profile(body)


def parse_row(line):
    match = ROW_RE.search(line)
    require(match is not None, ("malformed row", line))
    return int(match.group(1)), tuple(map(int, match.group(2).split(",")))


def component_audit():
    texts = []
    cases_by_component = []
    case_lines = []
    conclusions = Counter()
    manifest = []
    section_rows = []
    for index, (name, source_hash, output_hash) in enumerate(COMPONENTS):
        source = ROOT / f"04-computation/{name}.py"
        output = ROOT / f"05-knowledge/results/{name}.out"
        require(file_sha256(source) == source_hash, (name, "source changed"))
        require(file_sha256(output) == output_hash, (name, "output changed"))
        text = output.read_text(encoding="utf-8")
        require(text.count("all_exact_controls=PASS") == 1, (name, "exact controls missing"))
        require(
            f"consequence=projected k=2 first drift label z1<={EXPECTED_CONSEQUENCE_CAPS[index]}" in text,
            (name, "consequence handoff changed"),
        )
        if index:
            require(COMPONENTS[index - 1][1] in text and COMPONENTS[index - 1][2] in text, (name, "dependency hash handoff changed"))
        lines = tuple(line for line in text.splitlines() if line.startswith("CASE;"))
        require(len(lines) == EXPECTED_COMPONENT_CASE_COUNTS[index], (name, "CASE count changed"))
        rows = []
        for line in lines:
            row = parse_row(line)
            conclusion_match = CONCLUSION_RE.search(line)
            require(conclusion_match is not None, (name, "CASE conclusion missing", row))
            conclusion = conclusion_match.group(1)
            require(conclusion in dict(EXPECTED_CONCLUSIONS), (name, "unknown CASE conclusion", conclusion))
            rows.append(row)
            conclusions[conclusion] += 1
        for line in text.splitlines():
            if line.startswith("HIGH-SECTION;"):
                require("failures=0" in line and line.endswith("conclusion=EMPTY"), (name, "translated-section witness changed"))
                section_rows.append(parse_row(line))
        texts.append(text)
        cases_by_component.append(tuple(rows))
        case_lines.extend(lines)
        manifest.append((name, source_hash, output_hash))

    require(tuple(sorted(conclusions.items())) == EXPECTED_CONCLUSIONS, "closure conclusion partition changed")
    require(sha256(repr(tuple(case_lines)).encode()).hexdigest() == EXPECTED_CASE_LINE_SHA256, "CASE transcript digest changed")
    require(sha256(repr(tuple(manifest)).encode()).hexdigest() == EXPECTED_MANIFEST_SHA256, "component manifest changed")
    section_cases = tuple(
        parse_row(line)
        for text in texts
        for line in text.splitlines()
        if line.startswith("CASE;") and "conclusion=SECTION-CLOSURE-BELOW" in line
    )
    require(tuple(section_rows) == section_cases and len(section_rows) == 4, "translated-section CASE/witness join changed")

    atlas_rows = tuple(parse_row(line) for line in texts[0].splitlines() if line.startswith("ATLAS;"))
    closure_rows = tuple(sorted(row for rows in cases_by_component for row in rows))
    require(len(closure_rows) == len(set(closure_rows)) == 26, "component closure set changed")
    require(tuple(sorted(atlas_rows)) == closure_rows and len(atlas_rows) == 26, "pinned atlas/closure join changed")
    require("band_consequence=all projected k=2 scalar rows with 1580<=z1<=1599 are empty" in texts[-1], "terminal band consequence changed")
    return closure_rows, tuple(cases_by_component), tuple(sorted(conclusions.items())), tuple(section_rows)


def render(profiles):
    require(len(profiles) == comb(14, 6) == 3003, "six-body universe changed")
    candidate_rows = sum(row[7] for row in profiles)
    survivors = tuple(
        sorted(
            (survivor for row in profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    keys = tuple((row[1], row[0]) for row in survivors)
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0])
        for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "monolithic candidate count changed")
    require(len(keys) == len(set(keys)) == 26, "monolithic survivor set changed size")
    require(heights == EXPECTED_HEIGHTS, ("monolithic height profile changed", heights))
    require(len(high_keys) == 17, "HIGH-TAIL routing count changed")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    key_hash = sha256(repr(keys).encode()).hexdigest()
    require(profile_hash == EXPECTED_ATLAS_PROFILE_SHA256, "monolithic atlas profile digest changed")
    require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "monolithic survivor digest changed")
    require(key_hash == EXPECTED_KEY_SHA256, "monolithic key digest changed")

    closure_rows, cases_by_component, conclusions, section_rows = component_audit()
    require(keys == closure_rows, "monolithic atlas and component closure sets differ")

    semantic_payload = (
        EXPECTED_ATLAS_ENGINE_SHA256,
        EXPECTED_MANIFEST_SHA256,
        START,
        END,
        candidate_rows,
        heights,
        high_keys,
        EXPECTED_COMPONENT_CASE_COUNTS,
        EXPECTED_CONSEQUENCE_CAPS,
        conclusions,
        section_rows,
        profile_hash,
        survivor_hash,
        key_hash,
        1579,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "integrated semantic digest changed")

    lines = [
        "LRC14 projected k=2 integrated promotion audit for z1=1580..1599",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"atlas_engine_sha256={file_sha256(ATLAS_ENGINE)}",
        f"component_manifest_sha256={EXPECTED_MANIFEST_SHA256};component_files={len(COMPONENTS)}",
        "universe=six_body_roots:3003;body_labels:1..14;ordered distinct later nonaligned labels;scalar necessary condition",
        f"monolithic_atlas_band={START}..{END};candidate_rows={candidate_rows};survivors={len(survivors)};height_counts={heights};high_tail_rows={len(high_keys)}",
        f"closure_partition=component_cases:26;component_case_counts={EXPECTED_COMPONENT_CASE_COUNTS};conclusions={conclusions};translated_section_witnesses={len(section_rows)}",
        "set_join=monolithic_atlas:26;pinned_atlas:26;component_cases:26;duplicates:0;missing:0;extra:0",
        f"handoff_caps={EXPECTED_CONSEQUENCE_CAPS};terminal_cap=1579",
        "hostile_audit=four nonpositive exact one-high CASE rows match four zero-failure translated complete-section witnesses;all other CASE records are certified EMPTY or SCALAR-EMPTY",
        "semantic_boundary=finite-exact projected k=2 scalar necessary sector only;does not close z1<=1579, other k sectors, the full six-body/seven-tail rung, or LRC(14)",
        "promotion_state=all reserved THM-2995 computational inputs frozen;canonical theorem promotion remains a separate proof-text audit",
        f"atlas_profile_sha256={profile_hash}",
        f"atlas_survivor_sha256={survivor_hash}",
        f"closure_key_sha256={key_hash}",
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
        profiles = tuple(atlas_profile(body) for body in bodies)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=configure_atlas) as pool:
            profiles = tuple(pool.imap(atlas_profile, bodies, chunksize=4))
    payload = render(profiles)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
