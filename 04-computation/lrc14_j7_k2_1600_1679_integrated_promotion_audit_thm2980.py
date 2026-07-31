#!/usr/bin/env python3
"""Integrated exact promotion audit for the projected k=2 band 1600..1679.

This referee does not rerun the expensive row-closing engines.  Instead it
pins every component source/output transcript, verifies the handoff chain and
all 68 CASE records, and independently reruns the all-3003-body scalar atlas
as one monolithic 80-height band.  The monolithic survivor set must agree
exactly, row for row, with the disjoint union of the component closures.

The audited consequence is only the projected k=2 scalar necessary-sector
cap z1<=1599.  It is not a proof of the full six-body/seven-tail rung or of
LRC(14), and this computation does not itself change reserved THM-2980.
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
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_1600_1679_integrated_promotion_audit_thm2980.out"

EXPECTED_ATLAS_ENGINE_SHA256 = "a09b13e994ad6ab35e3324bd336e773b435b07859e6a3c924b84ec77f3e2aced"
EXPECTED_ATLAS_PROFILE_SHA256 = "8ef2eead4eb68eeff33e38152c09648032e5652886594c71214c54bd55e39d83"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "b94748a736bc0069efb3b091b112527ad42a146b7eb549bd15f101658d2249d6"
EXPECTED_KEY_SHA256 = "ee1fdaf1e6b60daeaf7c4402b542db1f702bb559588e4e1c31a662cb77a73245"
EXPECTED_CASE_LINE_SHA256 = "2f604c6ba4cedf42323eed3d7b898edeb133c8a19d254d3c5fedc827a8798290"
EXPECTED_MANIFEST_SHA256 = "5fcbe3af0d01135f528eede15a1db7a9c5611c2b4ce47c8b7956c980f839c0da"
EXPECTED_SEMANTIC_SHA256 = "857b3594c55da9401a575754bb6a5164e6fd1a54f0a7c1837c2a595fa6bf30ef"

START = 1600
END = 1679
EXPECTED_CANDIDATE_ROWS = 240_240
EXPECTED_HEIGHTS = (
    (1604, 1), (1608, 2), (1612, 19), (1616, 1), (1618, 1),
    (1620, 8), (1625, 1), (1630, 1), (1631, 1), (1642, 1),
    (1645, 1), (1650, 12), (1656, 5), (1660, 1), (1665, 1),
    (1668, 9), (1670, 2), (1672, 1),
)
EXPECTED_COMPONENT_CASE_COUNTS = (3, 10, 1, 5, 12, 2, 2, 1, 8, 2, 19, 3)
EXPECTED_CONSEQUENCE_CAPS = (1668, 1660, 1656, 1655, 1649, 1641, 1629, 1624, 1619, 1615, 1611, 1599)
EXPECTED_CONCLUSIONS = (
    ("EMPTY", 40),
    ("LITERAL-CLOSURE-BELOW", 1),
    ("LITERAL-TORSION-EMPTY", 1),
    ("SCALAR-EMPTY", 22),
    ("SECTION-CLOSURE-BELOW", 4),
)

# Ordered in proof-handoff order.  These hashes pin both the executable
# certificate producers and the exact transcripts whose CASE rows are joined.
COMPONENTS = (
    (
        "lrc14_j7_k2_z1672_z1670_exact_descent_closure_thm2941",
        "b87b7e0e73207e026323269cb2493cf0793456b8292eb6effe94cdd823bb11e4",
        "bfcf783f84037f9b1417b10c4a1a6ec71a2aeec54c8f20383a1752a3d3149134",
    ),
    (
        "lrc14_j7_k2_z1668_composite_closure_thm2941",
        "9019e9da5f6e42b69652e6328c13e37c1704b165316f14030987f23fb2734ed6",
        "e24fdb898c0c576e68d33d2a921a96e8cad251676010a94959fef32d8c86515e",
    ),
    (
        "lrc14_j7_k2_z1660_status_descent_closure_thm2941",
        "b7d3344888034c1abc4cc5da72b73c8a81992e287ba7fc843b35cfbc6d5b305a",
        "1a1a85083ebb382d478139ce5206a20131ed412c5e47c61b26b58975cb2fed2e",
    ),
    (
        "lrc14_j7_k2_z1656_boundary_closure_thm2980",
        "593988b1963fc9f8406d7f7d0ae0f5cc0e9a8bd2604a5554522787fd2d0fa6b3",
        "c93f0fd172bb3a2a2148686566d4645fdd2db2a73e77c5048ba302a7c1bb1207",
    ),
    (
        "lrc14_j7_k2_z1650_shell_closure_thm2980",
        "d38929a73bcdd36824273b727ca10a19c842f5c70e35310d1c041e8d3d23377e",
        "a1a3cfe48e4e8034760c61b7fdbae282a10e052e3d0374a837788fd08e0e0b26",
    ),
    (
        "lrc14_j7_k2_z1645_z1642_status_descent_thm2980",
        "f1b4ee77b3749b52be5e67a3b733282b2f1ccb8c70197240af5a820c7a3786c8",
        "e98282c6266f474378ddeb725b18f7d2542b82e58ca5b6a6512f8a7600a39148",
    ),
    (
        "lrc14_j7_k2_z1631_z1630_status_descent_thm2980",
        "8200ba7c4bb40db669aa9b830cdbef0707b6349412e9f66ee26a51cc05018574",
        "ebed197453c407e38fd3f07cafd09c4addf868866759234ee6d25afbad6486eb",
    ),
    (
        "lrc14_j7_k2_z1625_forced_high_closure_thm2980",
        "ec98feabf60fe138df56647fc436ddec54b3f36b0704043c46af21e648e43935",
        "85947083c3bfdaeeaafd02aee9dac4a651b16c97fd80c1b84759d9380217fe7a",
    ),
    (
        "lrc14_j7_k2_z1620_mixed_shell_closure_thm2980",
        "eef6a5c0e7e988783357460c8f4144f3ce3a1a4e0cd7a52b949c6575661ce1e8",
        "bb4e5df57439de20b40b72481257a20e9ebacd278d42a3a25b947b8d8f39504b",
    ),
    (
        "lrc14_j7_k2_z1618_z1616_status_descent_thm2980",
        "3d8a0a3daa4775826beeee7d98cd26c8f3483e410239274c4f6d41c4d82faf8a",
        "20fdd1070210a4a61686dd13dc3bb1bee0c1ebc896b47400853593f44fa38fd2",
    ),
    (
        "lrc14_j7_k2_z1612_shell_closure_thm2980",
        "8861a4d02bab7a49c8f19a5b5c118b914e57ce504cfafe5926b6a647574fedb3",
        "77c8e4d359c0649143bdca3a6cfcfedb69625251f4eeaaa31599ce2d5d6fa99b",
    ),
    (
        "lrc14_j7_k2_z1608_z1604_terminal_descent_thm2980",
        "23298a37076b7b291fa2ba5dbbe849ed98b1dff33e94be8f2e06cc5e47a83c8a",
        "65d6d371e2ddd846167f9452b6e41c7f26535a4f4b2c70d20b39610194e750ea",
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
ASPEC = spec_from_file_location("k2_1600_1679_integrated_atlas", ATLAS_ENGINE)
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
            previous_source = COMPONENTS[index - 1][1]
            previous_output = COMPONENTS[index - 1][2]
            require(previous_source in text and previous_output in text, (name, "dependency hash handoff changed"))
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
        texts.append(text)
        cases_by_component.append(tuple(rows))
        case_lines.extend(lines)
        manifest.append((name, source_hash, output_hash))

    require(tuple(sorted(conclusions.items())) == EXPECTED_CONCLUSIONS, "closure conclusion partition changed")
    require(sha256(repr(tuple(case_lines)).encode()).hexdigest() == EXPECTED_CASE_LINE_SHA256, "CASE transcript digest changed")
    require(sha256(repr(tuple(manifest)).encode()).hexdigest() == EXPECTED_MANIFEST_SHA256, "component manifest changed")

    top_atlas = tuple(parse_row(line) for line in texts[0].splitlines() if line.startswith("ATLAS;"))
    lower_atlas = tuple(parse_row(line) for line in texts[4].splitlines() if line.startswith("ATLAS;"))
    expected_top = tuple(sorted(cases_by_component[0] + tuple(row for row in cases_by_component[1] if row[0] == 1668)))
    expected_lower = tuple(sorted(row for rows in cases_by_component[4:] for row in rows))
    require(tuple(sorted(top_atlas)) == expected_top and len(top_atlas) == 12, "upper atlas/closure join changed")
    require(tuple(sorted(lower_atlas)) == expected_lower and len(lower_atlas) == 49, "lower atlas/closure join changed")

    upper = tuple(row for rows in cases_by_component[:4] for row in rows)
    lower = tuple(row for rows in cases_by_component[4:] for row in rows)
    require(len(upper) == 19 and len(lower) == 49, "19+49 component split changed")
    require(set(upper).isdisjoint(lower), "upper/lower closure sets overlap")
    all_rows = tuple(sorted(upper + lower))
    require(len(all_rows) == len(set(all_rows)) == 68, "closure rows are not a 68-element set")
    return all_rows, tuple(cases_by_component), tuple(sorted(conclusions.items()))


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
    require(len(keys) == len(set(keys)) == 68, "monolithic survivor set changed size")
    require(heights == EXPECTED_HEIGHTS, ("monolithic height profile changed", heights))
    require(len(high_keys) == 28, "HIGH-TAIL routing count changed")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    key_hash = sha256(repr(keys).encode()).hexdigest()
    require(profile_hash == EXPECTED_ATLAS_PROFILE_SHA256, "monolithic atlas profile digest changed")
    require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "monolithic survivor digest changed")
    require(key_hash == EXPECTED_KEY_SHA256, "monolithic key digest changed")

    closure_rows, cases_by_component, conclusions = component_audit()
    require(keys == closure_rows, "monolithic atlas and closure CASE sets differ")

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
        profile_hash,
        survivor_hash,
        key_hash,
        1599,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "integrated semantic digest changed")

    lines = [
        "LRC14 projected k=2 integrated promotion audit for z1=1600..1679",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"atlas_engine_sha256={file_sha256(ATLAS_ENGINE)}",
        f"component_manifest_sha256={EXPECTED_MANIFEST_SHA256};component_files={len(COMPONENTS)}",
        "universe=six_body_roots:3003;body_labels:1..14;ordered distinct later nonaligned labels;scalar necessary condition",
        f"monolithic_atlas_band={START}..{END};candidate_rows={candidate_rows};survivors={len(survivors)};height_counts={heights};high_tail_rows={len(high_keys)}",
        f"closure_partition=upper_thm2941:19;reserved_band_components:49;component_case_counts={EXPECTED_COMPONENT_CASE_COUNTS};conclusions={conclusions}",
        "set_join=monolithic_atlas:68;component_cases:68;duplicates:0;missing:0;extra:0;upper_lower_overlap:0",
        f"handoff_caps={EXPECTED_CONSEQUENCE_CAPS};terminal_cap=1599",
        "hostile_audit=positive and zero-high scalar controls remain typed inside their component transcripts;only certified EMPTY/SCALAR-EMPTY/LITERAL-CLOSURE-BELOW/LITERAL-TORSION-EMPTY/SECTION-CLOSURE-BELOW CASE records enter the set join",
        "semantic_boundary=finite-exact projected k=2 scalar necessary sector only;does not close z1<=1599, other k sectors, the full six-body/seven-tail rung, or LRC(14)",
        "promotion_state=all reserved THM-2980 computational inputs frozen;canonical theorem promotion remains a separate proof-text audit",
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
