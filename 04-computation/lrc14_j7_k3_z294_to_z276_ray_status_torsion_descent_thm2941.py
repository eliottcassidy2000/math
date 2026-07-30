#!/usr/bin/env python3
"""Close every occupied projected k=3 slice from z1=294 through z1=276.

The pinned z297 located-torsion theorem leaves the next occupied atlas height
at 294.  Exact periodic-ray enumeration and the common 16-cell Hunter table
close the occupied levels 294, 291, 288, 287, 285, 278, and 276 outright.
At level 286 they leave eight states, all on E=(1,8,10,11,12,14).

The z286 residuals inherit the projected high-label obligation.  A
duplicate-permitting ray upper bound excludes two or more high labels, while
the unrestricted zero-high hostile relaxation passes all eight states; hence
the gate is both essential and forces exactly one high label.  Literal
enumeration below the wall, including negative-amplitude labels, leaves eight
one-high/two-low cases and two low pairs.  In each case the distinct residues
of complete fixed-safe body cells satisfy |S|>d/r for r=2 or r=4.  Two cells
therefore have effective phase order at most four for every unit high ray.
Their strict radius-1/14 danger arcs are disjoint, so all eight residuals are
empty uniformly in the high label.  There is no finite label horizon.

Together the eight occupied levels contain 1,549 quotient states: 659 crude
fibre-capacity kills, 882 independently replayed exact Farkas kills, and the
eight z286 located-torsion kills.  The next occupied atlas height is 275.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Z297_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
)
Z297_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.out"
)
ATLAS_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z294_to_z276_ray_status_torsion_descent_thm2941.out"
)

EXPECTED_Z297_SHA256 = (
    "6147066a29e20145bec47d87f2eb1bb054fa5e6ffeaa78e977affe2796e1a0a0"
)
EXPECTED_Z297_OUTPUT_SHA256 = (
    "7b7da4228c105847d482f22b54d27beeea237a61758c3295f71cb5584aedfa25"
)
EXPECTED_ATLAS_SHA256 = (
    "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
)
EXPECTED_SEMANTIC_SHA256 = (
    "a553a7e6717d3d41a8bdc3b8e799e5378a17bac27c2e6c40c24265f2a28ab8e8"
)
EXPECTED_M_HISTOGRAM = (
    (2, 24),
    (3, 86),
    (4, 213),
    (5, 254),
    (6, 86),
    (7, 187),
    (12, 24),
    (24, 8),
)

LEVELS = (294, 291, 288, 287, 286, 285, 278, 276)
EXPECTED_LEVEL_TOTALS = {
    294: (2, 0, 2, 0),
    291: (3, 3, 0, 0),
    288: (439, 247, 192, 0),
    287: (6, 0, 6, 0),
    286: (354, 154, 192, 8),
    285: (16, 7, 9, 0),
    278: (15, 0, 15, 0),
    276: (714, 248, 466, 0),
}
EXPECTED_TOTALS = (1549, 659, 882, 8)
INHERITED_LEDGER = 375_758
FINAL_LEDGER = 375_713

EXPECTED_COUNTS = {
    (294, (1, 4, 8, 10, 12, 14)): (2, 0, 2, 0),
    (291, (1, 3, 6, 9, 10, 12)): (2, 2, 0, 0),
    (291, (1, 6, 8, 9, 10, 12)): (1, 1, 0, 0),
    (288, (1, 2, 8, 10, 12, 14)): (3, 2, 1, 0),
    (288, (1, 8, 10, 11, 12, 14)): (301, 173, 128, 0),
    (288, (1, 8, 10, 12, 13, 14)): (16, 12, 4, 0),
    (288, (2, 8, 10, 11, 12, 14)): (119, 60, 59, 0),
    (287, (1, 2, 7, 10, 12, 14)): (1, 0, 1, 0),
    (287, (1, 4, 7, 10, 12, 14)): (1, 0, 1, 0),
    (287, (1, 4, 8, 10, 12, 14)): (3, 0, 3, 0),
    (287, (1, 7, 8, 10, 12, 14)): (1, 0, 1, 0),
    (287, (1, 8, 10, 12, 13, 14)): (0, 0, 0, 0),
    (286, (1, 4, 10, 11, 12, 14)): (0, 0, 0, 0),
    (286, (1, 8, 10, 11, 12, 14)): (354, 154, 192, 8),
    (286, (2, 8, 10, 11, 12, 14)): (0, 0, 0, 0),
    (286, (3, 8, 10, 11, 12, 14)): (0, 0, 0, 0),
    (285, (1, 2, 6, 8, 12, 14)): (1, 0, 1, 0),
    (285, (1, 2, 6, 9, 12, 14)): (5, 1, 4, 0),
    (285, (1, 3, 6, 9, 10, 12)): (2, 2, 0, 0),
    (285, (1, 3, 6, 9, 12, 14)): (2, 2, 0, 0),
    (285, (1, 4, 6, 8, 12, 14)): (1, 0, 1, 0),
    (285, (1, 4, 6, 9, 12, 14)): (4, 1, 3, 0),
    (285, (1, 6, 7, 9, 12, 14)): (1, 1, 0, 0),
    (278, (1, 2, 6, 8, 12, 14)): (2, 0, 2, 0),
    (278, (2, 3, 6, 8, 12, 14)): (1, 0, 1, 0),
    (278, (2, 4, 6, 8, 12, 14)): (2, 0, 2, 0),
    (278, (2, 6, 8, 10, 12, 14)): (10, 0, 10, 0),
    (276, (1, 2, 4, 8, 10, 12)): (1, 1, 0, 0),
    (276, (1, 2, 4, 8, 12, 14)): (10, 2, 8, 0),
    (276, (1, 2, 8, 10, 12, 14)): (5, 3, 2, 0),
    (276, (1, 3, 4, 8, 10, 12)): (4, 4, 0, 0),
    (276, (1, 3, 4, 8, 12, 14)): (3, 0, 3, 0),
    (276, (1, 4, 5, 8, 12, 14)): (65, 24, 41, 0),
    (276, (1, 4, 6, 8, 12, 14)): (1, 0, 1, 0),
    (276, (1, 4, 8, 10, 12, 14)): (283, 141, 142, 0),
    (276, (1, 4, 8, 12, 13, 14)): (3, 3, 0, 0),
    (276, (1, 8, 10, 12, 13, 14)): (0, 0, 0, 0),
    (276, (2, 3, 4, 8, 12, 14)): (5, 0, 5, 0),
    (276, (2, 4, 5, 8, 12, 14)): (71, 7, 64, 0),
    (276, (2, 4, 6, 8, 12, 14)): (4, 0, 4, 0),
    (276, (2, 4, 8, 10, 12, 14)): (232, 56, 176, 0),
    (276, (2, 5, 8, 10, 12, 14)): (5, 0, 5, 0),
    (276, (2, 8, 10, 11, 12, 14)): (0, 0, 0, 0),
    (276, (3, 4, 8, 9, 10, 12)): (1, 1, 0, 0),
    (276, (4, 5, 8, 10, 12, 14)): (21, 6, 15, 0),
}

Z286_BODY = (1, 8, 10, 11, 12, 14)
Z286_RESIDUALS = (
    (784, 1848, 5390, 5880),
    (784, 4620, 5390, 5880),
    (1848, 3920, 5390, 5880),
    (1848, 5390, 5880, 8624),
    (1848, 5390, 5880, 43120),
    (3920, 4620, 5390, 5880),
    (4620, 5390, 5880, 8624),
    (4620, 5390, 5880, 43120),
)
Z286_TWO_HIGH_GAP = Q(31003267, 10140586456)
Z286_QUALIFYING_HISTOGRAM = ((2, 6), (4, 2))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


for path, expected in (
    (Z297_PATH, EXPECTED_Z297_SHA256),
    (Z297_OUTPUT_PATH, EXPECTED_Z297_OUTPUT_SHA256),
    (ATLAS_PATH, EXPECTED_ATLAS_SHA256),
):
    if expected is not None:
        require(file_sha256(path) == expected, ("dependency changed", path))

z297 = load("z294_to_z276_terminal_engine", Z297_PATH)


def atlas_data():
    counts = Counter()
    rows = {level: [] for level in LEVELS}
    z286_gate_rows = []
    pattern = re.compile(r"^row=E=([0-9,]+);.*;z1=([0-9]+);")
    for line in ATLAS_PATH.read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        first = int(match.group(2))
        counts[first] += 1
        if first in rows:
            rows[first].append(body)
        if first == 286 and body == Z286_BODY:
            require("suffix=HIGH-TAIL:" in line, "z286 projected high gate disappeared")
            z286_gate_rows.append(line)
    tasks = tuple((first, body) for first in LEVELS for body in rows[first])
    require(tasks == tuple(EXPECTED_COUNTS), ("atlas body order changed", tasks))
    expected_heights = {level: sum(first == level for first, _body in tasks) for level in LEVELS}
    for level in range(276, 297):
        require(counts[level] == expected_heights.get(level, 0), ("atlas gap changed", level, counts[level]))
    require(counts[275] == 10, ("next frontier changed", counts[275]))
    require(len(z286_gate_rows) == 1, ("z286 high-gate row changed", z286_gate_rows))
    return tasks, tuple(sorted(counts.items()))


def evaluate_task(task):
    first, body = task
    z297.FIRST = first
    z297.ray.FIRST = first
    z297.EXPECTED_COUNTS = {body: EXPECTED_COUNTS[task]}
    return first, z297.evaluate(body)


def terminal_record(records):
    rows = [row for first, row in records if first == 286 and row[10]]
    require(len(rows) == 1 and rows[0][0] == Z286_BODY, ("z286 residual body changed", rows))
    residuals = tuple(packet[0] for packet in rows[0][10])
    require(residuals == Z286_RESIDUALS, ("z286 residual ledger changed", residuals))
    z297.FIRST = 286
    z297.ray.FIRST = 286
    z297.EXPECTED_GAPS = {Z286_BODY: Z286_TWO_HIGH_GAP}
    z297.EXPECTED_TERMINAL = {Z286_BODY: (8, 8, 2)}
    terminal = z297.analyze_residual_body(Z286_BODY, residuals)
    require(terminal[5] == Z286_TWO_HIGH_GAP > 0, ("two-high gap changed", terminal[5]))
    require((terminal[7], terminal[8], terminal[9]) == (8, 8, 2), "finite terminal ledger changed")
    require(terminal[12] == Z286_QUALIFYING_HISTOGRAM, ("least-r histogram changed", terminal[12]))
    require(terminal[13] == Z286_QUALIFYING_HISTOGRAM, ("effective-order histogram changed", terminal[13]))
    require(z297.THREE_ALIGNED_CAP == Q(36, 91) < 1, "aligned completion cap changed")
    return terminal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")
    require(
        "all seven projected k3 z1=297 body rows are empty"
        in Z297_OUTPUT_PATH.read_text(),
        "upstream z297 conclusion changed",
    )
    tasks, atlas_counts = atlas_data()
    pair_rows, control_marginals, control_caps, control_instances = (
        z297.ray.local.controls()
    )
    if args.processes == 1:
        records = tuple(evaluate_task(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            records = tuple(pool.map(evaluate_task, tasks))

    totals = tuple(sum(len(row[index]) for _first, row in records) for index in (7, 8, 9, 10))
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    level_totals = {
        level: tuple(
            sum(len(row[index]) for first, row in records if first == level)
            for index in (7, 8, 9, 10)
        )
        for level in LEVELS
    }
    require(level_totals == EXPECTED_LEVEL_TOTALS, ("level totals changed", level_totals))
    m_histogram = Counter()
    for _first, row in records:
        m_histogram.update(dict(row[11]))
    m_histogram = tuple(sorted(m_histogram.items()))
    if EXPECTED_M_HISTOGRAM is not None:
        require(m_histogram == EXPECTED_M_HISTOGRAM, ("M histogram changed", m_histogram))
    terminal = terminal_record(records)
    require(INHERITED_LEDGER - len(tasks) == FINAL_LEDGER, "ledger arithmetic changed")

    semantic_payload = (
        tasks,
        records,
        terminal,
        totals,
        tuple(sorted(level_totals.items())),
        m_histogram,
        atlas_counts,
        pair_rows,
        control_marginals,
        control_caps,
        control_instances,
        INHERITED_LEDGER,
        FINAL_LEDGER,
        EXPECTED_Z297_OUTPUT_SHA256,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=294 through z1=276 exact ray/status and located-torsion descent",
        f"z297_source_sha256={file_sha256(Z297_PATH)}",
        f"z297_output_sha256={file_sha256(Z297_OUTPUT_PATH)}",
        f"atlas_output_sha256={file_sha256(ATLAS_PATH)}",
        "scope=all 45 exact scalar-atlas body rows at occupied heights 294,291,288,287,286,285,278,276;no finite label horizon",
        f"level_totals={level_totals}",
        f"totals=states:{totals[0]};crude_kills:{totals[1]};status_kills:{totals[2]};torsion_residuals:{totals[3]}",
        f"status_M_histogram={dict(m_histogram)}",
        "independent_exact_farkas_checks=882/882:PASS",
        f"pair_overlap_exhaustive_controls={pair_rows}",
        "common_table_controls=coverable positive instance accepted;incompatible hostile instance rejected",
    ]
    for first, row in records:
        (
            body, L, high, first_d, trials, checks, signs, states, crude,
            status, packets, mhist, verified_instance_digest,
            verified_instance_count,
        ) = row
        lines.append(
            f"BODY;z1={first};E={body};L={L};high={high};d1={first_d};trials={trials};"
            f"checks={checks};signs={dict(signs)};states={len(states)};crude={len(crude)};"
            f"status={len(status)};residual={len(packets)};M={dict(mhist)};"
            f"verified_farkas_instance_sha256={verified_instance_digest};"
            f"verified_exact_farkas_instances={verified_instance_count}"
        )
        for packet in packets:
            lines.append(f"RESIDUAL;z1={first};E={body};row={packet}")

    (
        body, L, high, required, residual_count, gap, gap_witness,
        zero_count, case_count, pair_count, low_signs, high_checks,
        qualifying, effective, cases,
    ) = terminal
    lines.extend(
        (
            f"TERMINAL;z1=286;E={body};L={L};high={high};required_suffix={z297.ftext(required)};"
            f"residual_states={residual_count};two_high_gap={z297.ftext(gap)};gap_witness={gap_witness};"
            f"zero_high_passes={zero_count};one_high_cases={case_count};low_pairs={pair_count};"
            f"low_sign_census={dict(low_signs)};high_ray_checks={high_checks};"
            f"least_r={dict(qualifying)};effective_orders={dict(effective)}",
            "z286_logic=the inherited HIGH-TAIL gate excludes zero-high packets;the exact duplicate-permitting >=2-high gap forces exactly one high label",
            "torsion_pigeonhole=distinct fixed-safe residues S modulo high denominator d satisfy |S|>d/r for r=2 or4;two residues in one mod-(d/r) coset have effective order<=r",
            "phase_law=z=(L/d)u+hL with gcd(u,d)=1;unit multiplication preserves effective order;phase gap>=1/r>=1/4>1/7,so strict radius-1/14 dangers are disjoint",
            "projected_consequence=all8 z286 cases retain full drift-safe mass1>three-aligned union cap36/91",
        )
    )
    for case in cases:
        lines.append(f"CASE;z1=286;E={body};row={case}")
    lines.extend(
        (
            "atlas_exact_check=heights293,292,290,289,284..279,277 empty;next occupied z275 with10 bodies",
            f"ledger_decrement={INHERITED_LEDGER}-{len(tasks)}={FINAL_LEDGER};decrement counts atlas body rows,not quotient states",
            "conclusion=all projected k3 rows at occupied heights294..276 are empty;with pinned z297 closure,projected k3 cap<=275;next exact frontier=275",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
