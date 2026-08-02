#!/usr/bin/env python3
"""Exact projected-k3 z1=232 screen and complete-cell descent.

The promoted THM-3078 evaluator is reused at task level with ``first=232``.
All three rows are recomputed from the pinned atlas rather than by changing a
top-level frontier constant.  The unique residual body is then passed through
the exact duplicate-two-high and one-high terminal.  Finally, a local hostile
audit rebuilds the fixed-safe carrier directly on the full Z/LZ grid, compares
it with both inherited implementations, and recomputes every residue support.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3078 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.py"
)
OUTPUT_3078 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.out"
)
SOURCE_3102 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.py"
)
OUTPUT_3102 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.out"
)

SOURCE_3078_SHA256 = "2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168"
OUTPUT_3078_SHA256 = "d5fc52e922b7e083f5cbe5aa5a6066304c4e5c8c7963904dd3b649477caf5e42"
SEMANTIC_3078_SHA256 = "749a2292da58925c9653ce919d5bab6b374f64f8b713ac52efa1f18fba74c918"
SOURCE_3102_SHA256 = "3db29155f4c1332c91605e5589dfdc19319eb81dc309a33ff8af161abb39a036"
OUTPUT_3102_SHA256 = "f358eae14c52783b561b8b799b02fb07f2988ef6b1b9f6cd33f3030ce177727d"
SEMANTIC_3102_SHA256 = "6cf01affce52a5dcc67a8634da815780e3d357e72b295a7c4951211c6f12b0da"

LEVEL = 232
NEXT_LEVEL = 231
NEIGHBOR_COUNTS = {
    233: (62, 45, 17),
    232: (3, 3, 0),
    231: (9, 5, 4),
}
EXPECTED_ROW_ORDER_SHA256 = (
    "9bda2b6c4582541b8303156c2c9a47bf8b44d6d0ca9ceeb0b30d78860eebb796"
)
EXPECTED_NEXT_ROW_ORDER_SHA256 = (
    "ff12b3d3c41a10aa38da0a1483bfb506b4c6a9a771c2325fb4df1df867088f5e"
)
EXPECTED_SCREEN = (147, 58, 65, 24)
EXPECTED_AUDIT = (0, 65)
EXPECTED_SCREEN_RECORD_SHA256 = (
    "3ca9c8569f05d8adf0ff04586347700597d7875794525452b13031e15acd909a"
)
EXPECTED_TERMINAL = (1, 1, 1, 23, 24, 23, 1, 0, 0, 0)
EXPECTED_TERMINAL_RECORD_SHA256 = (
    "56cd0b886afdb39ca3423c0cc5240f06298075138f6336ab31e65fbd65e2907c"
)
EXPECTED_HOSTILE_RECORD_SHA256 = (
    "828f0631b5743d0ecbf8f767a5f187bc373cb7aec20bbdf4850f3b81fc801772"
)
EXPECTED_CARRIER_SHA256 = (
    "656bb385012dd08368a49ef6d365dbaa34362d6e6707313afbfb20128e4945fa"
)
EXPECTED_SEMANTIC_SHA256 = (
    "a14adcd1e52323baf2b791f55d5846c4fbd422e3441fa2411138bdd98acca3d3"
)

LEDGER_BEFORE = 374325
LAYER_ROWS = 3
LEDGER_AFTER = 374322
NEXT_CAP = 231

RESIDUAL_BODY = (1, 9, 10, 11, 12, 14)
RESIDUAL_BANK = (
    (8, 9702, 10780, 24255),
    (24, 9702, 10780, 24255),
    (40, 9702, 10780, 24255),
    (56, 9702, 10780, 24255),
    (72, 9702, 10780, 24255),
    (88, 9702, 10780, 24255),
    (120, 9702, 10780, 24255),
    (168, 9702, 10780, 24255),
    (264, 9702, 10780, 24255),
    (280, 9702, 10780, 24255),
    (360, 9702, 10780, 24255),
    (440, 9702, 10780, 24255),
    (504, 9702, 10780, 24255),
    (616, 9702, 10780, 24255),
    (792, 9702, 10780, 24255),
    (840, 9702, 10780, 24255),
    (1320, 9702, 10780, 24255),
    (1848, 9702, 10780, 24255),
    (2520, 9702, 10780, 24255),
    (3080, 9702, 10780, 24255),
    (3960, 9702, 10780, 24255),
    (5544, 9702, 10780, 24255),
    (9240, 9702, 10780, 24255),
    (9702, 10780, 24255, 27720),
)
EXPECTED_RESIDUAL_SHA256 = (
    "588fa690d15f18a29083f40f422507462f895c8880df0fb75c6ea8c394189a66"
)
EXPECTED_GAP = "918906100793/174485171984760"
EXPECTED_CASE_SHA256 = (
    "f44f893f55db7160d032758170d876417de1fdd2ac1e69268fcea0eda295b495"
)
EXPECTED_ENDPOINT = ((234, 260), 37702, (0, 14850, 14, "L"))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    normalized = payload.replace(b"\r\n", b"\n")
    require(b"\r" not in normalized, f"bare CR in {path}")
    return hashlib.sha256(normalized).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(sha(SOURCE_3078) == SOURCE_3078_SHA256, "THM-3078 source changed")
require(sha(OUTPUT_3078) == OUTPUT_3078_SHA256, "THM-3078 output changed")
require(
    f"semantic_sha256={SEMANTIC_3078_SHA256}" in OUTPUT_3078.read_text(),
    "THM-3078 semantic changed",
)
require(sha(SOURCE_3102) == SOURCE_3102_SHA256, "THM-3102 source changed")
require(sha(OUTPUT_3102) == OUTPUT_3102_SHA256, "THM-3102 output changed")
require(
    f"semantic_sha256={SEMANTIC_3102_SHA256}" in OUTPUT_3102.read_text(),
    "THM-3102 semantic changed",
)
require(
    "promotion_consequence=ledger 374387-62=374325;"
    "projected_k3_cap:z1<=232;next_layer:z1=232_rows:3"
    in OUTPUT_3102.read_text(),
    "THM-3102 ledger consequence changed",
)
require(
    hashlib.sha256(repr(RESIDUAL_BANK).encode()).hexdigest()
    == EXPECTED_RESIDUAL_SHA256,
    "residual-bank transcription",
)

thm = load("thm3078_z232_descent_base", SOURCE_3078)
eng = thm.eng


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    row_lines = tuple(
        line for line in thm.ATLAS.read_text().splitlines() if line.startswith("row=")
    )
    require(len(row_lines) == 6060, ("atlas row-line count", len(row_lines)))
    census = defaultdict(lambda: [0, 0, 0])
    selected = []
    next_rows = []
    for line in row_lines:
        match = pattern.match(line)
        require(match is not None, ("unparsed atlas row", line))
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        derived_high = (
            thm.screen_engine.base.WALL.numerator
            * L
            // thm.screen_engine.base.WALL.denominator
            + 1
        )
        require(high == derived_high, (body, L, high, "high floor"))
        if first in NEIGHBOR_COUNTS:
            wall = first < high
            census[first][0] += 1
            census[first][1 if wall else 2] += 1
            row = (first, body, L, high, wall)
            if first == LEVEL:
                selected.append(row)
            elif first == NEXT_LEVEL:
                next_rows.append((body, L, high, wall))
    derived = {level: tuple(census[level]) for level in NEIGHBOR_COUNTS}
    require(derived == NEIGHBOR_COUNTS, ("neighbor census", derived))
    require(len(selected) == 3 and all(task[-1] for task in selected), selected)
    require(len(next_rows) == 9, next_rows)
    order = tuple((task[1], task[2], task[3], task[4]) for task in selected)
    require(
        hashlib.sha256(repr(order).encode()).hexdigest()
        == EXPECTED_ROW_ORDER_SHA256,
        "z232 row-order digest",
    )
    next_sha = hashlib.sha256(repr(tuple(next_rows)).encode()).hexdigest()
    if EXPECTED_NEXT_ROW_ORDER_SHA256 is not None:
        require(next_sha == EXPECTED_NEXT_ROW_ORDER_SHA256, next_sha)
    return tuple(selected), derived, order, tuple(next_rows), next_sha


def screen_worker(task):
    return thm.worker_screen(task)


def run_screen(tasks, processes):
    if processes <= 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        context = mp.get_context("spawn")
        with context.Pool(processes=processes) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost or duplicated screen task")
    return tuple(by_task[task] for task in tasks)


def screen_audit(rows):
    totals = tuple(sum(row[index] for row in rows) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    audit = (sum(row[19] for row in rows), sum(row[20] for row in rows))
    require(audit == EXPECTED_AUDIT, audit)
    require(all(row[16] == row[11] for row in rows), "unverified status row")
    record_sha = hashlib.sha256(repr(rows).encode()).hexdigest()
    require(record_sha == EXPECTED_SCREEN_RECORD_SHA256, record_sha)
    residual_rows = tuple(row for row in rows if row[12])
    require(len(residual_rows) == 1, residual_rows)
    require(
        residual_rows[0][1] == RESIDUAL_BODY
        and residual_rows[0][13] == RESIDUAL_BANK,
        residual_rows,
    )
    return totals, audit, record_sha, residual_rows


def run_terminals(residual_rows):
    terminals = tuple(
        thm.terminal_probe((LEVEL, row[1], row[13])) for row in residual_rows
    )
    totals = (
        len(terminals),
        sum(row[19] for row in terminals),
        sum(row[20] for row in terminals),
        sum(row[7] for row in terminals),
        sum(row[8] for row in terminals),
        sum(row[10] for row in terminals),
        sum(row[11] for row in terminals),
        sum(row[12] for row in terminals),
        sum(row[14] for row in terminals),
        sum(row[15] for row in terminals),
    )
    require(totals == EXPECTED_TERMINAL, totals)
    record_sha = hashlib.sha256(repr(terminals).encode()).hexdigest()
    require(record_sha == EXPECTED_TERMINAL_RECORD_SHA256, record_sha)
    terminal = terminals[0]
    require(terminal[1] == RESIDUAL_BODY and terminal[4] == 24, terminal)
    require(ftext(terminal[5]) == EXPECTED_GAP and terminal[5] > 0, terminal[5])
    require(
        (
            terminal[7],
            terminal[8],
            terminal[9],
            terminal[10],
            terminal[11],
            terminal[12],
            terminal[13],
            terminal[14],
            terminal[15],
            terminal[16],
            terminal[17],
            terminal[19],
            terminal[20],
        )
        == (23, 24, 1, 23, 1, 0, 24, 0, 0, 1, EXPECTED_CASE_SHA256, True, True),
        "terminal profile",
    )
    return terminals, totals, record_sha


def direct_full_grid_safe_cells(body, stream, labels):
    fixed = (*body, LEVEL, *labels)
    cells = []
    for cell in range(stream.L):
        clean = True
        for label in fixed:
            residue = cell * label % stream.L
            if 14 * residue < stream.L or 14 * (residue + label) > 13 * stream.L:
                clean = False
                break
        if clean:
            cells.append(cell)
    return tuple(cells)


def hostile_terminal_audit(residual_rows):
    screen_row = residual_rows[0]
    body = screen_row[1]
    residual = screen_row[13]
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(body)
    needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, _signs, _checks = eng.build_literal_tables(stream, needed)
    cases = eng.one_high_cases(stream, residual, low, high)
    require(len(cases) == 24, len(cases))

    cell_cache = {}
    carrier_records = []
    records = []
    for ds, high_d, low_rows, excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        if labels not in cell_cache:
            vector = tuple(map(int, thm.vector_fixed_safe_cells(stream, labels).tolist()))
            scalar = eng.fixed_safe_cells(stream, labels)
            direct = direct_full_grid_safe_cells(body, stream, labels)
            require(vector == scalar == direct, (body, labels, "three-way cells"))
            cell_cache[labels] = direct
            carrier_sha = hashlib.sha256(repr(direct).encode()).hexdigest()
            if EXPECTED_CARRIER_SHA256 is not None:
                require(carrier_sha == EXPECTED_CARRIER_SHA256, carrier_sha)
            carrier_records.append((body, labels, len(direct), carrier_sha))
        cells = cell_cache[labels]
        require(stream.L % high_d == 0, (body, high_d, "nondivisor"))
        support = tuple(sorted({cell % high_d for cell in cells}))
        kappa = (high_d + 6) // 7
        capacity = stream.L // high_d
        coarse = (len(cells) + capacity - 1) // capacity
        certificate = thm.hybrid_translated_certificate(
            np.array(cells, dtype=np.int64), stream.L, high_d
        )
        require(certificate[8], (body, ds, high_d, certificate))
        records.append(
            (
                ds,
                high_d,
                low_rows,
                excess,
                len(cells),
                capacity,
                coarse,
                len(support),
                kappa,
                len(support) - kappa,
                certificate[7],
                hashlib.sha256(repr(support).encode()).hexdigest(),
            )
        )

    endpoints = []
    for labels, cells in sorted(cell_cache.items()):
        best = None
        for cell in cells:
            for label in (*body, LEVEL, *labels):
                residue = cell * label % stream.L
                candidates = (
                    (14 * residue - stream.L, cell, label, "L"),
                    (13 * stream.L - 14 * (residue + label), cell, label, "R"),
                )
                for candidate in candidates:
                    if best is None or candidate < best:
                        best = candidate
        require(best is not None and best[0] >= 0, (body, labels, best))
        endpoints.append((labels, len(cells), best))

    records = tuple(records)
    carrier_records = tuple(carrier_records)
    endpoints = tuple(endpoints)
    hostile_sha = hashlib.sha256(repr(records).encode()).hexdigest()
    require(hostile_sha == EXPECTED_HOSTILE_RECORD_SHA256, hostile_sha)
    require(endpoints == (EXPECTED_ENDPOINT,), endpoints)
    require(
        (
            len(records),
            sum(row[-2] == "coarse-cardinality" for row in records),
            sum(row[-2] == "exact-cardinality" for row in records),
            min(row[-3] for row in records),
        )
        == (24, 23, 1, 6),
        "hostile terminal profile",
    )
    return records, carrier_records, endpoints, hostile_sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    tasks, neighbor_census, row_order, next_rows, next_sha = atlas_tasks()
    rows = run_screen(tasks, args.processes)
    screen_totals, audit_totals, screen_sha, residual_rows = screen_audit(rows)
    terminals, terminal_totals, terminal_sha = run_terminals(residual_rows)
    hostile_records, carriers, endpoints, hostile_sha = hostile_terminal_audit(
        residual_rows
    )
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z232-screen-terminal-v1",
        SOURCE_3078_SHA256,
        OUTPUT_3078_SHA256,
        SEMANTIC_3078_SHA256,
        SOURCE_3102_SHA256,
        OUTPUT_3102_SHA256,
        SEMANTIC_3102_SHA256,
        thm.ATLAS_SHA256,
        LEVEL,
        tuple(sorted(neighbor_census.items())),
        row_order,
        next_rows,
        rows,
        terminals,
        hostile_records,
        carriers,
        endpoints,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z232 exact screen and complete-cell cardinality descent",
        f"dependency=THM3078_source:{SOURCE_3078_SHA256};output:{OUTPUT_3078_SHA256};semantic:{SEMANTIC_3078_SHA256}",
        f"predecessor=THM3102_source:{SOURCE_3102_SHA256};output:{OUTPUT_3102_SHA256};semantic:{SEMANTIC_3102_SHA256};ledger:{LEDGER_BEFORE};cap:{LEVEL}",
        f"atlas=sha256:{thm.ATLAS_SHA256};rows:6060;neighbor_census:{tuple(sorted(neighbor_census.items()))}",
        f"layer=z1:{LEVEL};rows:{len(tasks)};wall:{NEIGHBOR_COUNTS[LEVEL][1]};order:{NEIGHBOR_COUNTS[LEVEL][2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256}",
        f"screen=states:{screen_totals[0]};crude:{screen_totals[1]};status:{screen_totals[2]};residual:{screen_totals[3]};residual_bodies:{len(residual_rows)};direct_farkas:{audit_totals[0]};legacy_farkas:{audit_totals[1]};screen_record_sha256:{screen_sha}",
        "order_control=rows:0;states:0;crude:0;status:0;residual:0;entire_layer_is_wall",
    ]
    for row in rows:
        lines.append(
            f"ROW;E={row[1]};L={row[2]};high={row[3]};branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.append(
        f"RESIDUAL;E={RESIDUAL_BODY};masks:{len(RESIDUAL_BANK)};sha256:{EXPECTED_RESIDUAL_SHA256};tuples:{RESIDUAL_BANK}"
    )
    lines.append(
        f"terminal=rows:{terminal_totals[0]};positive_two_high_gap:{terminal_totals[1]};closed:{terminal_totals[2]};zero_high_hostiles:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};coarse_cardinality:{terminal_totals[5]};exact_cardinality:{terminal_totals[6]};maxgap:{terminal_totals[7]};failures:{terminal_totals[8]};unit_checks:{terminal_totals[9]};terminal_record_sha256:{terminal_sha}"
    )
    terminal = terminals[0]
    lines.append(
        f"TERMINAL;E={terminal[1]};L={terminal[2]};high={terminal[3]};masks={terminal[4]};two_high_gap={ftext(terminal[5])};gap_witness={terminal[6]};zero_high_hostiles={terminal[7]};one_high_cases={terminal[8]};low_label_sets={terminal[9]};coarse_cardinality={terminal[10]};exact_cardinality={terminal[11]};maxgap={terminal[12]};failed={terminal[14]};minimum_certificate_slack={terminal[16]};case_sha256={terminal[17]}"
    )
    carrier = carriers[0]
    lines.append(
        f"hostile_terminal_audit=direct_full_grid_equals_scalar_equals_vector;cases:{len(hostile_records)};coarse:23;exact:1;minimum_actual_support_slack:6;record_sha256:{hostile_sha};carrier:E={carrier[0]},lows={carrier[1]},cells={carrier[2]},sha256={carrier[3]}"
    )
    for record in hostile_records:
        lines.append(
            f"CASE;ds={record[0]};high_d={record[1]};low_rows={record[2]};excess={record[3]};cells={record[4]};capacity={record[5]};coarse_lower_bound={record[6]};support={record[7]};kappa={record[8]};actual_slack={record[9]};method={record[10]};support_sha256={record[11]}"
        )
    lines.extend(
        [
            f"strict_open_control=E={RESIDUAL_BODY};lows={endpoints[0][0]};cells={endpoints[0][1]};minimum={endpoints[0][2]}",
            "direction_screen=ray_quotient_relaxes_global_suffix_order_but_preserves_fixed_first_label_distinct_label_witness_denominators_and_wall_high_gate;crude_and_exact_Farkas_status_exclusions_close_the_enlarged_state_set",
            "direction_terminal=positive_duplicate_permitting_two_high_gap_plus_wall_at_least_one_high_gate_forces_exactly_one_high;one_high_cases_uses_a_high_ray_supremum_and_is_a_superset_of_actual_literal_assignments",
            "translated_gate=complete_cell_support_uses_ceil(d/7),not_centered_beta;support_above_that_cap_gives_P_(E,Z)=T_for_every_high_label_on_the_denominator_ray;THM2941_completion_forces_P_subset_U_A_and_THM1166_gives_mu(U_A)<=36/91<1",
            f"zero_high_hostile=23_scalar_zero_high_passes_in_E={RESIDUAL_BODY}_are_excluded_by_the_inherited_wall_high_gate_and_are_not_counted_as_terminal_closures",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP};next_layer:z1={NEXT_LEVEL}_rows:{len(next_rows)}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{len(next_rows)};wall:{NEIGHBOR_COUNTS[NEXT_LEVEL][1]};order:{NEIGHBOR_COUNTS[NEXT_LEVEL][2]};row_order_sha256:{next_sha};status:occupied_unscouted_handoff",
        ]
    )
    for body, L, high, wall in next_rows:
        lines.append(
            f"NEXT;E={body};L={L};high={high};branch={'wall' if wall else 'order'}"
        )
    lines.extend(
        [
            "scope=projected_k3_necessary_atlas_only;z1=231_is_occupied_but_unscouted;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", newline="\n")
    print("\n".join(lines))


if __name__ == "__main__":
    mp.freeze_support()
    main()
