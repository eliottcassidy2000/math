#!/usr/bin/env python3
"""Exact projected-k3 z1=233 screen and complete-cell descent.

This companion reuses the promoted THM-3078 task-level evaluator rather than
changing a top-level frontier constant.  Every task carries ``first=233``;
the inherited evaluator sets both mutable FIRST locations before constructing
its exact ray quotient.  All 62 atlas rows are recomputed from scratch.

The four residual rows are then passed through THM-3078's physical terminal:
a duplicate-permitting two-high upper bound, exhaustive one-high low-pair
enumeration, and exact complete-cell translated-band cardinality.  A second
local audit compares every vectorized cell set with the scalar weak-endpoint
predicate and rebuilds every residue support.
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
SOURCE_3098 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.py"
)
OUTPUT_3098 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z233_exact_screen_complete_cell_cardinality_descent_thm3102.out"
)

SOURCE_3078_SHA256 = "2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168"
OUTPUT_3078_SHA256 = "d5fc52e922b7e083f5cbe5aa5a6066304c4e5c8c7963904dd3b649477caf5e42"
SEMANTIC_3078_SHA256 = "749a2292da58925c9653ce919d5bab6b374f64f8b713ac52efa1f18fba74c918"
SOURCE_3098_SHA256 = "744be96e68469afcf1864fadedfecfa237490b6ec3819f69bac16c52c5b558ed"
OUTPUT_3098_SHA256 = "722acc446f90a57df8ae84d50eb05d04250b1ffd688fb02e52ac0fcdd4111d07"
SEMANTIC_3098_SHA256 = "f5427dc8392ed287162f076df9b5e418afcec2a2f3a63a3b94a840c313a5b6bd"

LEVEL = 233
NEIGHBOR_COUNTS = {
    234: (381, 330, 51),
    233: (62, 45, 17),
    232: (3, 3, 0),
}
EXPECTED_ROW_ORDER_SHA256 = (
    "c503a3402cb2243ab586284d9d0773d07d0d0bddc95aae941c00ff7b2bf348ac"
)
EXPECTED_SCREEN = (1642, 570, 1039, 33)
EXPECTED_ORDER = (85, 53, 32, 0)
EXPECTED_AUDIT = (0, 1039)
EXPECTED_SCREEN_RECORD_SHA256 = (
    "2d24fe53a76095d32e7dfe0667ac1247fc01a6301a76f91dd8e5645be5c92313"
)
EXPECTED_TERMINAL = (4, 4, 4, 23, 33, 23, 10, 0, 0, 0)
EXPECTED_TERMINAL_RECORD_SHA256 = (
    "c1447cdbc2e1aabdf94f082a5bab5ac256480d25b6f236dd870f73240a989e0e"
)
EXPECTED_HOSTILE_RECORD_SHA256 = (
    "6db491bfaf2a70f0c12f0a690198cac8ab7a9adedfa70cb132606596ce6afa79"
)
EXPECTED_SEMANTIC_SHA256 = "6cf01affce52a5dcc67a8634da815780e3d357e72b295a7c4951211c6f12b0da"

LEDGER_BEFORE = 374387
LAYER_ROWS = 62
LEDGER_AFTER = 374325
NEXT_CAP = 232

RESIDUAL = {
    (1, 2, 5, 9, 12, 14): {
        "masks": 2,
        "mask_sha": "c9bfec7a7749796239d67be94074926cd2ac89e37515ea4753ac0e63357c4258",
        "gap": "21113341/6171345180",
        "zero": 0,
        "cases": 2,
        "coarse": 0,
        "exact": 2,
        "case_sha": "deb3af778b1d9d0103a013bf9a1a02bb64f10d88c4640c1d4272eb331663c7fe",
    },
    (1, 4, 5, 9, 12, 14): {
        "masks": 2,
        "mask_sha": "1532111bfb5cafba03c14ab3e116d1ed09a1fd194e453009fe03f4198311f2d9",
        "gap": "2583751/748041840",
        "zero": 0,
        "cases": 2,
        "coarse": 0,
        "exact": 2,
        "case_sha": "af4732c86fc10e42088c6f929c0024afd4438bbcd33c5735af39ee1ae8cd5cc3",
    },
    (1, 5, 6, 9, 12, 14): {
        "masks": 3,
        "mask_sha": "6203d6d71cedb53284763e70b2cf82fbdb17288eac6c3c5c6eff7ccbb9b40a6a",
        "gap": "233977/93059967",
        "zero": 0,
        "cases": 3,
        "coarse": 0,
        "exact": 3,
        "case_sha": "b07ec9fc8f3b7b00cd5f00f9cb4c2593661c9a78317c0b640c66f8e5a87526a4",
    },
    (1, 5, 9, 11, 12, 14): {
        "masks": 26,
        "mask_sha": "b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb",
        "gap": "39679901180/9741776933889",
        "zero": 23,
        "cases": 26,
        "coarse": 23,
        "exact": 3,
        "case_sha": "f6288fd9f7eceac090acee225c6a6d9f7e654f1bfa3d0ca1d80b5af62f3ef5d3",
    },
}

EXPECTED_ENDPOINTS = (
    ((1, 2, 5, 9, 12, 14), (234, 243), 3516, (0, 1350, 14, "L")),
    ((1, 4, 5, 9, 12, 14), (234, 243), 3388, (0, 1350, 14, "L")),
    ((1, 5, 6, 9, 12, 14), (234, 243), 3556, (0, 1350, 14, "L")),
    ((1, 5, 6, 9, 12, 14), (234, 324), 3648, (0, 1350, 14, "L")),
    ((1, 5, 6, 9, 12, 14), (234, 402), 3654, (0, 1350, 14, "L")),
    ((1, 5, 9, 11, 12, 14), (234, 243), 36550, (0, 14850, 14, "L")),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


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
require(sha(SOURCE_3098) == SOURCE_3098_SHA256, "THM-3098 source changed")
require(sha(OUTPUT_3098) == OUTPUT_3098_SHA256, "THM-3098 output changed")
require(
    f"semantic_sha256={SEMANTIC_3098_SHA256}" in OUTPUT_3098.read_text(),
    "THM-3098 semantic changed",
)
require(
    "promotion_consequence=ledger 374768-381=374387;"
    "projected_k3_cap:z1<=233;next_layer:z1=233_rows:62"
    in OUTPUT_3098.read_text(),
    "THM-3098 ledger consequence changed",
)

thm = load("thm3078_z233_descent_base", SOURCE_3078)
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
    tasks = []
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
        census[first][0] += 1
        census[first][1 if first < high else 2] += 1
        if first == LEVEL:
            tasks.append((first, body, L, high, first < high))
    derived = {level: tuple(census[level]) for level in NEIGHBOR_COUNTS}
    require(derived == NEIGHBOR_COUNTS, ("neighbor census", derived))
    require(len(tasks) == 62, len(tasks))
    order = tuple((task[1], task[2], task[3], task[4]) for task in tasks)
    require(
        hashlib.sha256(repr(order).encode()).hexdigest()
        == EXPECTED_ROW_ORDER_SHA256,
        "z233 row-order digest",
    )
    return tuple(tasks), derived, order


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
    order_rows = tuple(row for row in rows if not row[5])
    order_totals = tuple(
        sum(row[index] for row in order_rows) for index in (9, 10, 11, 12)
    )
    require(order_totals == EXPECTED_ORDER, order_totals)
    audit = (sum(row[19] for row in rows), sum(row[20] for row in rows))
    require(audit == EXPECTED_AUDIT, audit)
    require(all(row[16] == row[11] for row in rows), "unverified status row")
    record_sha = hashlib.sha256(repr(rows).encode()).hexdigest()
    require(record_sha == EXPECTED_SCREEN_RECORD_SHA256, record_sha)
    residual_rows = tuple(row for row in rows if row[12])
    require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), residual_rows)
    for row in residual_rows:
        expected = RESIDUAL[row[1]]
        require(row[12] == expected["masks"], (row[1], row[12]))
        require(
            hashlib.sha256(repr(row[13]).encode()).hexdigest()
            == expected["mask_sha"],
            (row[1], "mask digest"),
        )
    return totals, order_totals, audit, record_sha, residual_rows


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
    for row in terminals:
        expected = RESIDUAL[row[1]]
        require(row[4] == expected["masks"], (row[1], "terminal masks"))
        require(ftext(row[5]) == expected["gap"] and row[5] > 0,
                (row[1], "two-high gap", row[5]))
        require(
            (row[7], row[8], row[10], row[11], row[12], row[14], row[16], row[17])
            == (
                expected["zero"],
                expected["cases"],
                expected["coarse"],
                expected["exact"],
                0,
                0,
                1,
                expected["case_sha"],
            ),
            (row[1], "terminal profile"),
        )
        require(row[20], (row[1], "terminal not closed"))
    return terminals, totals, record_sha


def hostile_terminal_audit(residual_rows):
    all_records = []
    endpoint_records = []
    for screen_row in residual_rows:
        body = screen_row[1]
        residual = screen_row[13]
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(body)
        needed = {
            d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)
        }
        low, high, _signs, _checks = eng.build_literal_tables(stream, needed)
        cases = eng.one_high_cases(stream, residual, low, high)
        cell_cache = {}
        body_records = []
        for ds, high_d, low_rows, excess in cases:
            labels = tuple(sorted(label for _d, label in low_rows))
            if labels not in cell_cache:
                vector = tuple(
                    map(int, thm.vector_fixed_safe_cells(stream, labels).tolist())
                )
                scalar = eng.fixed_safe_cells(stream, labels)
                require(vector == scalar, (body, labels, "vector/scalar cells"))
                cell_cache[labels] = scalar
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
            record = (
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
            body_records.append(record)
            all_records.append((body, record))

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
            endpoint_records.append((body, labels, len(cells), best))

    all_records = tuple(all_records)
    endpoint_records = tuple(endpoint_records)
    hostile_sha = hashlib.sha256(repr(all_records).encode()).hexdigest()
    require(hostile_sha == EXPECTED_HOSTILE_RECORD_SHA256, hostile_sha)
    require(endpoint_records == EXPECTED_ENDPOINTS, endpoint_records)
    require(
        (
            len(all_records),
            sum(row[-2] == "coarse-cardinality" for _body, row in all_records),
            sum(row[-2] == "exact-cardinality" for _body, row in all_records),
            min(row[-3] for _body, row in all_records),
        )
        == (33, 23, 10, 1),
        "hostile terminal profile",
    )
    return all_records, endpoint_records, hostile_sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    tasks, neighbor_census, row_order = atlas_tasks()
    rows = run_screen(tasks, args.processes)
    screen_totals, order_totals, audit_totals, screen_sha, residual_rows = (
        screen_audit(rows)
    )
    terminals, terminal_totals, terminal_sha = run_terminals(residual_rows)
    hostile_records, endpoints, hostile_sha = hostile_terminal_audit(residual_rows)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z233-screen-terminal-v1",
        SOURCE_3078_SHA256,
        OUTPUT_3078_SHA256,
        SEMANTIC_3078_SHA256,
        SOURCE_3098_SHA256,
        OUTPUT_3098_SHA256,
        SEMANTIC_3098_SHA256,
        thm.ATLAS_SHA256,
        LEVEL,
        tuple(sorted(neighbor_census.items())),
        row_order,
        rows,
        terminals,
        hostile_records,
        endpoints,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z233 exact screen and complete-cell cardinality descent",
        f"dependency=THM3078_source:{SOURCE_3078_SHA256};output:{OUTPUT_3078_SHA256};semantic:{SEMANTIC_3078_SHA256}",
        f"predecessor=THM3098_source:{SOURCE_3098_SHA256};output:{OUTPUT_3098_SHA256};semantic:{SEMANTIC_3098_SHA256};ledger:{LEDGER_BEFORE};cap:{LEVEL}",
        f"atlas=sha256:{thm.ATLAS_SHA256};rows:6060;neighbor_census:{tuple(sorted(neighbor_census.items()))}",
        f"layer=z1:{LEVEL};rows:{len(tasks)};wall:{NEIGHBOR_COUNTS[LEVEL][1]};order:{NEIGHBOR_COUNTS[LEVEL][2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256}",
        f"screen=states:{screen_totals[0]};crude:{screen_totals[1]};status:{screen_totals[2]};residual:{screen_totals[3]};residual_bodies:{len(residual_rows)};direct_farkas:{audit_totals[0]};legacy_farkas:{audit_totals[1]};screen_record_sha256:{screen_sha}",
        f"order_control=rows:{NEIGHBOR_COUNTS[LEVEL][2]};states:{order_totals[0]};crude:{order_totals[1]};status:{order_totals[2]};residual:{order_totals[3]};order_rows_never_enter_wall_terminal",
    ]
    for row in rows:
        lines.append(
            f"ROW;E={row[1]};L={row[2]};high={row[3]};branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.append(
        f"terminal=rows:{terminal_totals[0]};positive_two_high_gap:{terminal_totals[1]};closed:{terminal_totals[2]};zero_high_hostiles:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};coarse_cardinality:{terminal_totals[5]};exact_cardinality:{terminal_totals[6]};maxgap:{terminal_totals[7]};failures:{terminal_totals[8]};unit_checks:{terminal_totals[9]};terminal_record_sha256:{terminal_sha}"
    )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};gap_witness={row[6]};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse_cardinality={row[10]};exact_cardinality={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_sha256={row[17]}"
        )
    lines.extend(
        [
            f"hostile_terminal_audit=scalar_vector_cell_sets_equal;cases:{len(hostile_records)};coarse:23;exact:10;minimum_support_slack:1;record_sha256:{hostile_sha}",
            "strict_open_controls=" + ";".join(
                f"E={body},lows={labels},cells={count},minimum={minimum}"
                for body, labels, count, minimum in endpoints
            ),
            "direction_screen=ray_quotient_relaxes_global_suffix_order_but_preserves_fixed_first_label_distinct_label_witness_denominators_and_wall_high_gate;crude_and_exact_Farkas_status_exclusions_close_the_enlarged_state_set",
            "direction_terminal=positive_duplicate_permitting_two_high_gap_plus_wall_at_least_one_high_gate_forces_exactly_one_high;one_high_cases_uses_a_high_ray_supremum_and_is_a_superset_of_actual_literal_assignments",
            "translated_gate=complete_cell_support_uses_ceil(d/7),not_centered_beta;support_above_that_cap_gives_P_(E,Z)=T_for_every_high_label_on_the_denominator_ray;THM2941_completion_forces_P_subset_U_A_and_THM1166_gives_mu(U_A)<=36/91<1",
            "zero_high_hostile=23_scalar_zero_high_passes_in_E=(1,5,9,11,12,14)_are_excluded_by_the_inherited_wall_high_gate_and_are_not_counted_as_terminal_closures",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP};next_layer:z1=232_rows:3",
            "scope=projected_k3_necessary_atlas_only;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", newline="\n")
    print("\n".join(lines))


if __name__ == "__main__":
    mp.freeze_support()
    main()
