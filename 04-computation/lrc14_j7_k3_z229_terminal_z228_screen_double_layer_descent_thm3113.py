#!/usr/bin/env python3
"""Exact projected-k3 z1=229 terminal and z1=228 screen descent.

The complete z1=229 and z1=228 atlas layers are recomputed through the
promoted THM-3111/THM-3109 screen.  Six z1=229 residual bodies are closed by
the inherited positive duplicate-two-high gap and translated complete-cell
cardinality terminal.  The five z1=228 rows have no residual state at all.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3111 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py"
)
OUTPUT_3111 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.out"
)

SOURCE_3111_SHA256 = "d86fcbf53f04dcc9e08d33f48c0682b063473a2d9d422a75c3ef37b58e9e1841"
OUTPUT_3111_SHA256 = "8a7b5a82e30cb84379655e50109803ea973c1df69a5affbcc7f29c0aeb76ee7e"
SEMANTIC_3111_SHA256 = "b4dd4e0e66d0429cac12dc96d1d769d3f05ab3c503e1dcbb1b94f5e6020c9b6c"

LEVELS = (229, 228)
NEXT_LEVEL = 227
EXPECTED_CENSUS = {229: (43, 36, 7), 228: (5, 5, 0), 227: (30, 28, 2)}
EXPECTED_ROW_ORDER_SHA256 = {
    229: "7449dd7ad70cf3c76c32edb2cc509e29989ac008c2e9a9968bceaabf3e2a2587",
    228: "ba3819e24debb3708459e0293ca4f35127ac533b1a7491a30212fb3fd9242cda",
    227: "17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c",
}
EXPECTED_SCREEN = {229: (1195, 571, 512, 112), 228: (99, 94, 5, 0)}
EXPECTED_ORDER = {229: (24, 24, 0, 0), 228: (0, 0, 0, 0)}
EXPECTED_FARKAS = {229: (0, 512), 228: (0, 5)}
EXPECTED_SCREEN_SHA256 = {
    229: "a3b822b2ba2413b498dd02b787b41171be51942e71e0392a7cbd53a0ff188059",
    228: "71cbcedb29695778b51252c711d9dbfc6ed83e5f886ed382c1fc9c8956e22928",
}
EXPECTED_TERMINAL = (6, 6, 6, 100, 112, 99, 13, 0, 0, 0)
EXPECTED_TERMINAL_SHA256 = "516b596b690bcb176afbbe02fee672171f2a8bbf359ba031d1983819ad87a9d0"
EXPECTED_CARRIER_SHA256 = "b6ed1106b0061b86909e11559c11244f267fd5112acfe8f550b978f4dccc6c51"
EXPECTED_DIRECT_SHA256 = "e7f794f8ce3ffbfbc014899b7fc7e7b22158ff7eb56d8ad02c93c6b0d7659e33"
EXPECTED_SEMANTIC_SHA256 = "a024972dce2c159096db81a751b49654ebc8fad18e2fb1229dbfa0f01784e7e4"

LEDGER_BEFORE = 374263
LAYER_ROWS = 48
LEDGER_AFTER = 374215
NEXT_CAP = 227

RESIDUAL = {
    (1, 4, 5, 9, 11, 14): {
        "masks": 14,
        "mask_sha": "0612bfd6334b896080174e3f633966df7c0ba1610b8d4a165acdda240b219959",
        "gap": "1346009063111/373226766637725",
        "zero": 11,
        "cases": 14,
        "coarse": 11,
        "exact": 3,
        "minimum": 1,
    },
    (1, 5, 6, 8, 9, 14): {
        "masks": 1,
        "mask_sha": "d7fc8f8c673d890f158cd1d55b10718afe3000a727c619aeeb130962b91ad43b",
        "gap": "60679854023/18484273855320",
        "zero": 1,
        "cases": 1,
        "coarse": 1,
        "exact": 0,
        "minimum": 127,
    },
    (1, 5, 9, 11, 12, 14): {
        "masks": 26,
        "mask_sha": "b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb",
        "gap": "145522417903/38298144512628",
        "zero": 23,
        "cases": 26,
        "coarse": 23,
        "exact": 3,
        "minimum": 1,
    },
    (1, 8, 10, 11, 12, 14): {
        "masks": 19,
        "mask_sha": "c284a76f62a2d6dbdea4b0f95ee3a6162b66ac1fe72cff43278cbe755ebb59ee",
        "gap": "14971895384/3550561219205",
        "zero": 16,
        "cases": 19,
        "coarse": 15,
        "exact": 4,
        "minimum": 1,
    },
    (2, 5, 8, 9, 11, 14): {
        "masks": 42,
        "mask_sha": "7f2a84a6ca43c97561c7d38f9773b09c6c4cc7be3903e7a1721293c2508d83ad",
        "gap": "26907523189147/7661147274160380",
        "zero": 39,
        "cases": 42,
        "coarse": 39,
        "exact": 3,
        "minimum": 1,
    },
    (2, 8, 10, 11, 12, 14): {
        "masks": 10,
        "mask_sha": "00415b63ca4ba804c327c3ebb302ddbaec054ca527ca8eb5195e8c66adc8c5ac",
        "gap": "396680527987/91537658972760",
        "zero": 10,
        "cases": 10,
        "coarse": 10,
        "exact": 0,
        "minimum": 14,
    },
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3111):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3111_worker_{task[0]}_{task[1]}").screen_worker(task)


def atlas_layers(base):
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    selected = {level: [] for level in (*LEVELS, NEXT_LEVEL)}
    for line in base.thm.ATLAS.read_text().splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        level = int(match.group(4))
        if level not in selected:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        high = int(match.group(3))
        selected[level].append((body, int(match.group(2)), high, level < high))
    answer = {}
    for level, rows in selected.items():
        rows = tuple(rows)
        census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
        require(census == EXPECTED_CENSUS[level], (level, census))
        require(
            hashlib.sha256(repr(rows).encode()).hexdigest()
            == EXPECTED_ROW_ORDER_SHA256[level],
            (level, rows),
        )
        answer[level] = rows
    return answer


def run_screen(level, rows, processes):
    tasks = tuple((level, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(processes=min(processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), (level, "lost or duplicated task"))
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN[level], (level, totals))
    order = tuple(row for row in screened if not row[5])
    order_totals = tuple(sum(row[index] for row in order) for index in (9, 10, 11, 12))
    require(order_totals == EXPECTED_ORDER[level], (level, order_totals))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS[level], (level, farkas))
    require(all(row[16] == row[11] for row in screened), (level, "status verification"))
    record_sha = hashlib.sha256(repr(screened).encode()).hexdigest()
    require(record_sha == EXPECTED_SCREEN_SHA256[level], (level, record_sha))
    residual_rows = tuple(row for row in screened if row[12])
    if level == 228:
        require(not residual_rows, residual_rows)
    else:
        require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), residual_rows)
        require(all(row[5] for row in residual_rows), "non-wall residual")
        for row in residual_rows:
            expected = RESIDUAL[row[1]]
            require(row[12] == expected["masks"], (row[1], row[12]))
            require(
                hashlib.sha256(repr(row[13]).encode()).hexdigest()
                == expected["mask_sha"],
                (row[1], "mask digest"),
            )
    return screened, totals, order_totals, farkas, record_sha, residual_rows


def run_terminals(base, residual_rows):
    terminals = tuple(
        base.thm.terminal_probe((229, row[1], row[13])) for row in residual_rows
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
    require(record_sha == EXPECTED_TERMINAL_SHA256, record_sha)
    for row in terminals:
        expected = RESIDUAL[row[1]]
        require(row[4] == expected["masks"], (row[1], "terminal masks"))
        require(ftext(row[5]) == expected["gap"] and row[5] > 0, (row[1], row[5]))
        require(
            (row[7], row[8], row[9], row[10], row[11], row[12], row[14], row[16], row[20])
            == (
                expected["zero"], expected["cases"], 1, expected["coarse"],
                expected["exact"], 0, 0, expected["minimum"], True,
            ),
            (row[1], "terminal profile"),
        )
    return terminals, totals, record_sha


def direct_carrier_audit(predecessor, base, residual_rows):
    carrier_records = []
    case_records = []
    for screen_row in residual_rows:
        body = screen_row[1]
        residual = screen_row[13]
        eng = base.eng
        eng.FIRST = 229
        eng.ray.FIRST = 229
        stream = eng.ray.Stream(body)
        needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
        low, high, _signs, _checks = eng.build_literal_tables(stream, needed)
        cases = eng.one_high_cases(stream, residual, low, high)
        require(len(cases) == RESIDUAL[body]["cases"], (body, "case count"))
        cache = {}
        for _ds, _high_d, low_rows, _excess in cases:
            labels = tuple(sorted(label for _d, label in low_rows))
            if labels in cache:
                continue
            direct = predecessor.direct_cells(stream.L, (*body, 229, *labels))
            scalar = eng.fixed_safe_cells(stream, labels)
            vector = tuple(map(int, base.thm.vector_fixed_safe_cells(stream, labels).tolist()))
            require(direct == scalar == vector, (body, labels, "carrier mismatch"))
            endpoint = None
            for cell in direct:
                for label in (*body, 229, *labels):
                    residue = cell * label % stream.L
                    for candidate in (
                        (14 * residue - stream.L, cell, label, "L"),
                        (13 * stream.L - 14 * (residue + label), cell, label, "R"),
                    ):
                        if endpoint is None or candidate < endpoint:
                            endpoint = candidate
            require(endpoint is not None and endpoint[0] >= 0, (body, labels, endpoint))
            cache[labels] = direct
            carrier_records.append(
                (body, labels, len(direct), hashlib.sha256(repr(direct).encode()).hexdigest(), endpoint)
            )

        for ds, high_d, low_rows, excess in cases:
            labels = tuple(sorted(label for _d, label in low_rows))
            cells = cache[labels]
            require(stream.L % high_d == 0, (body, high_d))
            capacity = stream.L // high_d
            coarse = (len(cells) + capacity - 1) // capacity
            support = tuple(sorted({cell % high_d for cell in cells}))
            kappa = (high_d + 6) // 7
            method = "coarse-cardinality" if coarse > kappa else "exact-cardinality"
            require(coarse > kappa or len(support) > kappa, (body, ds, high_d))
            if method == "exact-cardinality":
                require(support == tuple(range(high_d)), (body, high_d, support))
            case_records.append(
                (
                    body, ds, high_d, low_rows, excess, labels, len(cells), capacity,
                    coarse, len(support), kappa, len(support) - kappa, method,
                    hashlib.sha256(repr(support).encode()).hexdigest(),
                )
            )

    carrier_records = tuple(carrier_records)
    case_records = tuple(case_records)
    require(len(carrier_records) == 6, len(carrier_records))
    require(len(case_records) == 112, len(case_records))
    require(sum(row[-2] == "coarse-cardinality" for row in case_records) == 99, "coarse count")
    require(sum(row[-2] == "exact-cardinality" for row in case_records) == 13, "exact count")
    require(min(row[8] - row[10] for row in case_records if row[-2] == "coarse-cardinality") == 1, "coarse margin")
    require(min(row[-3] for row in case_records) == 1, "support slack")
    carrier_sha = hashlib.sha256(repr(carrier_records).encode()).hexdigest()
    direct_sha = hashlib.sha256(repr(case_records).encode()).hexdigest()
    require(carrier_sha == EXPECTED_CARRIER_SHA256, carrier_sha)
    require(direct_sha == EXPECTED_DIRECT_SHA256, direct_sha)
    return carrier_records, case_records, carrier_sha, direct_sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3111) == SOURCE_3111_SHA256, "THM-3111 source changed")
    require(sha(OUTPUT_3111) == OUTPUT_3111_SHA256, "THM-3111 output changed")
    require(f"semantic_sha256={SEMANTIC_3111_SHA256}" in OUTPUT_3111.read_text(), "THM-3111 semantic changed")

    driver = load("thm3113_driver")
    predecessor = driver.load("thm3113_predecessor")
    base = predecessor.load("thm3113_base")
    layers = atlas_layers(base)
    screens = {level: run_screen(level, layers[level], args.processes) for level in LEVELS}
    terminals, terminal_totals, terminal_sha = run_terminals(base, screens[229][5])
    carriers, cases, carrier_sha, direct_sha = direct_carrier_audit(predecessor, base, screens[229][5])
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z229-terminal-z228-screen-double-descent-v1",
        SOURCE_3111_SHA256, OUTPUT_3111_SHA256, SEMANTIC_3111_SHA256,
        base.thm.ATLAS_SHA256, layers, tuple(screens[level][0] for level in LEVELS),
        terminals, carriers, cases,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z229 terminal and z228 screen double-layer descent",
        f"dependency=THM3111_source:{SOURCE_3111_SHA256};output:{OUTPUT_3111_SHA256};semantic:{SEMANTIC_3111_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
    ]
    for level in LEVELS:
        screened, totals, order_totals, farkas, screen_sha, _residual = screens[level]
        census = EXPECTED_CENSUS[level]
        lines.append(
            f"layer=z1:{level};rows:{census[0]};wall:{census[1]};order:{census[2]};"
            f"row_order_sha256:{EXPECTED_ROW_ORDER_SHA256[level]};states:{totals[0]};"
            f"crude:{totals[1]};status:{totals[2]};residual:{totals[3]};"
            f"direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};"
            f"order_control:{order_totals}"
        )
        for row in screened:
            lines.append(
                f"ROW;z1={level};E={row[1]};L={row[2]};high={row[3]};"
                f"branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};"
                f"status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};"
                f"stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
            )
    lines.append(
        f"terminal=z1:229;rows:{terminal_totals[0]};positive_two_high_gap:{terminal_totals[1]};"
        f"closed:{terminal_totals[2]};zero_high_hostiles:{terminal_totals[3]};"
        f"one_high_cases:{terminal_totals[4]};coarse:{terminal_totals[5]};exact:{terminal_totals[6]};"
        f"maxgap:{terminal_totals[7]};failures:{terminal_totals[8]};unit_checks:{terminal_totals[9]};"
        f"terminal_record_sha256:{terminal_sha}"
    )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};"
            f"two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};"
            f"low_label_sets={row[9]};coarse={row[10]};exact={row[11]};failed={row[14]};"
            f"minimum_certificate_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.append(
        f"direct_carrier_audit=carriers:{len(carriers)};cases:{len(cases)};coarse:99;exact:13;"
        f"minimum_coarse_margin:1;minimum_support_slack:1;all_exact_supports:full;"
        f"carrier_record_sha256:{carrier_sha};case_record_sha256:{direct_sha}"
    )
    for body, labels, count, cells_sha, endpoint in carriers:
        lines.append(
            f"CARRIER;E={body};lows={labels};cells={count};cells_sha256={cells_sha};"
            f"minimum_weak_endpoint={endpoint}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment",
            "direction_terminal=at_z229_each_strictly_positive_duplicate_permitting_two_high_gap_together_with_the_wall_at_least_one_high_gate_forces_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=each_direct_complete_cell_is_an_inner_carrier_wholly_contained_in_the_strict_open_safe_set;the_coarse_bound_or_exact_full_projected_support_above_ceil(d/7)_forces_the_completed_carrier_contradiction",
            "transition=z229_requires_the_terminal_but_z228_has_no_residual_mask_and_dies_at_the_screen",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP};next_layer:z1={NEXT_LEVEL}_rows:{EXPECTED_CENSUS[NEXT_LEVEL][0]}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{EXPECTED_CENSUS[NEXT_LEVEL][0]};wall:{EXPECTED_CENSUS[NEXT_LEVEL][1]};order:{EXPECTED_CENSUS[NEXT_LEVEL][2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256[NEXT_LEVEL]};status:occupied_unscouted_handoff",
            "scope=projected_k3_necessary_atlas_only;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
