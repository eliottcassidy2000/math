#!/usr/bin/env python3
"""Exact projected-k3 z1=226 screen and terminal descent addendum."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3114 = ROOT / (
    "04-computation/lrc14_j7_k3_z227_exact_screen_zero_residual_descent_thm3114.py"
)
OUTPUT_3114 = ROOT / (
    "05-knowledge/results/lrc14_j7_k3_z227_exact_screen_zero_residual_descent_thm3114.out"
)
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z226_terminal_descent_thm3114.out"

SOURCE_3114_SHA256 = "a05dc531b8e32201e0a5dc3ef92bd0f89ecc5c7f25723275dce60d1f98703840"
OUTPUT_3114_SHA256 = "424bf01ca72a3d5a222eb2411ae58a718bde127109ff0d64182f75342b0c2a45"
SEMANTIC_3114_SHA256 = "17347efff69b101d47afc7d7097fcde498d210ed25737c1ef9935d6b3a5baacb"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 226
NEXT_LEVEL = 225
EXPECTED_CENSUS = (13, 13, 0)
EXPECTED_ROW_ORDER_SHA256 = "ec2f1218d61dd69d3811d68eb87a23254d276bb75647be2d4c883affe53a520e"
EXPECTED_SCREEN = (979, 526, 363, 90)
EXPECTED_ORDER = (0, 0, 0, 0)
EXPECTED_FARKAS = (0, 363)
EXPECTED_SCREEN_SHA256 = "3ff3c7ae650ed5e5f1ba67d49fa4cab1315a268498f9c4bd36d00443dcf8224e"
EXPECTED_NEXT_CENSUS = (78, 78, 0)
EXPECTED_NEXT_ROW_ORDER_SHA256 = "9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719"

RESIDUAL = {
    (1, 5, 6, 9, 12, 14): {
        "masks": 1,
        "mask_sha": "85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a",
        "gap": "2436193/902641740",
        "zero": 0,
        "cases": 1,
        "coarse": 0,
        "exact": 1,
    },
    (1, 5, 9, 11, 12, 14): {
        "masks": 65,
        "mask_sha": "997cc17a4be38886e37970060f6c0d92ee6904b345e973083f9b1b52096ad882",
        "gap": "15491357285/4724552761929",
        "zero": 61,
        "cases": 65,
        "coarse": 61,
        "exact": 4,
    },
    (1, 9, 10, 11, 12, 14): {
        "masks": 24,
        "mask_sha": "fe196c407eb6d7d1fc5a65f283a8ce078e534602152521e8b1081d908ed725ea",
        "gap": "1545772719847/373939773753546",
        "zero": 23,
        "cases": 24,
        "coarse": 23,
        "exact": 1,
    },
}

EXPECTED_TERMINAL = (3, 3, 3, 84, 90, 84, 6, 0, 0, 0)
EXPECTED_TERMINAL_SHA256 = "f343c72ad5b49b37eca86171374b274be6ef48bd9e28115107daf9ccb140c5a3"
EXPECTED_CARRIER_SHA256 = "f5060567a3cc0341902c541ec6910161d509b48d16f5303a467095c40fbab66f"
EXPECTED_CASE_SHA256 = "6eaf93061816e63631a0e5180697cc1be7fe69980705d2d87f0a4a822698cfa9"

LEDGER_BEFORE = 374185
LAYER_ROWS = 13
LEDGER_AFTER = 374172
NEXT_CAP = 225
EXPECTED_SEMANTIC_SHA256 = "cdd4f83bbd3a3244244bf14baad8935a43a4a33c0bc973acff5dcc4048a4dc52"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    payload = path.read_bytes()
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_layers(base) -> tuple[tuple, tuple]:
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    selected = {LEVEL: [], NEXT_LEVEL: []}
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
    rows = tuple(selected[LEVEL])
    next_rows = tuple(selected[NEXT_LEVEL])
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == EXPECTED_NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256,
            "row order")
    require(hashlib.sha256(repr(next_rows).encode()).hexdigest()
            == EXPECTED_NEXT_ROW_ORDER_SHA256, "next row order")
    return rows, next_rows


def screen_worker(task):
    prior = load(f"thm3114_z226_worker_{task[1]}", SOURCE_3114)
    return prior.screen_worker(task)


def run_screen(rows: tuple, processes: int) -> tuple:
    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost or duplicated task")
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    order = tuple(row for row in screened if not row[5])
    order_totals = tuple(sum(row[index] for row in order) for index in (9, 10, 11, 12))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(totals == EXPECTED_SCREEN, totals)
    require(order_totals == EXPECTED_ORDER, order_totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "status verification")
    screen_sha = hashlib.sha256(repr(screened).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual = tuple(row for row in screened if row[12])
    require(tuple(row[1] for row in residual) == tuple(RESIDUAL), "residual bodies")
    for row in residual:
        expected = RESIDUAL[row[1]]
        require(row[12] == expected["masks"], (row[1], "mask count"))
        require(hashlib.sha256(repr(row[13]).encode()).hexdigest()
                == expected["mask_sha"], (row[1], "mask digest"))
    return screened, totals, order_totals, farkas, screen_sha, residual


def run_terminals(base, residual: tuple) -> tuple:
    terminals = tuple(
        base.thm.terminal_probe((LEVEL, row[1], row[13])) for row in residual
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
        require(row[5] > 0 and ftext(row[5]) == expected["gap"], (row[1], row[5]))
        require((row[7], row[8], row[9], row[10], row[11], row[12], row[14], row[16], row[20])
                == (expected["zero"], expected["cases"], 1, expected["coarse"],
                    expected["exact"], 0, 0, 1, True),
                (row[1], "terminal profile"))
    return terminals, totals, record_sha


def direct_carrier_audit(predecessor, base, residual: tuple) -> tuple:
    carriers = []
    records = []
    for screen_row in residual:
        body = screen_row[1]
        masks = screen_row[13]
        eng = base.eng
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(body)
        needed = {d for ds in masks for d in eng.suffix_slots(ds, stream.first_d)}
        low, high, _signs, _checks = eng.build_literal_tables(stream, needed)
        cases = eng.one_high_cases(stream, masks, low, high)
        require(len(cases) == RESIDUAL[body]["cases"], (body, "case count"))
        cache = {}
        for _ds, _high_d, low_rows, _excess in cases:
            labels = tuple(sorted(label for _d, label in low_rows))
            if labels in cache:
                continue
            direct = predecessor.direct_cells(stream.L, (*body, LEVEL, *labels))
            scalar = eng.fixed_safe_cells(stream, labels)
            vector = tuple(map(int, base.thm.vector_fixed_safe_cells(stream, labels).tolist()))
            require(direct == scalar == vector, (body, labels, "carrier mismatch"))
            endpoint = None
            for cell in direct:
                for label in (*body, LEVEL, *labels):
                    residue = cell * label % stream.L
                    for candidate in (
                        (14 * residue - stream.L, cell, label, "L"),
                        (13 * stream.L - 14 * (residue + label), cell, label, "R"),
                    ):
                        if endpoint is None or candidate < endpoint:
                            endpoint = candidate
            require(endpoint is not None and endpoint[0] >= 0,
                    (body, labels, endpoint))
            cache[labels] = direct
            carriers.append((body, labels, len(direct),
                             hashlib.sha256(repr(direct).encode()).hexdigest(), endpoint))

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
            records.append((
                body, ds, high_d, low_rows, excess, labels, len(cells), capacity,
                coarse, len(support), kappa, len(support) - kappa, method,
                hashlib.sha256(repr(support).encode()).hexdigest(),
            ))
    carriers = tuple(carriers)
    records = tuple(records)
    require((len(carriers), len(records)) == (3, 90), "carrier/case count")
    require(sum(row[-2] == "coarse-cardinality" for row in records) == 84,
            "coarse count")
    require(sum(row[-2] == "exact-cardinality" for row in records) == 6,
            "exact count")
    require(min(row[8] - row[10] for row in records
                if row[-2] == "coarse-cardinality") == 1, "coarse margin")
    require(min(row[-3] for row in records) == 1, "support slack")
    carrier_sha = hashlib.sha256(repr(carriers).encode()).hexdigest()
    case_sha = hashlib.sha256(repr(records).encode()).hexdigest()
    require(carrier_sha == EXPECTED_CARRIER_SHA256, carrier_sha)
    require(case_sha == EXPECTED_CASE_SHA256, case_sha)
    return carriers, records, carrier_sha, case_sha


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3114) == SOURCE_3114_SHA256, "THM-3114 source changed")
    require(sha(OUTPUT_3114) == OUTPUT_3114_SHA256, "THM-3114 output changed")
    require(f"semantic_sha256={SEMANTIC_3114_SHA256}" in OUTPUT_3114.read_text(),
            "THM-3114 semantic changed")

    prior = load("thm3114_z226_prior", SOURCE_3114)
    screen_source = prior.load("thm3114_z226_screen_source", prior.SOURCE_3113)
    driver = screen_source.load("thm3114_z226_driver")
    predecessor = driver.load("thm3114_z226_predecessor")
    base = predecessor.load("thm3114_z226_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, "atlas semantic changed")
    rows, next_rows = parse_layers(base)
    screened, totals, order_totals, farkas, screen_sha, residual = run_screen(
        rows, args.processes
    )
    terminals, terminal_totals, terminal_sha = run_terminals(base, residual)
    carriers, cases, carrier_sha, case_sha = direct_carrier_audit(
        predecessor, base, residual
    )
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z226-terminal-descent-v1",
        SOURCE_3114_SHA256,
        OUTPUT_3114_SHA256,
        SEMANTIC_3114_SHA256,
        ATLAS_SHA256,
        rows,
        screened,
        terminals,
        carriers,
        cases,
        next_rows,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z226 terminal descent addendum to THM-3114",
        f"dependency=THM3114_source:{SOURCE_3114_SHA256};output:{OUTPUT_3114_SHA256};semantic:{SEMANTIC_3114_SHA256}",
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{EXPECTED_CENSUS[0]};wall:{EXPECTED_CENSUS[1]};order:{EXPECTED_CENSUS[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};order_control:{order_totals}",
    ]
    for row in screened:
        lines.append(
            f"ROW;z1={LEVEL};E={row[1]};L={row[2]};high={row[3]};"
            f"branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};"
            f"status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};"
            f"stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.append(
        f"terminal=z1:{LEVEL};rows:{terminal_totals[0]};positive_two_high_gap:{terminal_totals[1]};closed:{terminal_totals[2]};zero_high_hostiles:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};coarse:{terminal_totals[5]};exact:{terminal_totals[6]};maxgap:{terminal_totals[7]};failures:{terminal_totals[8]};unit_checks:{terminal_totals[9]};terminal_record_sha256:{terminal_sha}"
    )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};"
            f"two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};"
            f"one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};"
            f"exact={row[11]};failed={row[14]};minimum_certificate_slack={row[16]};"
            f"case_certificate_sha256={row[17]}"
        )
    lines.append(
        f"direct_carrier_audit=carriers:{len(carriers)};cases:{len(cases)};coarse:84;exact:6;minimum_coarse_margin:1;minimum_support_slack:1;all_exact_supports:full;carrier_record_sha256:{carrier_sha};case_record_sha256:{case_sha}"
    )
    for body, labels, count, cells_sha, endpoint in carriers:
        lines.append(
            f"CARRIER;E={body};lows={labels};cells={count};cells_sha256={cells_sha};"
            f"minimum_weak_endpoint={endpoint}"
        )
    lines.extend([
        "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment",
        "direction_terminal=each_positive_duplicate_permitting_two_high_gap_together_with_the_wall_at_least_one_high_gate_forces_exactly_one_high;the_high_ray_bank_enlarges_the_actual_assignment_set",
        "direction_carrier=each_direct_complete_cell_is_an_inner_carrier_in_the_strict_open_safe_set;coarse_or_exact_support_above_ceil(d/7)_forces_the_completed_carrier_contradiction",
        f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
        f"next_layer=z1:{NEXT_LEVEL};rows:{EXPECTED_NEXT_CENSUS[0]};wall:{EXPECTED_NEXT_CENSUS[1]};order:{EXPECTED_NEXT_CENSUS[2]};row_order_sha256:{EXPECTED_NEXT_ROW_ORDER_SHA256};status:occupied_unscouted_handoff",
        "scope=projected_k3_necessary_atlas_only;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ])
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
