#!/usr/bin/env python3
"""Exact projected-k3 z1=231 screen and complete-cell descent.

The nine atlas rows handed off by THM-3106 are recomputed with the promoted
THM-3078 screen.  The sole two-mask residue is then audited on the complete
Z/360360Z carrier, independently of the inherited vector implementation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3106 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py"
)
OUTPUT_3106 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.out"
)

SOURCE_3106_SHA256 = "9c38e808e22c9ac376217b9c76da69f198a14a4060cd4af4bf3e10b2c6a604f6"
OUTPUT_3106_SHA256 = "286bb3e31e1ef28af5640a8ebe5e0df5e57b3e7253aa4c1e5ac0483beb9d0e63"
SEMANTIC_3106_SHA256 = "a14adcd1e52323baf2b791f55d5846c4fbd422e3441fa2411138bdd98acca3d3"

LEVEL = 231
NEXT_LEVEL = 230
EXPECTED_ROW_ORDER_SHA256 = "ff12b3d3c41a10aa38da0a1483bfb506b4c6a9a771c2325fb4df1df867088f5e"
EXPECTED_NEXT_ROW_ORDER_SHA256 = "f5129616fe3d7f5884350abe7e77bb2684798e213600e266cd456bc9c582f058"
EXPECTED_LAYER_CENSUS = (9, 5, 4)
EXPECTED_NEXT_CENSUS = (50, 48, 2)
EXPECTED_SCREEN = (226, 127, 97, 2)
EXPECTED_FARKAS = (0, 97)
EXPECTED_SCREEN_SHA256 = "80cc1065fb84319d2d063b911cfe24576e088596ecf815a98d50e38a6f32a685"

RESIDUAL_BODY = (1, 5, 9, 11, 12, 13)
RESIDUAL_BANK = (
    (1560, 3080, 10296, 40040),
    (1560, 3080, 40040, 51480),
)
EXPECTED_RESIDUAL_SHA256 = "89963fe483156389680ed3a1bf1c6feb92fdc00f72ce6cffaccb3d6ca1e88cc8"
EXPECTED_RESIDUAL_RECORD_SHA256 = "b942954069fee35af4a5afb60815c61833ecb87292f2a64b4b26da9f226e7559"
EXPECTED_TERMINAL_SHA256 = "3ed64c990693b39b5064f7ceb6fcefd95ad1b2eea8c13bc159e6a79a7f130899"
EXPECTED_GAP = (53144999, 13513175949)
EXPECTED_CASE_DIGEST = "4fed6e77c92b1b163be8a0a2060039c002ca29025a42dc230a4f27801936a2c2"
EXPECTED_LITERAL_SIGNS = ((-1, 1319), (1, 1314))
EXPECTED_LITERAL_CHECKS = 26880
EXPECTED_LOW_LABELS = (243, 351)
EXPECTED_CARRIER_CELLS = 67700
EXPECTED_CARRIER_SHA256 = "02c3c649a0cdd6808479821077de9fe142fde6efc7fc6e3daa5380276d403b93"
EXPECTED_ENDPOINT = (0, 35100, 11, "L")
EXPECTED_DIRECT_CASE_SHA256 = "d75fcc91dfa1ed6c2e1240c609376f74e077457d61caac16da36a36d9ff89b91"
EXPECTED_DIRECT_RECORD_SHA256 = "a8996fee27e92e82100d42fab22908b78f1443c1f58c48db57d5636ded77c117"
EXPECTED_SUPPORTS = (
    (10296, 10296, 1471, 8825, "74673614815c2bcf6fbdb23d98e5e20eb4edbf37d3c40243eb3d3b199a39bf37"),
    (51480, 39624, 7355, 32269, "e28a6df29f568e4b3835555cfe11ab504413cea3202a19e7f5f575b86eb2fb5f"),
)

LEDGER_BEFORE = 374322
LAYER_ROWS = 9
LEDGER_AFTER = 374313
NEXT_CAP = 230
EXPECTED_SEMANTIC_SHA256 = "86fdf94cb22cd0353a28f8a47251a4c175193188de0fe78676c6e9e7d6ff899b"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name="thm3106_z231_base"):
    spec = importlib.util.spec_from_file_location(name, SOURCE_3106)
    require(spec is not None and spec.loader is not None, SOURCE_3106)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def isolated_screen_worker(task):
    return load(f"thm3106_z231_worker_{task[1]}").screen_worker(task)


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def atlas_rows(module):
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    selected = {LEVEL: [], NEXT_LEVEL: []}
    row_lines = tuple(
        line for line in module.thm.ATLAS.read_text().splitlines() if line.startswith("row=")
    )
    require(len(row_lines) == 6060, len(row_lines))
    for line in row_lines:
        match = pattern.match(line)
        require(match is not None, line)
        first = int(match.group(4))
        if first not in selected:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        L = int(match.group(2))
        high = int(match.group(3))
        selected[first].append((body, L, high, first < high))
    current = tuple(selected[LEVEL])
    following = tuple(selected[NEXT_LEVEL])
    current_census = (len(current), sum(x[3] for x in current), sum(not x[3] for x in current))
    next_census = (len(following), sum(x[3] for x in following), sum(not x[3] for x in following))
    require(current_census == EXPECTED_LAYER_CENSUS, current_census)
    require(next_census == EXPECTED_NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(current).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, current)
    require(hashlib.sha256(repr(following).encode()).hexdigest() == EXPECTED_NEXT_ROW_ORDER_SHA256, following)
    return current, following, current_census, next_census


def run_screen(rows, processes):
    tasks = tuple((LEVEL, body, L, high, wall) for body, L, high, wall in rows)
    if processes == 1:
        pairs = tuple(isolated_screen_worker(task) for task in tasks)
    else:
        context = mp.get_context("spawn")
        with context.Pool(processes=processes) as pool:
            pairs = tuple(pool.map(isolated_screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost or duplicated screen task")
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(totals == EXPECTED_SCREEN, totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "unverified status row")
    record_sha = hashlib.sha256(repr(screened).encode()).hexdigest()
    require(record_sha == EXPECTED_SCREEN_SHA256, record_sha)
    residual = tuple((row[1], row[13]) for row in screened if row[12])
    require(residual == ((RESIDUAL_BODY, RESIDUAL_BANK),), residual)
    require(hashlib.sha256(repr(residual).encode()).hexdigest() == EXPECTED_RESIDUAL_RECORD_SHA256, residual)
    return screened, totals, farkas, record_sha


def direct_cells(modulus, labels):
    cells = []
    for cell in range(modulus):
        if all(
            14 * (cell * label % modulus) >= modulus
            and 14 * ((cell * label % modulus) + label) <= 13 * modulus
            for label in labels
        ):
            cells.append(cell)
    return tuple(cells)


def terminal_audit(module):
    terminal = module.thm.terminal_probe((LEVEL, RESIDUAL_BODY, RESIDUAL_BANK))
    terminal_sha = hashlib.sha256(repr(terminal).encode()).hexdigest()
    require(terminal_sha == EXPECTED_TERMINAL_SHA256, terminal_sha)
    require((terminal[5].numerator, terminal[5].denominator) == EXPECTED_GAP, terminal[5])
    require(terminal[17] == EXPECTED_CASE_DIGEST and terminal[19] and terminal[20], terminal)
    require(
        (terminal[7], terminal[8], terminal[9], terminal[10], terminal[11], terminal[12], terminal[14], terminal[15], terminal[16])
        == (2, 2, 1, 2, 0, 0, 0, 0, 464),
        terminal,
    )

    eng = module.eng
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(RESIDUAL_BODY)
    needed = {d for ds in RESIDUAL_BANK for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, signs, checks = eng.build_literal_tables(stream, needed)
    require(signs == EXPECTED_LITERAL_SIGNS and checks == EXPECTED_LITERAL_CHECKS, (signs, checks))
    cases = eng.one_high_cases(stream, RESIDUAL_BANK, low, high)
    require(hashlib.sha256(repr(cases).encode()).hexdigest() == EXPECTED_DIRECT_CASE_SHA256, cases)
    label_sets = {tuple(sorted(label for _d, label in rows)) for _ds, _h, rows, _e in cases}
    require(label_sets == {EXPECTED_LOW_LABELS}, label_sets)

    labels = (*RESIDUAL_BODY, LEVEL, *EXPECTED_LOW_LABELS)
    cells = direct_cells(stream.L, labels)
    require(len(cells) == EXPECTED_CARRIER_CELLS, len(cells))
    carrier_sha = hashlib.sha256(repr(cells).encode()).hexdigest()
    require(carrier_sha == EXPECTED_CARRIER_SHA256, carrier_sha)
    require(cells == eng.fixed_safe_cells(stream, EXPECTED_LOW_LABELS), "direct/scalar carrier mismatch")
    vector = tuple(
        map(int, module.thm.vector_fixed_safe_cells(stream, EXPECTED_LOW_LABELS).tolist())
    )
    require(cells == vector, "direct/vector carrier mismatch")

    endpoint = None
    for cell in cells:
        for label in labels:
            residue = cell * label % stream.L
            for candidate in (
                (14 * residue - stream.L, cell, label, "L"),
                (13 * stream.L - 14 * (residue + label), cell, label, "R"),
            ):
                if endpoint is None or candidate < endpoint:
                    endpoint = candidate
    require(endpoint == EXPECTED_ENDPOINT, endpoint)

    records = []
    supports = []
    for ds, high_d, low_rows, excess in cases:
        support = tuple(sorted({cell % high_d for cell in cells}))
        capacity = stream.L // high_d
        coarse = (len(cells) + capacity - 1) // capacity
        kappa = (high_d + 6) // 7
        support_sha = hashlib.sha256(repr(support).encode()).hexdigest()
        records.append(
            (
                ds,
                high_d,
                low_rows,
                excess,
                capacity,
                coarse,
                len(support),
                kappa,
                len(support) - kappa,
                support_sha,
            )
        )
        supports.append((high_d, len(support), kappa, len(support) - kappa, support_sha))
    records = tuple(records)
    require(tuple(supports) == EXPECTED_SUPPORTS, supports)
    require(hashlib.sha256(repr(records).encode()).hexdigest() == EXPECTED_DIRECT_RECORD_SHA256, records)
    return terminal, terminal_sha, cases, cells, carrier_sha, endpoint, records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3106) == SOURCE_3106_SHA256, "THM-3106 source changed")
    require(sha(OUTPUT_3106) == OUTPUT_3106_SHA256, "THM-3106 output changed")
    require(f"semantic_sha256={SEMANTIC_3106_SHA256}" in OUTPUT_3106.read_text(), "THM-3106 semantic changed")
    require(hashlib.sha256(repr(RESIDUAL_BANK).encode()).hexdigest() == EXPECTED_RESIDUAL_SHA256, "residual transcription")

    module = load()
    rows, next_rows, census, next_census = atlas_rows(module)
    screened, totals, farkas, screen_sha = run_screen(rows, args.processes)
    terminal, terminal_sha, cases, cells, carrier_sha, endpoint, records = terminal_audit(module)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z231-screen-terminal-v1",
        SOURCE_3106_SHA256,
        OUTPUT_3106_SHA256,
        SEMANTIC_3106_SHA256,
        module.thm.ATLAS_SHA256,
        rows,
        next_rows,
        screened,
        terminal,
        cases,
        carrier_sha,
        endpoint,
        records,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z231 exact screen and complete-cell cardinality descent",
        f"dependency=THM3106_source:{SOURCE_3106_SHA256};output:{OUTPUT_3106_SHA256};semantic:{SEMANTIC_3106_SHA256}",
        f"atlas=sha256:{module.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256}",
        f"screen=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha}",
    ]
    for row in screened:
        lines.append(
            f"ROW;E={row[1]};L={row[2]};high={row[3]};branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.extend(
        [
            f"RESIDUAL;E={RESIDUAL_BODY};masks:{len(RESIDUAL_BANK)};sha256:{EXPECTED_RESIDUAL_SHA256};tuples:{RESIDUAL_BANK}",
            f"terminal=two_high_gap:{ftext(terminal[5])};gap_witness:{terminal[6]};zero_high_hostiles:{terminal[7]};one_high_cases:{terminal[8]};low_label_sets:{terminal[9]};coarse_cardinality:{terminal[10]};exact_cardinality:{terminal[11]};failures:{terminal[14]};minimum_certificate_slack:{terminal[16]};case_sha256:{terminal[17]};terminal_sha256:{terminal_sha}",
            f"hostile_terminal_audit=direct_full_grid_equals_scalar_equals_vector;labels:{EXPECTED_LOW_LABELS};cells:{len(cells)};carrier_sha256:{carrier_sha};minimum_weak_endpoint:{endpoint};record_sha256:{EXPECTED_DIRECT_RECORD_SHA256}",
        ]
    )
    for record in records:
        lines.append(
            f"CASE;ds={record[0]};high_d={record[1]};low_rows={record[2]};excess={ftext(record[3])};capacity={record[4]};coarse_lower_bound={record[5]};support={record[6]};kappa={record[7]};actual_slack={record[8]};support_sha256={record[9]}"
        )
    lines.extend(
        [
            "direction_screen=ray_quotient_and_Farkas_checks_close_a_superset_of_every_actual_projected_assignment",
            "direction_terminal=duplicate_permitting_two_high_gap_plus_wall_at_least_one_high_forces_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_therefore_enlarges_the_actual_assignment_set",
            "direction_carrier=weak_complete_cells_include_strict_open_safe_points;projected_support_above_ceil(d/7)_forces_the_THM2941_completed_carrier_and_THM1166_contradiction",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP};next_layer:z1={NEXT_LEVEL}_rows:{next_census[0]}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{EXPECTED_NEXT_ROW_ORDER_SHA256};status:occupied_unscouted_handoff",
            "scope=projected_k3_necessary_atlas_only;z1=230_is_occupied_but_unscouted;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", newline="\n")
    print("\n".join(lines))


if __name__ == "__main__":
    mp.freeze_support()
    main()
