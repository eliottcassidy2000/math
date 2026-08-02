#!/usr/bin/env python3
"""Exact projected-k3 z1=230 screen and compressed complete-cell descent.

All fifty atlas rows are recomputed through the promoted THM-3109/THM-3078
screen.  The six residual bodies are then closed by the inherited
duplicate-two-high/one-high terminal.  A separate scalar audit groups the
seventy-one one-high cases by their eight invariant low-label carriers and
rebuilds every complete cell and projected high-label support directly.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3109 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py"
)
OUTPUT_3109 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.out"
)

SOURCE_3109_SHA256 = "f40df5b485883111fa2245cbb76fc43d8f90418a1529ff022adfeb73a1972e83"
OUTPUT_3109_SHA256 = "9126830064e15e7e7f4c7128f8cbb968c95e660a8dbf9171657d23c3aeb87e66"
SEMANTIC_3109_SHA256 = "86fdf94cb22cd0353a28f8a47251a4c175193188de0fe78676c6e9e7d6ff899b"

LEVEL = 230
NEXT_LEVEL = 229
EXPECTED_CENSUS = (50, 48, 2)
EXPECTED_NEXT_CENSUS = (43, 36, 7)
EXPECTED_ROW_ORDER_SHA256 = "f5129616fe3d7f5884350abe7e77bb2684798e213600e266cd456bc9c582f058"
EXPECTED_NEXT_ROW_ORDER_SHA256 = "7449dd7ad70cf3c76c32edb2cc509e29989ac008c2e9a9968bceaabf3e2a2587"
EXPECTED_SCREEN = (4156, 2437, 1648, 71)
EXPECTED_ORDER = (8, 8, 0, 0)
EXPECTED_FARKAS = (0, 1648)
EXPECTED_SCREEN_SHA256 = "b74351d30dced85c2697bf8da85ca0ec99f500014b014ed13bce80921f96cfba"
EXPECTED_TERMINAL = (6, 6, 6, 68, 71, 68, 3, 0, 0, 0)
EXPECTED_TERMINAL_SHA256 = "90b1768cf0165683ceb086359180563081a0ff89a3876c88a47c7240ca036e37"
EXPECTED_CARRIER_SHA256 = "a8e749413d276397dd8ba521a62bad219e82639d7c53669030c41ec166ebeb5d"
EXPECTED_DIRECT_SHA256 = "c2640270e1af368ea4c40856ef9b8254d16214bea84ea2be6197f75f6d17e761"
EXPECTED_SEMANTIC_SHA256 = "b4dd4e0e66d0429cac12dc96d1d769d3f05ab3c503e1dcbb1b94f5e6020c9b6c"

LEDGER_BEFORE = 374313
LAYER_ROWS = 50
LEDGER_AFTER = 374263
NEXT_CAP = 229

RESIDUAL = {
    (1, 7, 9, 10, 11, 12): {
        "masks": 19,
        "mask_sha": "a4019536be949152e471d5c4944a0541795113a3910ff8f849931ea74ef4c6aa",
        "gap": "113487288713/32897729721780",
        "zero": 18,
        "cases": 19,
        "groups": 1,
        "coarse": 18,
        "exact": 1,
        "raw_case_sha": "4fa16025fbe6008dd2cdad07c2dd98e5b750abad9c381d58161e33d793e44c72",
    },
    (1, 7, 10, 11, 12, 13): {
        "masks": 12,
        "mask_sha": "6d05ac0e3fb722cafc60b7b1c878e91d70f8be2580cb2046da6466bf00f13164",
        "gap": "5173786082623/1879348592763460",
        "zero": 12,
        "cases": 12,
        "groups": 1,
        "coarse": 12,
        "exact": 0,
        "raw_case_sha": "077f42f51e555de11bca613749de70104ee8817e3eca171e44cd2e0be492d5a7",
    },
    (1, 8, 10, 11, 12, 14): {
        "masks": 4,
        "mask_sha": "ae54b5f019c3cf52cb420cea79de8aed2a3df7a71258e86e4f59c4d40a5b0786",
        "gap": "1349719633/341778176376",
        "zero": 4,
        "cases": 4,
        "groups": 1,
        "coarse": 4,
        "exact": 0,
        "raw_case_sha": "d129b445d8ab706ce0ae8332c5c669f73d965bcd232cbf5e3c1033d0dc4c0f9d",
    },
    (1, 9, 10, 11, 12, 14): {
        "masks": 34,
        "mask_sha": "4fcbf49394b967b6e3f16215480700de5aeb81ec8e258b5656db87462ce5db20",
        "gap": "8059115839/2508612085980",
        "zero": 34,
        "cases": 34,
        "groups": 3,
        "coarse": 34,
        "exact": 0,
        "raw_case_sha": "ea706e6fa52eb4dac3ee5651cc40dc4f677fe8c623701b3ee14fc4c5a278902f",
    },
    (2, 9, 10, 11, 12, 14): {
        "masks": 1,
        "mask_sha": "bc26567858a23a670523ab601601f5c9326e3113a76cb443b650a60891e03e98",
        "gap": "14389939/3132625860",
        "zero": 0,
        "cases": 1,
        "groups": 1,
        "coarse": 0,
        "exact": 1,
        "raw_case_sha": "3d193aba288bf77396a83a070191395675c58a68e758e03e3706612c96d1a205",
    },
    (3, 9, 10, 11, 12, 14): {
        "masks": 1,
        "mask_sha": "7e0bd979705f1f20992aef38f71b77110235cde264d01c9a97e57829d47d6f7a",
        "gap": "49827592/15730119405",
        "zero": 0,
        "cases": 1,
        "groups": 1,
        "coarse": 0,
        "exact": 1,
        "raw_case_sha": "5328d66ede80bcc47ae63277619eaaf774e45df07a852d3fa9b356fa3b1c72a9",
    },
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3109):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3109_z230_worker_{task[1]}").isolated_screen_worker(task)


def atlas_layers(predecessor, base):
    _old, current, _old_census, census = predecessor.atlas_rows(base)
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(current).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, current)

    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    following = []
    for line in base.thm.ATLAS.read_text().splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        first = int(match.group(4))
        if first != NEXT_LEVEL:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        high = int(match.group(3))
        following.append((body, int(match.group(2)), high, first < high))
    following = tuple(following)
    next_census = (
        len(following),
        sum(row[3] for row in following),
        sum(not row[3] for row in following),
    )
    require(next_census == EXPECTED_NEXT_CENSUS, next_census)
    require(
        hashlib.sha256(repr(following).encode()).hexdigest()
        == EXPECTED_NEXT_ROW_ORDER_SHA256,
        following,
    )
    return current, following, census, next_census


def run_screen(rows, processes):
    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        context = mp.get_context("spawn")
        with context.Pool(processes=processes) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost or duplicated screen task")
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    order = tuple(row for row in screened if not row[5])
    order_totals = tuple(sum(row[index] for row in order) for index in (9, 10, 11, 12))
    require(order_totals == EXPECTED_ORDER, order_totals)
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "unverified status row")
    record_sha = hashlib.sha256(repr(screened).encode()).hexdigest()
    if EXPECTED_SCREEN_SHA256 is not None:
        require(record_sha == EXPECTED_SCREEN_SHA256, record_sha)

    residual_rows = tuple(row for row in screened if row[12])
    require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), residual_rows)
    require(all(row[5] for row in residual_rows), "non-wall residual body")
    for row in residual_rows:
        expected = RESIDUAL[row[1]]
        require(row[12] == expected["masks"], (row[1], row[12]))
        require(
            hashlib.sha256(repr(row[13]).encode()).hexdigest() == expected["mask_sha"],
            (row[1], "residual digest"),
        )
    return screened, totals, order_totals, farkas, record_sha, residual_rows


def run_terminals(base, residual_rows):
    terminals = tuple(
        base.thm.terminal_probe((LEVEL, row[1], row[13])) for row in residual_rows
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
    if EXPECTED_TERMINAL_SHA256 is not None:
        require(record_sha == EXPECTED_TERMINAL_SHA256, record_sha)
    for row in terminals:
        expected = RESIDUAL[row[1]]
        require(row[4] == expected["masks"], (row[1], "terminal masks"))
        require(ftext(row[5]) == expected["gap"] and row[5] > 0, (row[1], row[5]))
        require(
            (row[7], row[8], row[9], row[10], row[11], row[12], row[14], row[20])
            == (
                expected["zero"],
                expected["cases"],
                expected["groups"],
                expected["coarse"],
                expected["exact"],
                0,
                0,
                True,
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
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(body)
        needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
        low, high, _signs, _checks = eng.build_literal_tables(stream, needed)
        cases = eng.one_high_cases(stream, residual, low, high)
        expected = RESIDUAL[body]
        require(len(cases) == expected["cases"], (body, "raw case count"))
        require(
            hashlib.sha256(repr(cases).encode()).hexdigest()
            == expected["raw_case_sha"],
            (body, "raw case digest"),
        )
        require(
            len({
                tuple(sorted(label for _d, label in low_rows))
                for _ds, _high_d, low_rows, _excess in cases
            })
            == expected["groups"],
            (body, "carrier group count"),
        )
        cache = {}
        for _ds, _high_d, low_rows, _excess in cases:
            labels = tuple(sorted(label for _d, label in low_rows))
            if labels in cache:
                continue
            direct = predecessor.direct_cells(stream.L, (*body, LEVEL, *labels))
            scalar = eng.fixed_safe_cells(stream, labels)
            vector = tuple(
                map(int, base.thm.vector_fixed_safe_cells(stream, labels).tolist())
            )
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
            require(endpoint is not None and endpoint[0] >= 0, (body, labels, endpoint))
            cache[labels] = direct
            carrier_records.append(
                (
                    body,
                    labels,
                    len(direct),
                    hashlib.sha256(repr(direct).encode()).hexdigest(),
                    endpoint,
                )
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
            case_records.append(
                (
                    body,
                    ds,
                    high_d,
                    low_rows,
                    excess,
                    labels,
                    len(cells),
                    capacity,
                    coarse,
                    len(support),
                    kappa,
                    len(support) - kappa,
                    method,
                    hashlib.sha256(repr(support).encode()).hexdigest(),
                )
            )

    carrier_records = tuple(carrier_records)
    case_records = tuple(case_records)
    require(len(carrier_records) == 8, len(carrier_records))
    require(len(case_records) == 71, len(case_records))
    require(sum(row[-2] == "coarse-cardinality" for row in case_records) == 68, "coarse count")
    require(sum(row[-2] == "exact-cardinality" for row in case_records) == 3, "exact count")
    require(
        min(
            row[8] - row[10]
            for row in case_records
            if row[-2] == "coarse-cardinality"
        )
        == 2,
        "minimum coarse lower-bound margin",
    )
    require(min(row[-3] for row in case_records) == 1, "minimum support slack")
    carrier_sha = hashlib.sha256(repr(carrier_records).encode()).hexdigest()
    direct_sha = hashlib.sha256(repr(case_records).encode()).hexdigest()
    if EXPECTED_CARRIER_SHA256 is not None:
        require(carrier_sha == EXPECTED_CARRIER_SHA256, carrier_sha)
    if EXPECTED_DIRECT_SHA256 is not None:
        require(direct_sha == EXPECTED_DIRECT_SHA256, direct_sha)
    return carrier_records, case_records, carrier_sha, direct_sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3109) == SOURCE_3109_SHA256, "THM-3109 source changed")
    require(sha(OUTPUT_3109) == OUTPUT_3109_SHA256, "THM-3109 output changed")
    require(f"semantic_sha256={SEMANTIC_3109_SHA256}" in OUTPUT_3109.read_text(), "THM-3109 semantic changed")

    predecessor = load("thm3109_z230_driver")
    base = predecessor.load("thm3106_z230_base")
    rows, next_rows, census, next_census = atlas_layers(predecessor, base)
    screened, screen_totals, order_totals, farkas, screen_sha, residual_rows = run_screen(rows, args.processes)
    terminals, terminal_totals, terminal_sha = run_terminals(base, residual_rows)
    carriers, cases, carrier_sha, direct_sha = direct_carrier_audit(predecessor, base, residual_rows)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z230-screen-compressed-carrier-v1",
        SOURCE_3109_SHA256,
        OUTPUT_3109_SHA256,
        SEMANTIC_3109_SHA256,
        base.thm.ATLAS_SHA256,
        rows,
        next_rows,
        screened,
        terminals,
        carriers,
        cases,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z230 exact screen and compressed complete-cell descent",
        f"dependency=THM3109_source:{SOURCE_3109_SHA256};output:{OUTPUT_3109_SHA256};semantic:{SEMANTIC_3109_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256}",
        f"screen=states:{screen_totals[0]};crude:{screen_totals[1]};status:{screen_totals[2]};residual:{screen_totals[3]};residual_bodies:{len(residual_rows)};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha}",
        f"order_control=rows:{census[2]};states:{order_totals[0]};crude:{order_totals[1]};status:{order_totals[2]};residual:{order_totals[3]}",
    ]
    for row in screened:
        lines.append(
            f"ROW;E={row[1]};L={row[2]};high={row[3]};branch={'wall' if row[5] else 'order'};states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};direct={row[19]};legacy={row[20]};stage_sha256={row[18]};residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.append(
        f"terminal=rows:{terminal_totals[0]};positive_two_high_gap:{terminal_totals[1]};closed:{terminal_totals[2]};zero_high_hostiles:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};coarse_cardinality:{terminal_totals[5]};exact_cardinality:{terminal_totals[6]};maxgap:{terminal_totals[7]};failures:{terminal_totals[8]};unit_checks:{terminal_totals[9]};terminal_record_sha256:{terminal_sha}"
    )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse_cardinality={row[10]};exact_cardinality={row[11]};failed={row[14]};minimum_certificate_slack={row[16]};raw_case_sha256={RESIDUAL[row[1]]['raw_case_sha']};case_certificate_sha256={row[17]}"
        )
    lines.append(
        f"compressed_carrier_audit=body_low_label_carriers:{len(carriers)};cases:{len(cases)};coarse:68;exact:3;minimum_coarse_margin:2;minimum_support_slack:1;carrier_record_sha256:{carrier_sha};case_record_sha256:{direct_sha}"
    )
    for body, labels, count, cells_sha, endpoint in carriers:
        lines.append(
            f"CARRIER;E={body};lows={labels};cells={count};cells_sha256={cells_sha};minimum_weak_endpoint={endpoint}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment",
            "direction_terminal=the_positive_duplicate_permitting_two_high_gap_and_the_wall_at_least_one_high_gate_force_exactly_one_high;zero_high_scalar_passes_are_hostiles_excluded_by_the_wall;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=each_direct_complete_cell_is_an_inner_carrier_wholly_contained_in_the_strict_open_safe_set;grouping_by_low_label_pair_reuses_the_same_carrier_without_identifying_distinct_high_cases;the_coarse_lower_bound_or_the_exact_projected_support_above_ceil(d/7)_forces_the_completed_carrier_contradiction",
            "boundary=the_compression_forgets_joint_two_high_compatibility_and_requires_a_separate_common_two_high_carrier_if_a_future_duplicate_two_high_gap_is_nonpositive",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP};next_layer:z1={NEXT_LEVEL}_rows:{next_census[0]}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{EXPECTED_NEXT_ROW_ORDER_SHA256};status:occupied_unscouted_handoff",
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
