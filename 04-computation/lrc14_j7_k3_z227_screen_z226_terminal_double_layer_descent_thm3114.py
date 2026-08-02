#!/usr/bin/env python3
"""Exact projected-k3 z1=227 screen and z1=226 terminal descent."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3113 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.py"
)
OUTPUT_3113 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.out"
)

SOURCE_3113_SHA256 = "1e23ec19fa147c55fb6d38a965eedae0132f5e069b9f820bfd5c300dce4d8f89"
OUTPUT_3113_SHA256 = "5b6329e5a47a41031c2d9d3f9def36deb25509d801dbb3543a1271de0bd2e3fb"
SEMANTIC_3113_SHA256 = "093ed7cea617666cae3463ae08c3429943c93dbd17b9af7a02d4601025692482"

LEVELS = (227, 226)
NEXT_LEVEL = 225
EXPECTED_CENSUS = {227: (30, 28, 2), 226: (13, 13, 0), 225: (78, 78, 0)}
EXPECTED_ROW_ORDER_SHA256 = {
    227: "17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c",
    226: "ec2f1218d61dd69d3811d68eb87a23254d276bb75647be2d4c883affe53a520e",
    225: "9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719",
}
EXPECTED_SCREEN = {227: (284, 143, 141, 0), 226: (979, 526, 363, 90)}
EXPECTED_ORDER = {227: (9, 9, 0, 0), 226: (0, 0, 0, 0)}
EXPECTED_FARKAS = {227: (0, 141), 226: (0, 363)}
EXPECTED_SCREEN_SHA256 = {
    227: "b4da9b80d12c969cf5561a0e30e36d282a955e77ff32a5f88fb9514cf78fb414",
    226: "5a64104b1ca45e7dd4fe08cfc5e0bcc5a3226a7f05740dc20f361d664a2d98a7",
}
EXPECTED_TERMINAL = (3, 3, 3, 84, 90, 84, 6, 0, 0, 0)
EXPECTED_TERMINAL_SHA256 = "f343c72ad5b49b37eca86171374b274be6ef48bd9e28115107daf9ccb140c5a3"
EXPECTED_CARRIER_SHA256 = "f5060567a3cc0341902c541ec6910161d509b48d16f5303a467095c40fbab66f"
EXPECTED_DIRECT_SHA256 = "6eaf93061816e63631a0e5180697cc1be7fe69980705d2d87f0a4a822698cfa9"
EXPECTED_SEMANTIC_SHA256 = "ad9c2724b6468d586ae70ab120f57519350c7de440c66a3765609bd1a7880d51"

LEDGER_BEFORE = 374215
LAYER_ROWS = 43
LEDGER_AFTER = 374172
NEXT_CAP = 225

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


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3113):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3114_worker_{task[0]}_{task[1]}").screen_worker(task)


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
            hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256[level],
            (level, "row order"),
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
    require(len(by_task) == len(tasks), (level, "lost task"))
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN[level], (level, totals))
    order = tuple(row for row in screened if not row[5])
    order_totals = tuple(sum(row[index] for row in order) for index in (9, 10, 11, 12))
    require(order_totals == EXPECTED_ORDER[level], (level, order_totals))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS[level], (level, farkas))
    require(all(row[16] == row[11] for row in screened), (level, "status verification"))
    # Exact certificate verification stays load-bearing, but MISTAKE-331
    # forbids hashing row[21], which binds a solver-selected dual basis and
    # contradiction scale.  THM-3078's first nineteen fields are the
    # canonical problem/result record.
    canonical_screened = tuple(row[:19] for row in screened)
    record_sha = hashlib.sha256(repr(canonical_screened).encode()).hexdigest()
    if EXPECTED_SCREEN_SHA256[level] is not None:
        require(record_sha == EXPECTED_SCREEN_SHA256[level], (level, record_sha))
    residual_rows = tuple(row for row in screened if row[12])
    if level == 227:
        require(not residual_rows, residual_rows)
    else:
        require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), residual_rows)
        require(all(row[5] for row in residual_rows), "non-wall residual")
        for row in residual_rows:
            expected = RESIDUAL[row[1]]
            require(row[12] == expected["masks"], (row[1], row[12]))
            require(
                hashlib.sha256(repr(row[13]).encode()).hexdigest() == expected["mask_sha"],
                (row[1], "mask digest"),
            )
    return screened, totals, order_totals, farkas, record_sha, residual_rows


def run_terminals(base, residual_rows):
    terminals = tuple(base.thm.terminal_probe((226, row[1], row[13])) for row in residual_rows)
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
                expected["zero"],
                expected["cases"],
                1,
                expected["coarse"],
                expected["exact"],
                0,
                0,
                1,
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
        eng.FIRST = 226
        eng.ray.FIRST = 226
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
            direct = predecessor.direct_cells(stream.L, (*body, 226, *labels))
            scalar = eng.fixed_safe_cells(stream, labels)
            vector = tuple(map(int, base.thm.vector_fixed_safe_cells(stream, labels).tolist()))
            require(direct == scalar == vector, (body, labels, "carrier mismatch"))
            endpoint = None
            for cell in direct:
                for label in (*body, 226, *labels):
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
    require(len(carrier_records) == 3 and len(case_records) == 90, "direct census")
    require(sum(row[-2] == "coarse-cardinality" for row in case_records) == 84, "coarse count")
    require(sum(row[-2] == "exact-cardinality" for row in case_records) == 6, "exact count")
    require(
        min(row[8] - row[10] for row in case_records if row[-2] == "coarse-cardinality") == 1,
        "coarse margin",
    )
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
    require(sha(SOURCE_3113) == SOURCE_3113_SHA256, "THM-3113 source changed")
    require(sha(OUTPUT_3113) == OUTPUT_3113_SHA256, "THM-3113 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3113_SHA256}" in OUTPUT_3113.read_text(),
        "THM-3113 semantic changed",
    )

    source = load("thm3114_source")
    driver = source.load("thm3114_driver")
    predecessor = driver.load("thm3114_predecessor")
    base = predecessor.load("thm3114_base")
    layers = atlas_layers(base)
    screens = {level: run_screen(level, layers[level], args.processes) for level in LEVELS}
    terminals, terminal_totals, terminal_sha = run_terminals(base, screens[226][5])
    carriers, cases, carrier_sha, direct_sha = direct_carrier_audit(
        predecessor, base, screens[226][5]
    )
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z227-screen-z226-terminal-double-descent-v2",
        SOURCE_3113_SHA256,
        OUTPUT_3113_SHA256,
        SEMANTIC_3113_SHA256,
        base.thm.ATLAS_SHA256,
        layers,
        tuple(tuple(row[:19] for row in screens[level][0]) for level in LEVELS),
        tuple((screens[level][1], screens[level][2], screens[level][3]) for level in LEVELS),
        terminals,
        carriers,
        cases,
        LEDGER_BEFORE,
        LEDGER_AFTER,
        NEXT_CAP,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z227 screen and z226 terminal double-layer descent",
        f"dependency=THM3113_source:{SOURCE_3113_SHA256};output:{OUTPUT_3113_SHA256};semantic:{SEMANTIC_3113_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
    ]
    for level in LEVELS:
        _screened, totals, order_totals, farkas, screen_sha, residual_rows = screens[level]
        census = EXPECTED_CENSUS[level]
        lines.append(
            f"layer=z1:{level};rows:{census[0]};wall:{census[1]};order:{census[2]};"
            f"row_order_sha256:{EXPECTED_ROW_ORDER_SHA256[level]};states:{totals[0]};"
            f"crude:{totals[1]};status:{totals[2]};residual:{totals[3]};"
            f"direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};"
            f"screen_record_sha256:{screen_sha};order_control:{order_totals};"
            f"residual_rows:{len(residual_rows)}"
        )
    lines.append(
        "terminal=z1:226;rows:3;positive_two_high_gap:3;closed:3;"
        f"zero_high_hostiles:{terminal_totals[3]};one_high_cases:{terminal_totals[4]};"
        f"coarse:{terminal_totals[5]};exact:{terminal_totals[6]};maxgap:0;failures:0;"
        f"terminal_record_sha256:{terminal_sha}"
    )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};"
            f"two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};"
            f"one_high_cases={row[8]};coarse={row[10]};exact={row[11]};"
            f"failed={row[14]};minimum_certificate_slack={row[16]};"
            f"case_certificate_sha256={row[17]}"
        )
    lines.append(
        "direct_carrier_audit=carriers:3;cases:90;coarse:84;exact:6;"
        f"minimum_coarse_margin:1;minimum_support_slack:1;all_exact_supports:full;"
        f"carrier_record_sha256:{carrier_sha};case_record_sha256:{direct_sha}"
    )
    for body, labels, count, cell_sha, endpoint in carriers:
        lines.append(
            f"CARRIER;E={body};lows={labels};cells={count};cells_sha256={cell_sha};"
            f"minimum_weak_endpoint={endpoint}"
        )
    exact_denominators = tuple(row[2] for row in cases if row[-2] == "exact-cardinality")
    exact_slacks = tuple(row[-3] for row in cases if row[-2] == "exact-cardinality")
    lines.extend(
        [
            f"exact_supports=denominators:{exact_denominators};support_minus_kappa:{exact_slacks}",
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment",
            "direction_terminal=at_z226_each_strictly_positive_duplicate_permitting_two_high_gap_together_with_the_wall_at_least_one_high_gate_forces_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=each_direct_complete_cell_is_an_inner_carrier_wholly_contained_in_the_strict_open_safe_set;the_coarse_bound_or_exact_full_projected_support_above_ceil(d/7)_forces_the_completed_carrier_contradiction",
            "evidence_boundary=all_returned_Farkas_certificates_are_verified_exactly;screen_and_semantic_digests_bind_only_the_canonical_19_field_problem_result_rows_plus_basis_invariant_branch_counts_and_never_the_solver_selected_dual_or_contradiction_magnitude",
            "transition=z227_dies_at_the_screen_but_z226_reenters_and_is_closed_by_the_terminal",
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
