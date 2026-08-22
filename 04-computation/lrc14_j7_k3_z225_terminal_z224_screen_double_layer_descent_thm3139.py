#!/usr/bin/env python3
"""Exact projected-k3 z225 terminal plus z224 screen double-layer descent."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3114 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.py"
)
OUTPUT_3114 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
)

SOURCE_3114_SHA256 = "8de1e3d03b5070a84b040ac13a173a598646107f85e7ba0defc2ca070808f162"
OUTPUT_3114_SHA256 = "c6cfd6e9b6bd9a89c2d679d86f9fe545d949269097698f51daa86b30350a543f"
SEMANTIC_3114_SHA256 = "ad9c2724b6468d586ae70ab120f57519350c7de440c66a3765609bd1a7880d51"

LEVEL = 225
SECOND_LEVEL = 224
EXPECTED_CENSUS = (78, 78, 0)
EXPECTED_ROW_ORDER_SHA256 = "9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719"
EXPECTED_SCREEN = (30615, 15904, 12813, 1898)
EXPECTED_FARKAS = (0, 12813)
EXPECTED_SCREEN_SHA256 = "75de10a09defef3d4d388b2eb21511d17f80ddffee274318677bc7cc7bc71cd7"
EXPECTED_TERMINAL = (24, 1898, 1754, 3025, 163, 2822, 203, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "483850de20c0a76be91c682ec252cb0970170527a75fcf34dbac179bef168669"
EXPECTED_CLOSURE_SHA256 = "54222f669aae83b2f217b667520ed9b9e27059ebbb8f0b31d177232cbade226b"
EXPECTED_CASE_VECTOR_SHA256 = "151c16219d8b079f910963e19ff3b813f29584fdae2ecc8ec19969c469528cb6"
EXPECTED_SECOND_CENSUS = (4, 3, 1)
EXPECTED_SECOND_ROW_ORDER_SHA256 = "b42d3f6d42e76a10f89ba7d1fda6f76884081cba9cf32ca4fcedba4757811cb4"
EXPECTED_SECOND_SCREEN = (12, 2, 10, 0)
EXPECTED_SECOND_FARKAS = (0, 10)
EXPECTED_SECOND_SCREEN_SHA256 = "f6928f10e464ee694ebce756f7a6a2cb43360a7cd3bd121be4b6f5b98f91b21d"
EXPECTED_SEMANTIC_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"

LEDGER_BEFORE = 374172
LAYER_ROWS = 78 + 4
LEDGER_AFTER = 374090
NEXT_CAP = 223
NEXT_LEVEL = 223
NEXT_CENSUS = (65, 54, 11)
NEXT_ROW_ORDER_SHA256 = "58454a33efe0b74fddeddb498d1518033ccf74064b19ec611d63a223f6065566"

RESIDUAL = {
    (1, 2, 9, 10, 12, 14): (5, "0c52d36e000a66b3d3b4528ece4285e5a585b85e0b5b8a8cebe85395bfba0f80"),
    (1, 2, 9, 11, 12, 14): (6, "0c64bd346ca527b4d5bd94fdb11a6ee4c5456db8ce4a00265aafd341a2de9cfd"),
    (1, 4, 5, 9, 11, 13): (12, "0c2c73dd3c7f53a2354b5bf4c14b1b2bf469bf7884349d559f297d350b4c8496"),
    (1, 4, 5, 9, 11, 14): (54, "0ce76273248a7e4a6f7d18b3dc6a3822ad0e087dd1ae21a942870bc54afc99f7"),
    (1, 4, 9, 10, 11, 13): (11, "29572f4c4d1a69798ae75e294be29c9f527277a61314e78abd0339de943c7208"),
    (1, 4, 9, 10, 11, 14): (33, "72b55c092f3715d1f7bb921f7023659d17a1897a981f25731eefc52fc6aa531a"),
    (1, 4, 9, 10, 12, 14): (2, "ed5b319db1ec967d1661c6b58a2dd513eb5a8e10c34c59baf71be1b187410100"),
    (1, 4, 9, 11, 12, 13): (10, "b2863149069402179ad4b29e69a516ac183b3eac60809bcc19e53bdc0eaa4382"),
    (1, 4, 9, 11, 12, 14): (49, "9dacfc2c4516461ea0dd2aa7da763511e7d86a020b82c2c17b29fcbd8d0d1690"),
    (1, 4, 9, 11, 13, 14): (241, "beea643b0b32460fac60238314465c4497f922d22197d3e4620d20bbdb96469f"),
    (1, 5, 9, 11, 12, 14): (592, "531dabbf04975055646c5ea283799941abc05e7e5a4f1f12a428d4864fd82858"),
    (1, 5, 9, 11, 13, 14): (33, "c2c5b6c404710f18f04dd14caf5e5d79f792778358d4caae9d7771f99a50c07d"),
    (1, 9, 10, 11, 12, 14): (175, "eb3145a48626fe08baed48b6bc579682341f984c26f0e1cc7c3935e3b9df0a8f"),
    (1, 9, 11, 12, 13, 14): (306, "a8c28371ce2f0881d8e1df0090dbc70f96b1182803a76502de1185cf8bf0a8d0"),
    (2, 4, 5, 9, 11, 13): (2, "d8b37ccd3571f2fdc12b3e7694d16b9043d88be0063781c94664ea1f46aa5514"),
    (2, 4, 9, 11, 12, 13): (1, "1b0d31ed668cc4bdab490d69c3b6820ee741458f1a3506929ae2fedd3163a419"),
    (2, 4, 9, 11, 13, 14): (5, "8a2c5e88ac2358b2dfeb0bc6baa0ec072aaa5bc8ab1e944fa3afa76e71402e1c"),
    (2, 5, 9, 11, 12, 14): (87, "4c0a737eca24cccf631d02d454612583ba0105497e4d5bab79436e226446953b"),
    (2, 5, 9, 11, 13, 14): (27, "0f233e1bbbf67929d328322cbc9de2fd13cd6f6dcbd2d2a8eca1312353a4f51d"),
    (2, 8, 9, 11, 12, 14): (16, "a333b97261c634d74fcaea3ab46fca4fb7427cdd0a87e6357b9e6e898415d5a9"),
    (2, 9, 10, 11, 12, 14): (27, "f0712792992ce7ad341e5b2c05bdf55fa5b53d2a62336a3a945bc05c3a9098d5"),
    (2, 9, 11, 12, 13, 14): (182, "c5a98abd79fa5b9b541522c9637a55535c30a184d90219c1094f71163f2b7aa4"),
    (4, 5, 9, 11, 13, 14): (19, "27aafd15fb1f6c87502eb3a753d0b11d6175430d2a73b8ea81ad24b5df8ebab1"),
    (4, 9, 11, 12, 13, 14): (3, "6b74105488639c2740b870c9a4f31e7e8c485e8c48eb529740f3ab7d1ed503e6"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3114):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3139_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    index, body, bank = task
    source = load(f"thm3139_terminal_source_{index}")
    driver = source.load(f"thm3139_terminal_driver_{index}")
    predecessor = driver.load(f"thm3139_terminal_predecessor_{index}")
    bridge = predecessor.load(f"thm3139_terminal_bridge_{index}")
    base = bridge.load(f"thm3139_terminal_base_{index}")
    return index, base.thm.terminal_probe((LEVEL, body, bank))


def atlas_rows(base, level):
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    rows = []
    for line in base.thm.ATLAS.read_text().splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(4)) != level:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        high = int(match.group(3))
        rows.append((body, int(match.group(2)), high, level < high))
    return tuple(rows)


def terminal_semantic(row):
    # row[6] is a chosen duplicate-gap witness and is deliberately omitted.
    return (*row[:6], *row[7:])


def closure_row(row):
    return (
        row[1],
        row[4],
        (row[5].numerator, row[5].denominator),
        row[7],
        row[8],
        row[9],
        row[10],
        row[11],
        row[12],
        row[14],
        row[16],
        row[17],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=12)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3114) == SOURCE_3114_SHA256, "THM-3114 source changed")
    require(sha(OUTPUT_3114) == OUTPUT_3114_SHA256, "THM-3114 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3114_SHA256}" in OUTPUT_3114.read_text(),
        "THM-3114 semantic changed",
    )

    source = load("thm3139_source")
    driver = source.load("thm3139_driver")
    predecessor = driver.load("thm3139_predecessor")
    bridge = predecessor.load("thm3139_bridge")
    base = bridge.load("thm3139_base")

    rows = atlas_rows(base, LEVEL)
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_CENSUS, census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, "row order")
    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if args.processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost screen task")
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[5] and row[16] == row[11] for row in screened), "wall/status control")
    canonical_screen = tuple(row[:19] for row in screened)
    screen_sha = hashlib.sha256(repr(canonical_screen).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual_rows = tuple(row for row in screened if row[12])
    require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), "residual bodies")
    for row in residual_rows:
        count, mask_sha = RESIDUAL[row[1]]
        require(row[12] == count, (row[1], row[12]))
        require(hashlib.sha256(repr(row[13]).encode()).hexdigest() == mask_sha, (row[1], "mask sha"))

    terminal_tasks = tuple((index, row[1], row[13]) for index, row in enumerate(residual_rows))
    if args.processes == 1:
        terminal_pairs = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(terminal_tasks))) as pool:
            terminal_pairs = tuple(pool.map(terminal_worker, terminal_tasks))
    terminal_by_index = {index: row for index, row in terminal_pairs}
    terminals = tuple(terminal_by_index[index] for index in range(len(terminal_tasks)))
    terminal_summary = (
        len(terminals),
        sum(row[4] for row in terminals),
        sum(row[7] for row in terminals),
        sum(row[8] for row in terminals),
        sum(row[9] for row in terminals),
        sum(row[10] for row in terminals),
        sum(row[11] for row in terminals),
        sum(row[12] for row in terminals),
        sum(row[14] for row in terminals),
        min(row[16] for row in terminals),
        all(row[19] for row in terminals),
        all(row[20] for row in terminals),
    )
    require(terminal_summary == EXPECTED_TERMINAL, terminal_summary)
    terminal_semantics = tuple(terminal_semantic(row) for row in terminals)
    require(
        hashlib.sha256(repr(terminal_semantics).encode()).hexdigest()
        == EXPECTED_TERMINAL_SEMANTIC_SHA256,
        "terminal semantic",
    )
    closure_rows = tuple(closure_row(row) for row in terminals)
    closure_packet = ("z225-terminal-closure-solver-witness-free-v2", closure_rows)
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_vector_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_vector_sha == EXPECTED_CASE_VECTOR_SHA256, case_vector_sha)

    second_rows = atlas_rows(base, SECOND_LEVEL)
    second_census = (
        len(second_rows),
        sum(row[3] for row in second_rows),
        sum(not row[3] for row in second_rows),
    )
    require(second_census == EXPECTED_SECOND_CENSUS, second_census)
    require(
        hashlib.sha256(repr(second_rows).encode()).hexdigest()
        == EXPECTED_SECOND_ROW_ORDER_SHA256,
        "second row order",
    )
    second_tasks = tuple(
        (SECOND_LEVEL, body, ruler, high, wall)
        for body, ruler, high, wall in second_rows
    )
    if args.processes == 1:
        second_pairs = tuple(screen_worker(task) for task in second_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(second_tasks))) as pool:
            second_pairs = tuple(pool.map(screen_worker, second_tasks))
    second_by_task = {task: row for task, row in second_pairs}
    require(len(second_by_task) == len(second_tasks), "lost second screen task")
    second_screened = tuple(second_by_task[task] for task in second_tasks)
    second_totals = tuple(
        sum(row[index] for row in second_screened)
        for index in (9, 10, 11, 12)
    )
    require(second_totals == EXPECTED_SECOND_SCREEN, second_totals)
    second_farkas = (
        sum(row[19] for row in second_screened),
        sum(row[20] for row in second_screened),
    )
    require(second_farkas == EXPECTED_SECOND_FARKAS, second_farkas)
    # This layer contains three wall rows and one order row.  The atlas census
    # above verifies that split; every row still needs its exact status
    # certificate checked, but the order row must not carry the wall flag.
    require(
        all(row[16] == row[11] for row in second_screened),
        "second status control",
    )
    second_canonical_screen = tuple(row[:19] for row in second_screened)
    second_screen_sha = hashlib.sha256(
        repr(second_canonical_screen).encode()
    ).hexdigest()
    require(second_screen_sha == EXPECTED_SECOND_SCREEN_SHA256, second_screen_sha)
    require(all(row[12] == 0 for row in second_screened), "second residual")

    next_rows = atlas_rows(base, NEXT_LEVEL)
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(next_census == NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(next_rows).encode()).hexdigest() == NEXT_ROW_ORDER_SHA256, "next row order")
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z225-terminal-z224-screen-double-layer-cap223-v1",
        SOURCE_3114_SHA256,
        OUTPUT_3114_SHA256,
        SEMANTIC_3114_SHA256,
        base.thm.ATLAS_SHA256,
        (EXPECTED_ROW_ORDER_SHA256, census, totals, farkas, screen_sha),
        (EXPECTED_TERMINAL, closure_sha, case_vector_sha),
        (
            EXPECTED_SECOND_ROW_ORDER_SHA256,
            second_census,
            second_totals,
            second_farkas,
            second_screen_sha,
        ),
        LEDGER_BEFORE,
        LEDGER_AFTER,
        NEXT_CAP,
        (NEXT_LEVEL, next_census, NEXT_ROW_ORDER_SHA256),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z225 terminal plus z224 screen double-layer descent",
        f"dependency=THM3114_source:{SOURCE_3114_SHA256};output:{OUTPUT_3114_SHA256};semantic:{SEMANTIC_3114_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};residual_rows:{len(residual_rows)}",
        f"terminal=z1:{LEVEL};rows:{terminal_summary[0]};masks:{terminal_summary[1]};positive_two_high_gap:{sum(row[19] for row in terminals)};closed:{sum(row[20] for row in terminals)};zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};closure_sha256:{closure_sha};case_vector_sha256:{case_vector_sha}",
        f"second_screen=z1:{SECOND_LEVEL};rows:{second_census[0]};wall:{second_census[1]};order:{second_census[2]};row_order_sha256:{EXPECTED_SECOND_ROW_ORDER_SHA256};states:{second_totals[0]};crude:{second_totals[1]};status:{second_totals[2]};residual:{second_totals[3]};direct_farkas:{second_farkas[0]};legacy_farkas:{second_farkas[1]};screen_record_sha256:{second_screen_sha}",
    ]
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};exact={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment",
            "direction_terminal=the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "direction_second_screen=the_same_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_z224_projected_assignment;all_four_rows_have_zero_residual",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_each_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;the_top_level_semantic_binds_only_these_solver_witness_free_digests_and_basis_invariant_counts",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied_prime_full_order_handoff",
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
