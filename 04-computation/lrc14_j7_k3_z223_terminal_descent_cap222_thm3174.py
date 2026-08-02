#!/usr/bin/env python3
"""Exact projected-k3 z223 screen and terminal descent to cap222."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3139 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
)
OUTPUT_3139 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
)
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.out"

SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"

LEVEL = 223
EXPECTED_CENSUS = (65, 54, 11)
EXPECTED_ROW_ORDER_SHA256 = "58454a33efe0b74fddeddb498d1518033ccf74064b19ec611d63a223f6065566"
EXPECTED_SCREEN = (2986, 1094, 1767, 125)
EXPECTED_FARKAS = (0, 1767)
EXPECTED_SCREEN_SHA256 = "88220270f2b1b43aeb0fe237472b31e971ad0a275ef904301f08a89cd748d2eb"
EXPECTED_RESIDUAL_SHA256 = "822fcd1cd347acbec8773b95a04bb5e17c598a8225cc9ebc4ab6ce7aeca96b83"
EXPECTED_TERMINAL = (7, 125, 115, 125, 13, 114, 11, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "54ca859fd914c26f52b78b93df89c29b2fe85b13983935de372d63fb847ec824"
EXPECTED_CLOSURE_SHA256 = "35b165deda53bd549fec3cf6ad4b70ce8f7c02b26c639d386cde8f80b6dfe071"
EXPECTED_CASE_VECTOR_SHA256 = "11ba29c9a3f933a495709a4449146d9fbce064b9978160724bdf324a52623774"
EXPECTED_SEMANTIC_SHA256 = "f6493459d78f6af58b999bd54c8b72e2a8f989c533a1e5a91d5e6bc337352ec7"

LEDGER_BEFORE = 374090
LAYER_ROWS = 65
LEDGER_AFTER = 374025
NEXT_CAP = 222
NEXT_LEVEL = 222
NEXT_CENSUS = (219, 199, 20)
NEXT_ROW_ORDER_SHA256 = "68ea5b7128d3deafd4bb9a14f6d1d09ba867a5f117645a3e1f4f3d03e761de73"

RESIDUAL = {
    (1, 5, 6, 9, 12, 14): (3, "6203d6d71cedb53284763e70b2cf82fbdb17288eac6c3c5c6eff7ccbb9b40a6a"),
    (1, 5, 9, 11, 12, 14): (26, "b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb"),
    (1, 6, 8, 10, 13, 14): (21, "d056b140d4b86a0ccca4e0ffa0027a7aacf47bc3936b2ad2d264e53b96d21d92"),
    (1, 8, 10, 11, 12, 14): (35, "a3f1566c773889925e5bf327a4d40da538a1cb4054b43895fecb3f71c0baba9c"),
    (1, 8, 10, 12, 13, 14): (29, "640d7f8226f73a963de2e88d78cce8e984533b8f64fd073a9c9fe7ed2ed7f408"),
    (2, 6, 8, 10, 12, 14): (1, "da13fe1fc019d89a16fc88a6e281d249673eb4f07475243634be873613648ff7"),
    (2, 8, 10, 11, 12, 14): (10, "00415b63ca4ba804c327c3ebb302ddbaec054ca527ca8eb5195e8c66adc8c5ac"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3139):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3174_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    index, body, bank = task
    entry = load(f"thm3174_terminal_entry_{index}")
    source = entry.load(f"thm3174_terminal_source_{index}")
    driver = source.load(f"thm3174_terminal_driver_{index}")
    predecessor = driver.load(f"thm3174_terminal_predecessor_{index}")
    bridge = predecessor.load(f"thm3174_terminal_bridge_{index}")
    base = bridge.load(f"thm3174_terminal_base_{index}")
    return index, base.thm.terminal_probe((LEVEL, body, bank))


def fixed_modulus(task, row=None):
    level, body, ruler, high, wall = task
    require(level == LEVEL, level)
    require(ruler == 14 * lcm(*body), (body, ruler))
    require(gcd(level, ruler) == 1, (body, ruler, "coprimality"))
    first_d = ruler // gcd(level, ruler)
    require(first_d == ruler, (body, first_d, ruler))
    if row is None:
        return
    require(row[:6] == (level, body, ruler, high, first_d, wall), (body, "header"))
    require(all(lcm(*ds) == ruler for ds in row[13]), (body, "residual modulus"))


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


def metadata(entry, base, rows, next_rows):
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, "row order")
    require(hashlib.sha256(repr(next_rows).encode()).hexdigest() == NEXT_ROW_ORDER_SHA256, "next row order")
    packet = (
        "lrc14-k3-z223-terminal-cap222-v1",
        SOURCE_3139_SHA256,
        OUTPUT_3139_SHA256,
        SEMANTIC_3139_SHA256,
        base.thm.ATLAS_SHA256,
        (
            EXPECTED_ROW_ORDER_SHA256,
            EXPECTED_CENSUS,
            EXPECTED_SCREEN,
            EXPECTED_FARKAS,
            EXPECTED_SCREEN_SHA256,
            EXPECTED_RESIDUAL_SHA256,
        ),
        (
            EXPECTED_TERMINAL,
            EXPECTED_TERMINAL_SEMANTIC_SHA256,
            EXPECTED_CLOSURE_SHA256,
            EXPECTED_CASE_VECTOR_SHA256,
        ),
        LEDGER_BEFORE,
        LEDGER_AFTER,
        NEXT_CAP,
        (NEXT_LEVEL, NEXT_CENSUS, NEXT_ROW_ORDER_SHA256),
    )
    return census, next_census, hashlib.sha256(repr(packet).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=12)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3139) == SOURCE_3139_SHA256, "THM-3139 source changed")
    require(sha(OUTPUT_3139) == OUTPUT_3139_SHA256, "THM-3139 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3139_SHA256}" in OUTPUT_3139.read_text(),
        "THM-3139 semantic changed",
    )

    entry = load("thm3174_entry")
    source = entry.load("thm3174_source")
    driver = source.load("thm3174_driver")
    predecessor = driver.load("thm3174_predecessor")
    bridge = predecessor.load("thm3174_bridge")
    base = bridge.load("thm3174_base")
    rows = entry.atlas_rows(base, LEVEL)
    next_rows = entry.atlas_rows(base, NEXT_LEVEL)
    census, next_census, semantic = metadata(entry, base, rows, next_rows)
    if args.metadata_only:
        print("SEMANTIC_SHA256", semantic)
        return
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    for task in tasks:
        fixed_modulus(task)
    if args.processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost screen task")
    screened = tuple(by_task[task] for task in tasks)
    for task, row in zip(tasks, screened):
        fixed_modulus(task, row)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "status control")
    canonical_screen = tuple(row[:19] for row in screened)
    screen_sha = hashlib.sha256(repr(canonical_screen).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual_rows = tuple(row for row in screened if row[12])
    residual = tuple((row[1], row[13]) for row in residual_rows)
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(residual_sha == EXPECTED_RESIDUAL_SHA256, residual_sha)
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
    closure_packet = (
        "z223-terminal-closure-solver-witness-free-v1",
        tuple(closure_row(row) for row in terminals),
    )
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_vector_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_vector_sha == EXPECTED_CASE_VECTOR_SHA256, case_vector_sha)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    lines = [
        "LRC14 projected k3 z223 terminal descent and cap222",
        f"dependency=THM3139_source:{SOURCE_3139_SHA256};output:{OUTPUT_3139_SHA256};semantic:{SEMANTIC_3139_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};residual_rows:{len(residual_rows)};residual_sha256:{residual_sha}",
        f"terminal=z1:{LEVEL};rows:{terminal_summary[0]};masks:{terminal_summary[1]};positive_two_high_gap:{sum(row[19] for row in terminals)};closed:{sum(row[20] for row in terminals)};zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};closure_sha256:{closure_sha};case_vector_sha256:{case_vector_sha}",
    ]
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};exact={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;all_order_rows_have_zero_residual",
            "direction_terminal=each_residual_is_a_wall_row;the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "fixed_modulus=gcd(223,L)=1_for_every_row_so_first_d=L;every_residual_status_tuple_has_lcm(ds)=L",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_the_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;the_top_level_semantic_binds_only_solver_witness_free_digests_and_basis_invariant_counts",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied",
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
