#!/usr/bin/env python3
"""Exact projected-k3 z219 common-gcd-three screen and terminal descent."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3218 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.py"
)
OUTPUT_3218 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.out"
)

SOURCE_3218_SHA256 = "97c55e0a9acc9b42c7c7a856889bd2b8fbd49854ac39a1fe974a488f63a3cff8"
OUTPUT_3218_SHA256 = "af7fe69ba68579611f171b9b2252f9a9923bdd2785720a6cb1b208a611ee0565"
SEMANTIC_3218_SHA256 = "ed67ff55840d07965a6cf8a478c0fc49cab637511ad3d684b2de8831dfcfcd94"

LEVEL = 219
EXPECTED_CENSUS = (16, 15, 1)
EXPECTED_ROW_ORDER_SHA256 = "6a39c83d47b8080fd52e64479b886f1b6a887150f32cfbf30a676fbb0aaa54ba"
EXPECTED_G_STRATA = ((3, 16, 15, 1),)
EXPECTED_SCREEN = (2107, 1022, 661, 424)
EXPECTED_WALL_SCREEN = (2106, 1021, 661, 424)
EXPECTED_ORDER_SCREEN = (1, 1, 0, 0)
EXPECTED_FARKAS = (0, 661)
EXPECTED_WALL_FARKAS = (0, 661)
EXPECTED_ORDER_FARKAS = (0, 0)
EXPECTED_SCREEN_SHA256 = "868d5cce16d5f704f783987b4ecba4e55fda009bce1d48611bb8f3dec647e248"
EXPECTED_RESIDUAL_SHA256 = "69918f14f67d818749581af98f8671d7567da951bf0be52e3274fa558a186881"
EXPECTED_RESIDUAL_INDEX_SHA256 = "d5fbefd0839e5eb1f47b42c8c41b02bfd5073ae84fc91354d6c2d754460c8157"
EXPECTED_QUOTIENT_COUNTS = ((3, 3, 424),)
EXPECTED_TERMINAL = (10, 424, 413, 424, 13, 411, 13, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "ae07780990ec507a6a3ff883363dd223af4f11cafdc48543d8e5f616c9896049"
EXPECTED_CLOSURE_SHA256 = "14b8e92a82cfde0daddde1859237f0ec724858a9b993e51208fae1d1b18a0641"
EXPECTED_CASE_VECTOR_SHA256 = "127d5e241f7c6aa4493b158a316890b9d8455205f793dafbdd56725de2c516f5"
EXPECTED_MIN_GAP = ((1, 5, 7, 11, 12, 13), 176513570591453, 67539934271289480)
EXPECTED_SEMANTIC_SHA256 = "999e0186c706441074e50ca1a2e0d689a9f19c6fa2c13137abc730eb51f92a49"

LEDGER_BEFORE = 373427
LAYER_ROWS = 16
LEDGER_AFTER = 373411
NEXT_CAP = 218
NEXT_LEVEL = 218
NEXT_CENSUS = (119, 117, 2)
NEXT_ROW_ORDER_SHA256 = "d416b484148e008e9f58273f8533cc9d05a171e5ba4963214d3fc73d7fc20bc7"

RESIDUAL = {
    (1, 5, 7, 9, 12, 13): (24, "451a4d0903a0b4c390e6bb832635a67291db39cfa551f0b19a9ececa25ceb47d"),
    (1, 5, 7, 11, 12, 13): (45, "28c7a040cb15b60967761f7bb1f0ccb5ec2e46ce83b13752808d4cd5d5c3c4a5"),
    (1, 5, 9, 11, 12, 14): (48, "1bdeb0ce2e6585cda0312e29b8a2ca3ec95e08a700da594a46fe804a8bbea1ef"),
    (1, 7, 9, 11, 12, 13): (24, "4373397faf1b72e6f954fa42c267c445ae9f2df02d368762fcf8ef5c22cd6589"),
    (1, 8, 10, 12, 13, 14): (18, "ae3a12b2d7b3879f1ab20968616d86adf5373ddc5a35c46e4555070d03b32177"),
    (1, 9, 10, 11, 12, 14): (59, "a0dfcb2dd55f7eeb02cd97e2e36b045568e0df48480bd697624d229fc41712fd"),
    (2, 5, 7, 9, 11, 13): (85, "f37f19e9c2f077467930812e3388e0d028b78b5af939eda310dae5b4f72cec2a"),
    (2, 5, 9, 11, 12, 14): (24, "18588ef4ffe140a6fd053612c91ac51118e320984c912e224c0b6e292ea38bd3"),
    (2, 5, 9, 11, 13, 14): (77, "0acc575c1f7da65a5351532cb540d20b2a1f9fda2e07b4325750fb7594c9dffd"),
    (2, 5, 9, 12, 13, 14): (20, "3f643d26a2ba3a13f0f87fae140c245b907db3b0d9826820e647b5d2df979eca"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3218):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3230_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    successor = load(f"thm3230_terminal_successor_{task[0]}")
    successor.LEVEL = LEVEL
    return successor.terminal_worker(task)


def modulus_data(task, row=None):
    level, body, ruler, high, wall = task
    require(level == LEVEL, level)
    require(ruler == 14 * lcm(*body), (body, ruler, "body ruler"))
    g = gcd(level, ruler)
    require(g == 3, (body, ruler, g))
    first_d = ruler // g
    if row is not None:
        require(row[:6] == (level, body, ruler, high, first_d, wall), (body, "header"))
        for ds in row[13]:
            d = lcm(*ds)
            require(d % first_d == 0, (body, ds, first_d))
            q = d // first_d
            require(q == 3, (body, g, q))
    return g, first_d


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


def component_counts(base):
    pattern = re.compile(
        r"^row=E=([0-9,]+);h=[^;]+;r=([0-9]+);"
        r"L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    counts = []
    for line in base.thm.ATLAS.read_text(encoding="utf-8").splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(5)) == LEVEL:
            counts.append(int(match.group(2)))
    require(len(counts) == EXPECTED_CENSUS[0], len(counts))
    return tuple(counts)


def heavy_first(tasks, counts):
    require(len(tasks) == len(counts), (len(tasks), len(counts)))
    indexed = tuple(enumerate(tasks))
    return tuple(
        task
        for index, task in sorted(
            indexed,
            key=lambda item: (-item[1][2] * counts[item[0]], -item[1][2], item[0]),
        )
    )


def metadata(base, rows, next_rows):
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, "row order")
    require(
        hashlib.sha256(repr(next_rows).encode()).hexdigest() == NEXT_ROW_ORDER_SHA256,
        "next row order",
    )
    strata = tuple(
        (
            g,
            sum(gcd(LEVEL, row[1]) == g for row in rows),
            sum(gcd(LEVEL, row[1]) == g and row[3] for row in rows),
            sum(gcd(LEVEL, row[1]) == g and not row[3] for row in rows),
        )
        for g in (3,)
    )
    require(strata == EXPECTED_G_STRATA, strata)
    packet = (
        "lrc14-k3-z219-common-gcd-three-terminal-cap218-v1",
        (SOURCE_3218_SHA256, OUTPUT_3218_SHA256, SEMANTIC_3218_SHA256, base.thm.ATLAS_SHA256),
        (
            LEVEL,
            EXPECTED_CENSUS,
            EXPECTED_ROW_ORDER_SHA256,
            EXPECTED_G_STRATA,
            EXPECTED_SCREEN,
            EXPECTED_WALL_SCREEN,
            EXPECTED_ORDER_SCREEN,
            EXPECTED_FARKAS,
            EXPECTED_WALL_FARKAS,
            EXPECTED_ORDER_FARKAS,
            EXPECTED_SCREEN_SHA256,
            EXPECTED_RESIDUAL_SHA256,
            EXPECTED_RESIDUAL_INDEX_SHA256,
            EXPECTED_QUOTIENT_COUNTS,
            tuple(RESIDUAL.items()),
        ),
        (
            EXPECTED_TERMINAL,
            EXPECTED_TERMINAL_SEMANTIC_SHA256,
            EXPECTED_CLOSURE_SHA256,
            EXPECTED_CASE_VECTOR_SHA256,
            EXPECTED_MIN_GAP,
        ),
        LEDGER_BEFORE,
        LEDGER_AFTER,
        NEXT_CAP,
        (NEXT_LEVEL, NEXT_CENSUS, NEXT_ROW_ORDER_SHA256),
    )
    return census, next_census, strata, hashlib.sha256(repr(packet).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=12)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3218) == SOURCE_3218_SHA256, "THM-3218 source changed")
    require(sha(OUTPUT_3218) == OUTPUT_3218_SHA256, "THM-3218 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3218_SHA256}" in OUTPUT_3218.read_text(encoding="utf-8"),
        "THM-3218 semantic changed",
    )
    residual_index_sha = hashlib.sha256(repr(tuple(RESIDUAL.items())).encode()).hexdigest()
    require(residual_index_sha == EXPECTED_RESIDUAL_INDEX_SHA256, residual_index_sha)

    successor = load("thm3230_successor")
    canonical = successor.load("thm3230_canonical")
    wrapper = canonical.load("thm3230_wrapper")
    entry = wrapper.load("thm3230_entry")
    source = entry.load("thm3230_source")
    driver = source.load("thm3230_driver")
    predecessor = driver.load("thm3230_predecessor")
    bridge = predecessor.load("thm3230_bridge")
    base = bridge.load("thm3230_base")
    rows = entry.atlas_rows(base, LEVEL)
    next_rows = entry.atlas_rows(base, NEXT_LEVEL)
    census, next_census, strata, semantic = metadata(base, rows, next_rows)
    if args.metadata_only:
        print("SEMANTIC_SHA256", semantic)
        return
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if args.processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        scheduled = heavy_first(tasks, component_counts(base))
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            pairs = tuple(pool.imap_unordered(screen_worker, scheduled, chunksize=1))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost screen task")
    screened = tuple(by_task[task] for task in tasks)
    for task, row in zip(tasks, screened):
        modulus_data(task, row)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    wall_rows = tuple(row for row in screened if row[5])
    order_rows = tuple(row for row in screened if not row[5])
    wall_totals = tuple(sum(row[index] for row in wall_rows) for index in (9, 10, 11, 12))
    order_totals = tuple(sum(row[index] for row in order_rows) for index in (9, 10, 11, 12))
    require(wall_totals == EXPECTED_WALL_SCREEN, wall_totals)
    require(order_totals == EXPECTED_ORDER_SCREEN, order_totals)
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    wall_farkas = (sum(row[19] for row in wall_rows), sum(row[20] for row in wall_rows))
    order_farkas = (sum(row[19] for row in order_rows), sum(row[20] for row in order_rows))
    require(farkas == EXPECTED_FARKAS, farkas)
    require(wall_farkas == EXPECTED_WALL_FARKAS, wall_farkas)
    require(order_farkas == EXPECTED_ORDER_FARKAS, order_farkas)
    require(all(row[16] == row[11] for row in screened), "status control")
    canonical_screen = tuple(row[:19] for row in screened)
    screen_sha = hashlib.sha256(repr(canonical_screen).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual_rows = tuple(row for row in screened if row[12])
    require(all(row[5] for row in residual_rows), "order residual")
    residual = tuple((row[1], row[13]) for row in residual_rows)
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(residual_sha == EXPECTED_RESIDUAL_SHA256, residual_sha)
    require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), "residual bodies")
    quotient_counts = Counter()
    for row in residual_rows:
        count, mask_sha = RESIDUAL[row[1]]
        require(row[12] == count, (row[1], row[12]))
        require(hashlib.sha256(repr(row[13]).encode()).hexdigest() == mask_sha, (row[1], "mask sha"))
        g = gcd(LEVEL, row[2])
        for ds in row[13]:
            quotient_counts[(g, lcm(*ds) // row[4])] += 1
    quotient_packet = tuple((g, q, count) for (g, q), count in sorted(quotient_counts.items()))
    require(quotient_packet == EXPECTED_QUOTIENT_COUNTS, quotient_packet)

    terminal_tasks = tuple((index, row[1], row[13]) for index, row in enumerate(residual_rows))
    if args.processes == 1:
        terminal_pairs = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(terminal_tasks))) as pool:
            terminal_pairs = tuple(pool.imap_unordered(terminal_worker, terminal_tasks, chunksize=1))
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
    terminal_semantic_sha = hashlib.sha256(repr(terminal_semantics).encode()).hexdigest()
    require(terminal_semantic_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256, terminal_semantic_sha)
    closure_packet = (
        "z219-terminal-closure-solver-witness-free-v1",
        tuple(closure_row(row) for row in terminals),
    )
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_vector_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_vector_sha == EXPECTED_CASE_VECTOR_SHA256, case_vector_sha)
    min_terminal = min(terminals, key=lambda row: row[5])
    min_gap = (min_terminal[1], min_terminal[5].numerator, min_terminal[5].denominator)
    require(min_gap == EXPECTED_MIN_GAP, min_gap)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    lines = [
        "LRC14 projected k3 z219 common-gcd-three terminal descent and cap218",
        f"dependency=THM3218_source:{SOURCE_3218_SHA256};output:{OUTPUT_3218_SHA256};semantic:{SEMANTIC_3218_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};residual_rows:{len(residual_rows)};residual_sha256:{residual_sha};residual_index_sha256:{residual_index_sha}",
        f"wall_screen=states:{wall_totals[0]};crude:{wall_totals[1]};status:{wall_totals[2]};residual:{wall_totals[3]};direct_farkas:{wall_farkas[0]};legacy_farkas:{wall_farkas[1]}",
        f"order_screen=states:{order_totals[0]};crude:{order_totals[1]};status:{order_totals[2]};residual:{order_totals[3]};direct_farkas:{order_farkas[0]};legacy_farkas:{order_farkas[1]}",
        f"valuation_strata={strata};first_d_equals_L_over_g:1;residual_quotient_counts:{quotient_packet}",
        f"terminal=z1:{LEVEL};rows:{terminal_summary[0]};masks:{terminal_summary[1]};positive_two_high_gap:{sum(row[19] for row in terminals)};closed:{sum(row[20] for row in terminals)};zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};minimum_gap:{EXPECTED_MIN_GAP[1]}/{EXPECTED_MIN_GAP[2]};minimum_gap_body:{EXPECTED_MIN_GAP[0]};terminal_semantic_sha256:{terminal_semantic_sha};closure_sha256:{closure_sha};case_vector_sha256:{case_vector_sha}",
    ]
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};exact={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;the_unique_order_row_has_zero_residual",
            "direction_modulus=first_d=L/g_for_g=gcd(219,L)=3;every_residual_D_over_first_d_equals_3;this_is_one_occupied_pointed_divisor_vertex_not_a_C3_action_or_a_proved_common_carrier_transition",
            "direction_terminal=each_residual_is_a_wall_row;the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_the_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;all_10_residual_mask_banks_are_bound_individually",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied",
            "scope=projected_k3_necessary_atlas_only;the_common_gcd_three_fibre_is_a_pointed_divisor_condition_not_a_C3_torsor_or_PSL2_action;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
