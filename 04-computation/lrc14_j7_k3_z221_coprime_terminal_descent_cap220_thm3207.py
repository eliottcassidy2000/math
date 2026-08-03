#!/usr/bin/env python3
"""Exact projected-k3 z221 coprime screen and terminal descent to cap220."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3174 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py"
)
OUTPUT_3174 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.out"
)

SOURCE_3174_SHA256 = "3ac95e58861078828671e606e1c97705ac4def2728019ac21bd78eba0b9f1c18"
OUTPUT_3174_SHA256 = "4edc6c7efb64ce1062aa863a5f4c17e4735f42bf7b4ac371c6bca660346abe43"
SEMANTIC_3174_SHA256 = "f6493459d78f6af58b999bd54c8b72e2a8f989c533a1e5a91d5e6bc337352ec7"

LEVEL = 221
EXPECTED_CENSUS = (90, 83, 7)
EXPECTED_ROW_ORDER_SHA256 = "ab4f5b7fe3b5e1330d6e51dfb150527fa27cac56b755ed5e2ae2522f0f27ceb4"
EXPECTED_SCREEN = (8705, 3898, 4420, 387)
EXPECTED_FARKAS = (0, 4420)
EXPECTED_SCREEN_SHA256 = "820ff7d9bc0d593434b672f4bbdbd08d7c7a164b063503e2e0a98d793ce5f82a"
EXPECTED_RESIDUAL_SHA256 = "c8b5bb5bcda2034f95968022374e5b455c73e487f81873624a0a4d4c9020f38d"
EXPECTED_RESIDUAL_INDEX_SHA256 = "8f778fe45cf89af0ad1e441d996e2c77fd754b1bc61ab72a150ab5711060fccd"
EXPECTED_TERMINAL = (22, 387, 319, 392, 45, 318, 74, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "93e09fce6debb80536e4c9f40e3c6f6cb79f993f542cd939884757794ac4e924"
EXPECTED_CLOSURE_SHA256 = "6623b2b236e2fb96f003e13e549e66d621b9ad9596c3bcf182ab7dc03b71df03"
EXPECTED_CASE_VECTOR_SHA256 = "6c0bb93eee6fb06bed9ff192b54d3eb7a3d508161d0b0313e3917834ce18fb7b"
EXPECTED_SEMANTIC_SHA256 = "aad27ae040933935c6eacfd2801cef1644df3de3a8f031aebf7d1b238ac35e74"

LEDGER_BEFORE = 373806
LAYER_ROWS = 90
LEDGER_AFTER = 373716
NEXT_CAP = 220
NEXT_LEVEL = 220
NEXT_CENSUS = (289, 249, 40)
NEXT_ROW_ORDER_SHA256 = "bc75fef8aebdc6b350957ea96842b38a24bad2174200646ebaf2817844515799"

RESIDUAL = {
    (1, 4, 5, 9, 11, 14): (14, "0612bfd6334b896080174e3f633966df7c0ba1610b8d4a165acdda240b219959"),
    (1, 4, 10, 11, 12, 14): (19, "6961f6de4513c05009609877c88eb1f0ef9a22d15c2ea43ab39aca9c2b12ffd4"),
    (1, 5, 6, 9, 12, 14): (3, "6203d6d71cedb53284763e70b2cf82fbdb17288eac6c3c5c6eff7ccbb9b40a6a"),
    (1, 5, 9, 11, 12, 14): (26, "b4ed40d958f3692bf6bbbd8b0122c8770ee7da4dff0016f90a2743329bed4cfb"),
    (1, 6, 7, 9, 10, 12): (1, "92cb1a4f703fcb8557cc14eedb215e43d960d8c22d9fe2b2d2eef2ec76262ee1"),
    (1, 6, 8, 10, 11, 14): (14, "3ad46d2ebc96664252de030dc30fed1a1f9d9c9254e9bab66a23cd36dfcdc350"),
    (1, 6, 9, 10, 12, 14): (1, "92cb1a4f703fcb8557cc14eedb215e43d960d8c22d9fe2b2d2eef2ec76262ee1"),
    (1, 8, 10, 11, 12, 14): (121, "88bea840ce316ce78d7a268a9d9935fee33408c06baf8dffd6fe513368d262d2"),
    (1, 9, 10, 11, 12, 14): (59, "291eca59ffd2704800bb5bc6bd346110aa933fb953ee913fe8d877fd6141a67a"),
    (2, 3, 8, 10, 11, 14): (28, "37cf38def76bd82ab3acb21c35cb74c253457def3d673f3fe7497712b3f832f4"),
    (2, 5, 9, 11, 12, 14): (3, "2bef8aea8bb1d91721454f7cc7d2adda97b5ab2e47385d9669d0bba1b884e0d5"),
    (2, 6, 8, 10, 11, 14): (7, "6a024800c5993df185e6c9087b0bc54e651f41fe0df46ca564079eae83b082b6"),
    (2, 6, 8, 10, 12, 14): (1, "da13fe1fc019d89a16fc88a6e281d249673eb4f07475243634be873613648ff7"),
    (2, 6, 10, 11, 12, 14): (2, "b783f164c6da8760d8f3eb801a6b78a4e186482a8b91df03c1061cfcbb1f78e8"),
    (2, 8, 10, 11, 12, 14): (38, "d0e9f8635b4fb45829417a8abbd176c3eaba6ca23624d0c039e25338518298ef"),
    (2, 9, 10, 11, 12, 14): (37, "8bab55fb4fd2c25fc788631e9cc4245a20361b34ad66b18c5a79a00d027c63df"),
    (4, 5, 6, 9, 11, 14): (1, "e3e9ac257e24e5b81fcafee67b6829a990dd17f451dce91904d871e674243399"),
    (4, 6, 8, 10, 11, 14): (2, "e76b6ca6d2faa5bbaea9c9677ab14bcdec93e949165195a5554cf968d54f8daa"),
    (4, 6, 9, 10, 11, 14): (1, "3cca42015fdfd04e21d980d876d55e6330f6919081d5b9a4f3b8b509c07a896c"),
    (4, 8, 10, 11, 12, 14): (6, "47b84188aa3d63fdec07e61b952d741866bcf73a695bfd96ec3d775859d523f3"),
    (4, 9, 10, 11, 12, 14): (2, "aaf98e72b92248c70315e2326fc86160346dde26f514361a55526721c7d3bbd9"),
    (6, 8, 10, 11, 12, 14): (1, "0b820006d9374e7e90af0f5a5e864a6221a7edefe8d145a78d9631ebbd55ccbd"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3174):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3207_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    index, body, bank = task
    wrapper = load(f"thm3207_terminal_wrapper_{index}")
    entry = wrapper.load(f"thm3207_terminal_entry_{index}")
    source = entry.load(f"thm3207_terminal_source_{index}")
    driver = source.load(f"thm3207_terminal_driver_{index}")
    predecessor = driver.load(f"thm3207_terminal_predecessor_{index}")
    bridge = predecessor.load(f"thm3207_terminal_bridge_{index}")
    base = bridge.load(f"thm3207_terminal_base_{index}")
    return index, base.thm.terminal_probe((LEVEL, body, bank))


def modulus_data(task, row=None):
    level, body, ruler, high, wall = task
    require(level == LEVEL, level)
    require(ruler == 14 * lcm(*body), (body, ruler, "body ruler"))
    require(gcd(level, ruler) == 1, (body, ruler, "coprimality"))
    first_d = ruler
    if row is not None:
        require(row[:6] == (level, body, ruler, high, first_d, wall), (body, "header"))
        require(all(lcm(*ds) == ruler for ds in row[13]), (body, "residual modulus"))
    return first_d


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
    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    require(all(modulus_data(task) == task[2] for task in tasks), "fixed modulus")
    packet = (
        "lrc14-k3-z221-coprime-terminal-cap220-v1",
        (SOURCE_3174_SHA256, OUTPUT_3174_SHA256, SEMANTIC_3174_SHA256, base.thm.ATLAS_SHA256),
        (
            LEVEL,
            EXPECTED_CENSUS,
            EXPECTED_ROW_ORDER_SHA256,
            EXPECTED_SCREEN,
            EXPECTED_FARKAS,
            EXPECTED_SCREEN_SHA256,
            EXPECTED_RESIDUAL_SHA256,
            EXPECTED_RESIDUAL_INDEX_SHA256,
            tuple(RESIDUAL.items()),
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
    require(sha(SOURCE_3174) == SOURCE_3174_SHA256, "THM-3174 source changed")
    require(sha(OUTPUT_3174) == OUTPUT_3174_SHA256, "THM-3174 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3174_SHA256}" in OUTPUT_3174.read_text(encoding="utf-8"),
        "THM-3174 semantic changed",
    )
    residual_index_sha = hashlib.sha256(repr(tuple(RESIDUAL.items())).encode()).hexdigest()
    if EXPECTED_RESIDUAL_INDEX_SHA256 is not None:
        require(residual_index_sha == EXPECTED_RESIDUAL_INDEX_SHA256, "residual index changed")

    wrapper = load("thm3207_wrapper")
    entry = wrapper.load("thm3207_entry")
    source = entry.load("thm3207_source")
    driver = source.load("thm3207_driver")
    predecessor = driver.load("thm3207_predecessor")
    bridge = predecessor.load("thm3207_bridge")
    base = bridge.load("thm3207_base")
    rows = entry.atlas_rows(base, LEVEL)
    next_rows = entry.atlas_rows(base, NEXT_LEVEL)
    census, next_census, semantic = metadata(base, rows, next_rows)
    if args.metadata_only:
        print("SEMANTIC_SHA256", semantic)
        return
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if args.processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost screen task")
    screened = tuple(by_task[task] for task in tasks)
    for task, row in zip(tasks, screened):
        modulus_data(task, row)
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
    terminal_semantic_sha = hashlib.sha256(repr(terminal_semantics).encode()).hexdigest()
    require(terminal_semantic_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256, terminal_semantic_sha)
    closure_packet = (
        "z221-terminal-closure-solver-witness-free-v1",
        tuple(closure_row(row) for row in terminals),
    )
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    if EXPECTED_CLOSURE_SHA256 is not None:
        require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_vector_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_vector_sha == EXPECTED_CASE_VECTOR_SHA256, case_vector_sha)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    lines = [
        "LRC14 projected k3 z221 coprime terminal descent and cap220",
        f"dependency=THM3174_source:{SOURCE_3174_SHA256};output:{OUTPUT_3174_SHA256};semantic:{SEMANTIC_3174_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};residual_rows:{len(residual_rows)};residual_sha256:{residual_sha};residual_index_sha256:{residual_index_sha}",
        f"fixed_modulus=gcd(221,L):1;first_d_equals_L:1;every_residual_lcm_equals_L:1;order_rows:7;order_residual:0",
        f"terminal=z1:{LEVEL};rows:{terminal_summary[0]};masks:{terminal_summary[1]};positive_two_high_gap:{sum(row[19] for row in terminals)};closed:{sum(row[20] for row in terminals)};zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};terminal_semantic_sha256:{terminal_semantic_sha};closure_sha256:{closure_sha};case_vector_sha256:{case_vector_sha}",
    ]
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};exact={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;all_seven_order_rows_have_zero_residual",
            "direction_modulus=gcd(221,L)=1_for_all_90_rows_so_first_d=L;every_residual_status_tuple_has_lcm(ds)=L;the_z222_Boolean_divisor_square_has_disappeared",
            "direction_terminal=each_residual_is_a_wall_row;the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_the_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;all_22_residual_mask_banks_are_bound_individually",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied",
            "scope=projected_k3_necessary_atlas_only;coprimality_does_not_itself_close_wall_rows;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
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
