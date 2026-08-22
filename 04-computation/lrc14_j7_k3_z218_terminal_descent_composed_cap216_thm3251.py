#!/usr/bin/env python3
"""Exact projected-k3 z218 terminal descent and composed cap216."""

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
SOURCE_3139 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
)
OUTPUT_3139 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out"
)
SOURCE_3230 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.py"
)
OUTPUT_3230 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.out"
)
SOURCE_3242 = ROOT / (
    "04-computation/lrc14_j7_k3_z217_exact_status_annihilation_thm3242.py"
)
OUTPUT_3242 = ROOT / (
    "05-knowledge/results/lrc14_j7_k3_z217_exact_status_annihilation_thm3242.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z218_terminal_descent_composed_cap216_thm3251.out"
)

SOURCE_3139_SHA256 = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
OUTPUT_3139_SHA256 = "4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5"
SEMANTIC_3139_SHA256 = "1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a"
SOURCE_3230_SHA256 = "0b14b1a86e2a08aae20433ee12d601128090f04658c0934fc6b263fb31d9e723"
OUTPUT_3230_SHA256 = "d8a4e87382da55879ee5f517f79277e4c3913823d5f968931ded22f84ace3bcf"
SEMANTIC_3230_SHA256 = "999e0186c706441074e50ca1a2e0d689a9f19c6fa2c13137abc730eb51f92a49"
SOURCE_3242_SHA256 = "165f89184856aae9c4d784e060123c3b639d232d135ff3c744cdfefb5580454a"
OUTPUT_3242_SHA256 = "10741cdf762edb9b370f91bbfe8739208d349d1d2cbaa37c5b68f3f1a37d34f3"
SEMANTIC_3242_SHA256 = "a28872885af6aef23f3aa00fe026b5f4c23b6d48302924fab37279cc2a4daf50"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 218
EXPECTED_CENSUS = (119, 117, 2)
EXPECTED_ROW_SHA256 = "d416b484148e008e9f58273f8533cc9d05a171e5ba4963214d3fc73d7fc20bc7"
EXPECTED_G_STRATA = ((2, 119, 117, 2),)
EXPECTED_TOTALS = (12833, 5375, 6522, 936)
EXPECTED_WALL = (12827, 5370, 6521, 936)
EXPECTED_ORDER = (6, 5, 1, 0)
EXPECTED_FARKAS = (0, 6522)
EXPECTED_WALL_FARKAS = (0, 6521)
EXPECTED_ORDER_FARKAS = (0, 1)
EXPECTED_SCREEN_SHA256 = "d89494bb6dca6e5328765782965c0f61d402f8eb05b580dbcc65a1ea665630c3"
EXPECTED_RESIDUAL_SHA256 = "e53f23cdce4d26addd82ab9a8321a031a70693c26f3b941ca6a0c86ce033790e"
EXPECTED_RESIDUAL_INDEX_SHA256 = "fa8cfed8aebd06d4f9382a29d8c026d4333465457f770c5ae081bfe6144b0cd7"
EXPECTED_QUOTIENT_COUNTS = ((2, 2, 936),)

EXPECTED_TERMINAL = (24, 936, 867, 1107, 63, 1031, 76, 0, 0, 1, 24, 24)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = (
    "18646b406cbe94aa45a1a3bd1abaad338c089c6bf1bf6459bee4df60b03325b7"
)
EXPECTED_CLOSURE_SHA256 = "6e09995ebdd03326120531c8d522efec6eb0483242303e6ba4d14cecdfde4272"
EXPECTED_CASE_VECTOR_SHA256 = "689e2709c2de478ffe8b137b7fb1c9f1b90eac9bb304dc0106a18f097ec76da8"
EXPECTED_MIN_GAP = ((1, 5, 6, 9, 12, 14), 730589, 522413892)
EXPECTED_SEMANTIC_SHA256 = "fe1f0e0df1c652be79b1e9d1e964bc4db64f1930e25b8fb8a6f7d3c5f47ccedb"

LEDGER_BEFORE = 373411
LAYER_ROWS = 119
INTERMEDIATE_LEDGER = 373292
INTERMEDIATE_CAP = 217
COMPOSED_ROWS = 8
LEDGER_AFTER = 373284
NEXT_CAP = 216
NEXT_LEVEL = 216
NEXT_CENSUS = (480, 447, 33)
NEXT_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"

# Ordered by the promoted atlas.  Each bank is bound by its independent mask hash.
RESIDUAL = (
    (8, (1, 2, 5, 9, 12, 14), 4, "b2ea24078384825da92ca7bfba74042ff3cf59c474e2900f4d86d169827434a7"),
    (23, (1, 3, 5, 9, 12, 14), 1, "85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a"),
    (32, (1, 4, 5, 6, 9, 14), 1, "85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a"),
    (34, (1, 4, 5, 9, 11, 14), 33, "f0c9c8b6c2a829316266541d275a37ac074aac4a8b4d4e0f60e03144073f795b"),
    (35, (1, 4, 5, 9, 12, 14), 1, "8bc196a9a680f512f95f366572b354132266820bfa047cf741a53845cd36a510"),
    (47, (1, 5, 6, 9, 11, 14), 2, "ae23d45b243e90361761d833aefcaa4c315ba2830558275ba4514872655c623f"),
    (48, (1, 5, 6, 9, 12, 14), 1, "85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a"),
    (51, (1, 5, 9, 10, 12, 14), 2, "79153270d573d2c59b4a80ef5f54f09a205e612c5cb053278a8525094b9efaa5"),
    (52, (1, 5, 9, 11, 12, 14), 217, "6fd60e94d3c706bb37f9bcd5d2fc6902a8fd3a25f95fb55339d87cb9b220d63d"),
    (61, (1, 7, 8, 10, 12, 13), 6, "015c36861016dcd443d6d3991c438e4558d7d046820e5960583f826980649b4f"),
    (63, (1, 8, 10, 11, 12, 14), 16, "ce36929944846eab8ed6fb77a132c5b571b7ba6f5a2302387f2685ec3f292f55"),
    (64, (1, 8, 10, 12, 13, 14), 16, "a5735baaedea8704219771659d029c6ef4063f3a5ba39e90166c6397d929a88d"),
    (65, (1, 9, 10, 11, 12, 14), 24, "fe196c407eb6d7d1fc5a65f283a8ce078e534602152521e8b1081d908ed725ea"),
    (71, (2, 3, 5, 9, 12, 14), 1, "85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a"),
    (74, (2, 3, 8, 10, 11, 14), 8, "fb3f066c58817a8ec0fbcc41b1b282f43cce2054279f600eb4400eea1b651ecb"),
    (79, (2, 4, 5, 9, 11, 14), 7, "366472439f6c971c7d89427f834d1d6b393f7d071f7f4c139405c29e190288d2"),
    (90, (2, 5, 9, 10, 12, 14), 2, "79153270d573d2c59b4a80ef5f54f09a205e612c5cb053278a8525094b9efaa5"),
    (91, (2, 5, 9, 11, 12, 14), 365, "4f47b092d529d7d685e24ad7644fd653f27424d47d7f19541660dfbbeab87abb"),
    (99, (2, 8, 10, 11, 12, 14), 2, "ea5092fa03963a603902fd1e857bf0d87f6ed908ff4577594f44f18ecec4574c"),
    (100, (2, 9, 10, 11, 12, 14), 81, "75bba523c3513b5383b0cd26808d32083dd7e64e29d4d57cb4c5b239a2eb08f4"),
    (104, (3, 5, 9, 11, 12, 14), 1, "344c16510bf13cfcf43e93edfe5840c4dceb813c8fa63d4bf411dcc8f46ceeb0"),
    (105, (3, 9, 10, 11, 12, 14), 1, "8f9fb2957e82b3458338dbd98f3818ff923539b9f722dbf03b437fd20d14abd0"),
    (107, (4, 5, 9, 11, 12, 14), 5, "2c83791d9f942b923e8793853de226b39e207407fa062203c93587c119fcbe1b"),
    (108, (4, 5, 9, 11, 13, 14), 139, "ab3a7d7bf776af96bcbc8f314e8f48f72f02fdec609a28266424405b488aa198"),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
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
    return load(f"thm3251_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    local_index, atlas_index, body, bank = task
    entry = load(f"thm3251_terminal_entry_{atlas_index}")
    source = entry.load(f"thm3251_terminal_source_{atlas_index}")
    driver = source.load(f"thm3251_terminal_driver_{atlas_index}")
    predecessor = driver.load(f"thm3251_terminal_predecessor_{atlas_index}")
    bridge = predecessor.load(f"thm3251_terminal_bridge_{atlas_index}")
    base = bridge.load(f"thm3251_terminal_base_{atlas_index}")
    row = base.thm.terminal_probe((LEVEL, body, bank))
    require(row[0] == LEVEL and row[1] == body and row[4] == len(bank), (atlas_index, row[:5]))
    return local_index, atlas_index, row


def terminal_semantic(row):
    # row[6] is a deterministic duplicate-gap maximizer witness.
    return (*row[:6], *row[7:])


def closure_row(row):
    return (
        row[1], row[4], (row[5].numerator, row[5].denominator), row[7], row[8],
        row[9], row[10], row[11], row[12], row[14], row[16], row[17], row[19], row[20],
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
    return tuple(
        task
        for index, task in sorted(
            enumerate(tasks),
            key=lambda item: (-item[1][2] * counts[item[0]], -item[1][2], item[0]),
        )
    )


def check_dependency(path, expected_sha, semantic=None):
    require(lf_sha(path) == expected_sha, (path, "hash changed"))
    if semantic is not None:
        require(
            f"semantic_sha256={semantic}" in path.read_text(encoding="utf-8"),
            (path, "semantic changed"),
        )


def metadata(base, rows, next_rows):
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows), sum(row[3] for row in next_rows), sum(not row[3] for row in next_rows)
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256, "row order")
    require(
        hashlib.sha256(repr(next_rows).encode()).hexdigest() == NEXT_ROW_SHA256,
        "next row order",
    )
    strata = tuple(
        (
            g,
            sum(gcd(LEVEL, row[1]) == g for row in rows),
            sum(gcd(LEVEL, row[1]) == g and row[3] for row in rows),
            sum(gcd(LEVEL, row[1]) == g and not row[3] for row in rows),
        )
        for g in (2,)
    )
    require(strata == EXPECTED_G_STRATA, strata)
    residual_index_sha = hashlib.sha256(repr(RESIDUAL).encode()).hexdigest()
    if EXPECTED_RESIDUAL_INDEX_SHA256 is not None:
        require(residual_index_sha == EXPECTED_RESIDUAL_INDEX_SHA256, residual_index_sha)
    packet = (
        "lrc14-k3-z218-terminal-composed-cap216-v1",
        (
            SOURCE_3139_SHA256, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256,
            SOURCE_3230_SHA256, OUTPUT_3230_SHA256, SEMANTIC_3230_SHA256,
            SOURCE_3242_SHA256, OUTPUT_3242_SHA256, SEMANTIC_3242_SHA256,
            ATLAS_SHA256,
        ),
        (
            LEVEL, EXPECTED_CENSUS, EXPECTED_ROW_SHA256, EXPECTED_G_STRATA,
            EXPECTED_TOTALS, EXPECTED_WALL, EXPECTED_ORDER, EXPECTED_FARKAS,
            EXPECTED_WALL_FARKAS, EXPECTED_ORDER_FARKAS, EXPECTED_SCREEN_SHA256,
            EXPECTED_RESIDUAL_SHA256, residual_index_sha, EXPECTED_QUOTIENT_COUNTS,
            RESIDUAL,
        ),
        (
            EXPECTED_TERMINAL, EXPECTED_TERMINAL_SEMANTIC_SHA256,
            EXPECTED_CLOSURE_SHA256, EXPECTED_CASE_VECTOR_SHA256, EXPECTED_MIN_GAP,
        ),
        (
            LEDGER_BEFORE, LAYER_ROWS, INTERMEDIATE_LEDGER, INTERMEDIATE_CAP,
            COMPOSED_ROWS, LEDGER_AFTER, NEXT_CAP,
        ),
        (NEXT_LEVEL, NEXT_CENSUS, NEXT_ROW_SHA256),
    )
    return census, next_census, strata, residual_index_sha, hashlib.sha256(repr(packet).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=12)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    check_dependency(SOURCE_3139, SOURCE_3139_SHA256)
    check_dependency(OUTPUT_3139, OUTPUT_3139_SHA256, SEMANTIC_3139_SHA256)
    check_dependency(SOURCE_3230, SOURCE_3230_SHA256)
    check_dependency(OUTPUT_3230, OUTPUT_3230_SHA256, SEMANTIC_3230_SHA256)
    check_dependency(SOURCE_3242, SOURCE_3242_SHA256)
    check_dependency(OUTPUT_3242, OUTPUT_3242_SHA256, SEMANTIC_3242_SHA256)

    entry = load("thm3251_entry")
    source = entry.load("thm3251_source")
    driver = source.load("thm3251_driver")
    predecessor = driver.load("thm3251_predecessor")
    bridge = predecessor.load("thm3251_bridge")
    base = bridge.load("thm3251_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    rows = entry.atlas_rows(base, LEVEL)
    next_rows = entry.atlas_rows(base, NEXT_LEVEL)
    census, next_census, strata, residual_index_sha, semantic = metadata(base, rows, next_rows)
    if args.metadata_only:
        print("RESIDUAL_INDEX_SHA256", residual_index_sha)
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
        level, body, ruler, high, wall = task
        require(row[:6] == (level, body, ruler, high, ruler // 2, wall), (body, "header"))
        require(ruler == 14 * lcm(*body), (body, ruler))
        require(gcd(level, ruler) == 2, (body, ruler))
        require(row[16] == row[11], (body, "status verification"))

    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    wall_rows = tuple(row for row in screened if row[5])
    order_rows = tuple(row for row in screened if not row[5])
    wall_totals = tuple(sum(row[index] for row in wall_rows) for index in (9, 10, 11, 12))
    order_totals = tuple(sum(row[index] for row in order_rows) for index in (9, 10, 11, 12))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    wall_farkas = (sum(row[19] for row in wall_rows), sum(row[20] for row in wall_rows))
    order_farkas = (sum(row[19] for row in order_rows), sum(row[20] for row in order_rows))
    require(totals == EXPECTED_TOTALS, totals)
    require(wall_totals == EXPECTED_WALL, wall_totals)
    require(order_totals == EXPECTED_ORDER, order_totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(wall_farkas == EXPECTED_WALL_FARKAS, wall_farkas)
    require(order_farkas == EXPECTED_ORDER_FARKAS, order_farkas)
    canonical_screen = tuple(row[:19] for row in screened)
    screen_sha = hashlib.sha256(repr(canonical_screen).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)

    residual_with_index = tuple((index, row) for index, row in enumerate(screened) if row[12])
    require(all(row[5] for _index, row in residual_with_index), "order residual")
    residual = tuple((row[1], row[13]) for _index, row in residual_with_index)
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(residual_sha == EXPECTED_RESIDUAL_SHA256, residual_sha)
    require(len(residual_with_index) == len(RESIDUAL), len(residual_with_index))
    quotient_counts = Counter()
    for (atlas_index, row), expected in zip(residual_with_index, RESIDUAL):
        expected_index, expected_body, expected_count, expected_mask_sha = expected
        require((atlas_index, row[1], row[12]) == (expected_index, expected_body, expected_count), expected)
        require(hashlib.sha256(repr(row[13]).encode()).hexdigest() == expected_mask_sha, expected_body)
        for ds in row[13]:
            divisor = lcm(*ds)
            require(divisor % row[4] == 0 and row[2] % divisor == 0, (row[1], ds))
            quotient_counts[(gcd(LEVEL, row[2]), divisor // row[4])] += 1
    quotient_packet = tuple((g, q, count) for (g, q), count in sorted(quotient_counts.items()))
    require(quotient_packet == EXPECTED_QUOTIENT_COUNTS, quotient_packet)

    terminal_tasks = tuple(
        (local, atlas_index, row[1], row[13])
        for local, (atlas_index, row) in enumerate(residual_with_index)
    )
    if args.processes == 1:
        terminal_pairs = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(terminal_tasks))) as pool:
            terminal_pairs = tuple(pool.imap_unordered(terminal_worker, terminal_tasks, chunksize=1))
    terminal_by_index = {local: (atlas_index, row) for local, atlas_index, row in terminal_pairs}
    require(len(terminal_by_index) == len(terminal_tasks), "lost terminal task")
    terminals = tuple(terminal_by_index[index][1] for index in range(len(terminal_tasks)))
    terminal_summary = (
        len(terminals), sum(row[4] for row in terminals), sum(row[7] for row in terminals),
        sum(row[8] for row in terminals), sum(row[9] for row in terminals),
        sum(row[10] for row in terminals), sum(row[11] for row in terminals),
        sum(row[12] for row in terminals), sum(row[14] for row in terminals),
        min(row[16] for row in terminals), sum(row[19] for row in terminals),
        sum(row[20] for row in terminals),
    )
    require(terminal_summary == EXPECTED_TERMINAL, terminal_summary)
    terminal_semantic_sha = hashlib.sha256(
        repr(tuple(terminal_semantic(row) for row in terminals)).encode()
    ).hexdigest()
    require(terminal_semantic_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256, terminal_semantic_sha)
    closure_packet = (
        "z218-terminal-scout-v1",
        tuple(closure_row(row) for row in terminals),
    )
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_sha == EXPECTED_CASE_VECTOR_SHA256, case_sha)
    min_terminal = min(terminals, key=lambda row: row[5])
    min_gap = (min_terminal[1], min_terminal[5].numerator, min_terminal[5].denominator)
    require(min_gap == EXPECTED_MIN_GAP, min_gap)

    require(LEDGER_BEFORE - LAYER_ROWS == INTERMEDIATE_LEDGER, "z218 ledger arithmetic")
    require(INTERMEDIATE_LEDGER - COMPOSED_ROWS == LEDGER_AFTER, "z217 composition arithmetic")

    lines = [
        "LRC14 THM3251 projected-k3 z218 terminal descent and composed cap216",
        (
            f"engine=THM3139_source:{SOURCE_3139_SHA256};output:{OUTPUT_3139_SHA256};"
            f"semantic:{SEMANTIC_3139_SHA256}"
        ),
        (
            f"dependencies=THM3230_source:{SOURCE_3230_SHA256};output:{OUTPUT_3230_SHA256};"
            f"semantic:{SEMANTIC_3230_SHA256};THM3242_source:{SOURCE_3242_SHA256};"
            f"output:{OUTPUT_3242_SHA256};semantic:{SEMANTIC_3242_SHA256}"
        ),
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        (
            f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};"
            f"row_sha256:{EXPECTED_ROW_SHA256};states:{totals[0]};crude:{totals[1]};"
            f"status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};"
            f"legacy_farkas:{farkas[1]};screen_sha256:{screen_sha};"
            f"residual_rows:{len(residual_with_index)};residual_sha256:{residual_sha};"
            f"residual_index_sha256:{residual_index_sha}"
        ),
        (
            f"wall_screen=states:{wall_totals[0]};crude:{wall_totals[1]};"
            f"status:{wall_totals[2]};residual:{wall_totals[3]};"
            f"direct_farkas:{wall_farkas[0]};legacy_farkas:{wall_farkas[1]}"
        ),
        (
            f"order_screen=states:{order_totals[0]};crude:{order_totals[1]};"
            f"status:{order_totals[2]};residual:{order_totals[3]};"
            f"direct_farkas:{order_farkas[0]};legacy_farkas:{order_farkas[1]}"
        ),
        (
            f"modulus=strata:{strata};first_d_equals_L_over_2:1;"
            f"residual_quotients:{quotient_packet}"
        ),
        (
            f"terminal=rows:{terminal_summary[0]};masks:{terminal_summary[1]};"
            f"positive_two_high_gap:{terminal_summary[10]};closed:{terminal_summary[11]};"
            f"zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};"
            f"low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};"
            f"exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};"
            f"failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};"
            f"minimum_gap:{EXPECTED_MIN_GAP[1]}/{EXPECTED_MIN_GAP[2]};"
            f"minimum_gap_body:{EXPECTED_MIN_GAP[0]};"
            f"semantic_sha256:{terminal_semantic_sha};closure_sha256:{closure_sha};"
            f"case_sha256:{case_sha}"
        ),
    ]
    for local, ((atlas_index, _screen_row), row) in enumerate(zip(residual_with_index, terminals)):
        lines.append(
            "TERMINAL;"
            f"local={local};atlas={atlas_index};E={row[1]};L={row[2]};high={row[3]};"
            f"masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};"
            f"one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};"
            f"exact={row[11]};maxgap={row[12]};failed={row[14]};"
            f"minimum_slack={row[16]};case_sha256={row[17]}"
        )
    lines.extend(
        (
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;both_order_rows_have_zero_residual",
            "direction_modulus=every_row_has_gcd(218,L)=2_and_first_d=L/2;every_residual_D_over_first_d_equals_2;this_is_pointed_divisor_metadata_not_a_C2_action_or_carrier_transition",
            "direction_terminal=every_residual_is_a_wall_row;the_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_complete_one_high_bank_enlarges_the_actual_assignment_set",
            "direction_composition=only_after_z218_closes_does_promoted_THM3242_remove_the_disjoint_z217_layer;the_two_gcd_strata_create_no_common_C2_or_C7_carrier",
            "evidence_boundary=the_companion_freshly_replays_all_screen_and_terminal_workers;the_prior_interleaved_checkpoint_is_not_a_dependency_or_ordering_witness;all_24_mask_banks_are_bound_individually",
            f"intermediate_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={INTERMEDIATE_LEDGER};projected_k3_cap:z1<={INTERMEDIATE_CAP}",
            f"promotion_consequence=ledger {INTERMEDIATE_LEDGER}-{COMPOSED_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_sha256:{NEXT_ROW_SHA256};status:occupied",
            "scope=projected_k3_necessary_atlas_only;no_C2_action_or_modular_free_factor_or_physical_cover_classification;no_k_le_1_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
