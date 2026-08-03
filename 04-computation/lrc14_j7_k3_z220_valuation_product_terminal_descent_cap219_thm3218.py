#!/usr/bin/env python3
"""Exact projected-k3 z220 valuation-product screen and terminal descent."""

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
SOURCE_3207 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.py"
)
OUTPUT_3207 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.out"
)

SOURCE_3207_SHA256 = "b3d1f0c451087017c1363d42be8789df78d5eec7db7a05b49dc5ca9e194f2091"
OUTPUT_3207_SHA256 = "1e4b01f6e5bfb179ad5fe6ba786124ff2c9fdf3833eed4917c0b3ce0abb7b76d"
SEMANTIC_3207_SHA256 = "aad27ae040933935c6eacfd2801cef1644df3de3a8f031aebf7d1b238ac35e74"

LEVEL = 220
EXPECTED_CENSUS = (289, 249, 40)
EXPECTED_ROW_ORDER_SHA256 = "bc75fef8aebdc6b350957ea96842b38a24bad2174200646ebaf2817844515799"
EXPECTED_G_STRATA = (
    (4, 72, 52, 20),
    (20, 183, 163, 20),
    (44, 5, 5, 0),
    (220, 29, 29, 0),
)
EXPECTED_SCREEN = (96780, 47054, 49356, 370)
EXPECTED_WALL_SCREEN = (95673, 46043, 49260, 370)
EXPECTED_ORDER_SCREEN = (1107, 1011, 96, 0)
EXPECTED_FARKAS = (12, 49344)
EXPECTED_WALL_FARKAS = (12, 49248)
EXPECTED_ORDER_FARKAS = (0, 96)
EXPECTED_SCREEN_SHA256 = "f7168e67f3eb21529edcb73cae97815804827cfdb0b2fe7df66dbd37a5e8ba07"
EXPECTED_RESIDUAL_SHA256 = "71525cd0a7e90047f6e391e769142b57c37f09506d656a52d0dbde57da1502f7"
EXPECTED_RESIDUAL_INDEX_SHA256 = "db97f96ac8dc62f61bbc912068e1e1d2f525efad732d8fb88e4cefa57389d259"
EXPECTED_QUOTIENT_COUNTS = (
    (4, 4, 2),
    (20, 10, 15),
    (20, 20, 295),
    (44, 44, 5),
    (220, 220, 53),
)
EXPECTED_TERMINAL = (31, 370, 326, 644, 163, 589, 55, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "4254a3d491e57cbdb3e1a9fda66a002e81bbf3989d39f877e5d2b86357073f11"
EXPECTED_CLOSURE_SHA256 = "afc495948dcc329795cad66cd0aee6a8a616971a99acd67e3b3a477e353d7839"
EXPECTED_CASE_VECTOR_SHA256 = "41500dcc3c8273b9ec8e0d545ca1d6a9f8a52cfb03f15088c143d2d7473e25cc"
EXPECTED_MIN_GAP = ((1, 4, 9, 10, 12, 14), 28369051, 33158840715)
EXPECTED_SEMANTIC_SHA256 = "ed67ff55840d07965a6cf8a478c0fc49cab637511ad3d684b2de8831dfcfcd94"

LEDGER_BEFORE = 373716
LAYER_ROWS = 289
LEDGER_AFTER = 373427
NEXT_CAP = 219
NEXT_LEVEL = 219
NEXT_CENSUS = (16, 15, 1)
NEXT_ROW_ORDER_SHA256 = "6a39c83d47b8080fd52e64479b886f1b6a887150f32cfbf30a676fbb0aaa54ba"

RESIDUAL = {
    (1, 2, 9, 10, 12, 14): (2, "4f32da16d53bcf6e1a28ccc785dc2f42930f948fde68a7dec2ff26b34acc5ff2"),
    (1, 3, 4, 7, 9, 12): (2, "26fc9dfbeee2d04254823f79f0a9603b1a0ea69c4ec833a2e39c2770bfaec21c"),
    (1, 3, 8, 10, 12, 14): (2, "c5aaa76dc3ca5eb08e385bb77263f203626f4f76b58477c29581390bae4166ed"),
    (1, 4, 5, 9, 12, 14): (2, "34d2977ab7eb9178fccda734367927d5e57d4aba63b77163e7e78af2c3e9940e"),
    (1, 4, 7, 9, 10, 12): (3, "05125c59e1eb6dd5b01bc0609f2a199281c0b9fb9285a3d36bc8f86bb2043c45"),
    (1, 4, 9, 10, 12, 14): (50, "c09bbbb524f65617f30bdc220c8caecb0c25637d5c848f7e9ba49377ee0d848e"),
    (1, 4, 10, 11, 12, 14): (2, "78e9a4573080110874d551f5965a7d19365a2ea142a4c2607c60ad8bcd8ee1ca"),
    (1, 4, 10, 12, 13, 14): (36, "89e3286355495a6b58cf97dcd849124efe3829bdb97b4ded7b30f3c3c9114087"),
    (1, 5, 8, 9, 12, 14): (2, "d3c369865b6310b8dd8679315bacef642347f1a7630bf7d28a14b3acd640c0f7"),
    (1, 6, 7, 9, 10, 12): (5, "bb541cfc99af507b0fd6ae9b5050ebd80a839ef8cf51fa1cc603b5bf88e4f009"),
    (1, 6, 8, 9, 10, 14): (3, "a97f4ac31802dcc61a98b9259cfa88510edf77cd081d6e538f0672aa7f61980e"),
    (1, 6, 8, 10, 13, 14): (8, "6af5b1c38e1558ebb8c0d1d64e2959e7a99ce96164d394ca0ddb025c866abb35"),
    (1, 6, 9, 10, 12, 14): (2, "dea84d9609e712661b3fc08bb406a27fccd782ca316662251b7962bfe3a2e1ff"),
    (1, 6, 10, 12, 13, 14): (2, "a54886566611f49975a1ea9a44fb1d2e36a6a1c5d3343ca84af7a7c420ed9605"),
    (1, 7, 8, 10, 12, 13): (6, "706ae5baaeacd5e653681ea73c377338a838872950ade1055280816e92fbb2b4"),
    (1, 8, 9, 10, 12, 13): (2, "9619e38cb9ed84c94ce376a15c8d56ec2470e4a828568db8468aefb8fae34c17"),
    (1, 8, 9, 10, 12, 14): (15, "335ca240eae9f0059bb4f5ca01843860407133e474d4eb2da8cc4cd2645bc3c7"),
    (1, 8, 10, 11, 12, 14): (22, "6cc492ca868235bd81c0af043dcd5c606e80674559d39afd99d1278831593f8b"),
    (1, 8, 10, 12, 13, 14): (142, "bc9c05b856056eeab2c8f900c84492785c17d16c4b6b5b6ab81c1d0ed45b87d6"),
    (2, 3, 8, 10, 12, 14): (2, "53930e30ba25041a0e0a16ba9a10ecaa35c3189d6c100c453b18052b564725b6"),
    (2, 5, 8, 9, 12, 14): (3, "47605b11f3487c8a3094fe93280a254eafc567eac1bfbbb43b82a84c2c4cc62a"),
    (2, 5, 9, 11, 12, 14): (12, "321380fdfd2869d0cdfb9feed3ac745a5b89838deb2e57dafb44343a7da4276a"),
    (2, 6, 8, 10, 12, 14): (15, "ddd9ccb50ddd00659925b3c0e0dce223ab3e778d47d7bd0b7c4d27fda9738969"),
    (2, 6, 9, 10, 12, 14): (1, "93b75b16d47aefd555383fb728281698d77e437b7b57aba3be2a1d1fa0cc37b6"),
    (2, 8, 9, 10, 12, 14): (3, "0c0b74adfd998b31987f8894837600e5af9fa4e38be6ebd59c66f4439308e281"),
    (2, 8, 9, 11, 12, 14): (5, "06e2cbd2101f0e91e001d8e0f8c4fefac8919f5034c3b22ea925d4630380b9d1"),
    (2, 8, 10, 11, 12, 14): (14, "e04a6174354a2e58bf5a24e3e60615b60e10ed83a611e0ae3fd7414a86bda08c"),
    (2, 8, 10, 12, 13, 14): (1, "bb13f55da3a6960a55c06c6b126a32c57c284063edf32369379f1a6e2a2724b9"),
    (2, 9, 10, 11, 12, 14): (3, "d2e93d5a182b82261ab3c62dc13c3fd5f12f4f29d2b236902dbf2bef64ed95c3"),
    (3, 4, 9, 10, 12, 14): (2, "3dec485c5a4be6168762ef0daa39d1ed7337b2113b1775c25a862b10bb450a81"),
    (3, 8, 9, 10, 12, 14): (1, "caab6a60429a426e91c9e990eaa57a1b51eec2149c158aa842440b759052fbc0"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3207):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def screen_worker(task):
    return load(f"thm3218_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    index, body, bank = task
    canonical = load(f"thm3218_terminal_canonical_{index}")
    wrapper = canonical.load(f"thm3218_terminal_wrapper_{index}")
    entry = wrapper.load(f"thm3218_terminal_entry_{index}")
    source = entry.load(f"thm3218_terminal_source_{index}")
    driver = source.load(f"thm3218_terminal_driver_{index}")
    predecessor = driver.load(f"thm3218_terminal_predecessor_{index}")
    bridge = predecessor.load(f"thm3218_terminal_bridge_{index}")
    base = bridge.load(f"thm3218_terminal_base_{index}")
    return index, base.thm.terminal_probe((LEVEL, body, bank))


def modulus_data(task, row=None):
    level, body, ruler, high, wall = task
    require(level == LEVEL, level)
    require(ruler == 14 * lcm(*body), (body, ruler, "body ruler"))
    g = gcd(level, ruler)
    require(g in (4, 20, 44, 220), (body, ruler, g))
    first_d = ruler // g
    if row is not None:
        require(row[:6] == (level, body, ruler, high, first_d, wall), (body, "header"))
        for ds in row[13]:
            d = lcm(*ds)
            require(d % first_d == 0, (body, ds, first_d))
            q = d // first_d
            require(q == g or (g, q) == (20, 10), (body, g, q))
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
        for g in (4, 20, 44, 220)
    )
    require(strata == EXPECTED_G_STRATA, strata)
    packet = (
        "lrc14-k3-z220-valuation-product-terminal-cap219-v1",
        (SOURCE_3207_SHA256, OUTPUT_3207_SHA256, SEMANTIC_3207_SHA256, base.thm.ATLAS_SHA256),
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
    require(sha(SOURCE_3207) == SOURCE_3207_SHA256, "THM-3207 source changed")
    require(sha(OUTPUT_3207) == OUTPUT_3207_SHA256, "THM-3207 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3207_SHA256}" in OUTPUT_3207.read_text(encoding="utf-8"),
        "THM-3207 semantic changed",
    )
    residual_index_sha = hashlib.sha256(repr(tuple(RESIDUAL.items())).encode()).hexdigest()
    require(residual_index_sha == EXPECTED_RESIDUAL_INDEX_SHA256, residual_index_sha)

    canonical = load("thm3218_canonical")
    wrapper = canonical.load("thm3218_wrapper")
    entry = wrapper.load("thm3218_entry")
    source = entry.load("thm3218_source")
    driver = source.load("thm3218_driver")
    predecessor = driver.load("thm3218_predecessor")
    bridge = predecessor.load("thm3218_bridge")
    base = bridge.load("thm3218_base")
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
        "z220-terminal-closure-solver-witness-free-v1",
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
        "LRC14 projected k3 z220 valuation-product terminal descent and cap219",
        f"dependency=THM3207_source:{SOURCE_3207_SHA256};output:{OUTPUT_3207_SHA256};semantic:{SEMANTIC_3207_SHA256}",
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
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;all_forty_order_rows_have_zero_residual",
            "direction_modulus=first_d=L/g_for_g=gcd(220,L);every_residual_D_over_first_d_is_the_top_divisor_g_except_fifteen_g20_masks_at_q10;the_occupied_values_10_and_20_are_vertices_of_the_ambient_top_two_cover_but_no_common_carrier_transition_is_proved",
            "direction_terminal=each_residual_is_a_wall_row;the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_the_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;all_31_residual_mask_banks_are_bound_individually",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied",
            "scope=projected_k3_necessary_atlas_only;the_divisor_product_is_a_pointed_join_semilattice_not_a_V4_torsor_or_PSL2_action;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
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
