#!/usr/bin/env python3
"""Exact projected-k3 z217 status annihilation."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from math import gcd, lcm
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
    "lrc14_j7_k3_z217_exact_status_annihilation_thm3242.out"
)

SOURCE_3109_SHA256 = "1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c"
OUTPUT_3109_SHA256 = "e6d56fa6a419ffe229b8334090e02d98c9d2cdf5f2fa5e24baedd5f722dadf70"
SEMANTIC_3109_SHA256 = "5be5c2fc680d6873600e77227f51264d74a7cd353652795e5ef74215e0fda843"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 217
EXPECTED_CENSUS = (8, 7, 1)
EXPECTED_ROW_SHA256 = "885fa20869364ee4a808324612ad373278464a339bf47e183418c1731af4ea6d"
EXPECTED_TOTALS = (66, 0, 66, 0)
EXPECTED_WALL = (64, 0, 64, 0)
EXPECTED_ORDER = (2, 0, 2, 0)
EXPECTED_FARKAS = (0, 66)
EXPECTED_WALL_FARKAS = (0, 64)
EXPECTED_ORDER_FARKAS = (0, 2)
EXPECTED_SCREEN_SHA256 = "69e3333a4b1dfbbcdd8530d933aa787b7a9b295f01308fcb613f3e65e8eed11b"
EMPTY_SHA256 = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
EXPECTED_SEMANTIC_SHA256 = "a28872885af6aef23f3aa00fe026b5f4c23b6d48302924fab37279cc2a4daf50"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name):
    spec = importlib.util.spec_from_file_location(name, SOURCE_3109)
    require(spec is not None and spec.loader is not None, SOURCE_3109)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_rows(base):
    pattern = re.compile(
        r"^row=E=([0-9,]+);h=[^;]+;r=([0-9]+);"
        r"L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    rows = []
    for line in base.thm.ATLAS.read_text(encoding="utf-8").splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(5)) != LEVEL:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        ruler = int(match.group(3))
        high = int(match.group(4))
        rows.append((body, ruler, high, LEVEL < high))
    return tuple(rows)


def screen_worker(indexed_task):
    index, task = indexed_task
    source = load(f"thm3242_worker_{index}")
    returned_task, row = source.isolated_screen_worker(task)
    require(returned_task == task, (index, returned_task, task))
    return index, task, tuple(row[:19]), row[19], row[20]


def summarize(records, rows):
    require(tuple(sorted(records)) == tuple(range(len(rows))), "index coverage")
    canonical = tuple(records[index][2] for index in range(len(rows)))
    require(
        all(records[index][1] == (LEVEL, *rows[index]) for index in range(len(rows))),
        "task order",
    )
    for row in canonical:
        require(row[0] == LEVEL, row[0])
        require(row[2] == 14 * lcm(*row[1]), (row[1], row[2]))
        require(gcd(LEVEL, row[2]) == 7, (row[1], row[2]))
        require(row[4] == row[2] // 7, (row[1], row[4], row[2] // 7))
        require(row[16] == row[11], (row[1], "status verification"))
        require(row[12] == 0 and row[13] == (), (row[1], "residual"))

    totals = tuple(sum(row[index] for row in canonical) for index in (9, 10, 11, 12))
    wall_rows = tuple(row for row in canonical if row[5])
    order_rows = tuple(row for row in canonical if not row[5])
    wall = tuple(sum(row[index] for row in wall_rows) for index in (9, 10, 11, 12))
    order = tuple(sum(row[index] for row in order_rows) for index in (9, 10, 11, 12))
    farkas = (
        sum(records[index][3] for index in records),
        sum(records[index][4] for index in records),
    )
    wall_farkas = (
        sum(records[index][3] for index in records if records[index][2][5]),
        sum(records[index][4] for index in records if records[index][2][5]),
    )
    order_farkas = (
        sum(records[index][3] for index in records if not records[index][2][5]),
        sum(records[index][4] for index in records if not records[index][2][5]),
    )
    digest = hashlib.sha256(repr(canonical).encode()).hexdigest()
    residual = tuple((row[1], row[13]) for row in canonical if row[12])
    residual_digest = hashlib.sha256(repr(residual).encode()).hexdigest()

    require(totals == EXPECTED_TOTALS, totals)
    require(wall == EXPECTED_WALL, wall)
    require(order == EXPECTED_ORDER, order)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(wall_farkas == EXPECTED_WALL_FARKAS, wall_farkas)
    require(order_farkas == EXPECTED_ORDER_FARKAS, order_farkas)
    require(digest == EXPECTED_SCREEN_SHA256, digest)
    require(residual == () and residual_digest == EMPTY_SHA256, residual_digest)
    return canonical, totals, wall, order, farkas, wall_farkas, order_farkas, digest


def build_lines(summary, rows, census, row_sha, semantic, records):
    canonical, totals, wall, order, farkas, wall_farkas, order_farkas, digest = summary
    lines = [
        "LRC14 THM3242 projected-k3 z217 exact-status annihilation",
        (
            f"dependency=THM3109_source:{SOURCE_3109_SHA256};"
            f"output:{OUTPUT_3109_SHA256};semantic:{SEMANTIC_3109_SHA256};"
            f"atlas:{ATLAS_SHA256}"
        ),
        (
            f"layer=z1:{LEVEL};census:{census};row_sha256:{row_sha};"
            "all_gcd_seven:1;first_d_equals_L_over_7:1"
        ),
        (
            f"screen=total:{totals};wall:{wall};order:{order};farkas:{farkas};"
            f"wall_farkas:{wall_farkas};order_farkas:{order_farkas};"
            f"screen_sha256:{digest}"
        ),
        f"residual=bodies:0;masks:0;sha256:{EMPTY_SHA256}",
    ]
    for index, row in enumerate(canonical):
        lines.append(
            "ROW;"
            f"index={index};E={','.join(map(str, row[1]))};wall={int(row[5])};"
            f"L={row[2]};high={row[3]};first_d={row[4]};"
            f"states={row[9]};crude={row[10]};status={row[11]};"
            f"residual={row[12]};direct_farkas={records[index][3]};"
            f"legacy_farkas={records[index][4]}"
        )
    lines.extend(
        (
            "direction=all_66_states_are_excluded_by_raw_verified_legacy_Farkas_certificates;no_terminal_step_is_needed",
            "consequence=z217_is_empty_but_the_current_cap_remains218_and_ledger373411_until_z218_is_independently_discharged",
            "scope=z217_layer_only;gcd7_is_a_pointed_divisor_condition_not_a_C7_action_or_common_carrier;no_k_le_1_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(lf_sha(SOURCE_3109) == SOURCE_3109_SHA256, "THM3109 source changed")
    require(lf_sha(OUTPUT_3109) == OUTPUT_3109_SHA256, "THM3109 output changed")

    source = load("thm3242_source")
    base = source.load("thm3242_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, base.thm.ATLAS_SHA256)
    rows = atlas_rows(base)
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    row_sha = hashlib.sha256(repr(rows).encode()).hexdigest()
    require(census == EXPECTED_CENSUS, census)
    require(row_sha == EXPECTED_ROW_SHA256, row_sha)

    tasks = tuple((index, (LEVEL, *row)) for index, row in enumerate(rows))
    if args.processes == 1:
        completed = tuple(map(screen_worker, tasks))
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            completed = tuple(pool.imap_unordered(screen_worker, tasks, chunksize=1))
    records = {record[0]: record for record in completed}
    require(len(records) == len(completed), "duplicate worker index")
    summary = summarize(records, rows)

    packet = (
        "z217-exact-layer-referee-v1",
        SOURCE_3109_SHA256,
        ATLAS_SHA256,
        EXPECTED_CENSUS,
        EXPECTED_ROW_SHA256,
        EXPECTED_TOTALS,
        EXPECTED_WALL,
        EXPECTED_ORDER,
        EXPECTED_FARKAS,
        EXPECTED_SCREEN_SHA256,
        EMPTY_SHA256,
    )
    semantic = hashlib.sha256(repr(packet).encode()).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)
    lines = build_lines(summary, rows, census, row_sha, semantic, records)
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.output is not None:
        args.output.write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    mp.freeze_support()
    main()
