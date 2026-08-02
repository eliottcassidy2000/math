#!/usr/bin/env python3
"""Exact projected-k3 z1=227 inherited-screen descent.

The complete z1=227 atlas layer is rebuilt from the pinned THM-2941 atlas and
passed through the independently promoted THM-3113/THM-3111 screen.  This
script deliberately stops at the screen: the residual bank is empty.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
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
    "lrc14_j7_k3_z227_exact_screen_zero_residual_descent_thm3114.out"
)

SOURCE_3113_SHA256 = "59922de78673b94711d0cfe45c9954b8b2056b85f06e904d690ef18f389d111e"
OUTPUT_3113_SHA256 = "a872fe57e12f878b8e237ba71611a39c1f03ba893d45d3c3ae4dd6480e660595"
SEMANTIC_3113_SHA256 = "a024972dce2c159096db81a751b49654ebc8fad18e2fb1229dbfa0f01784e7e4"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"

LEVEL = 227
NEXT_LEVEL = 226
EXPECTED_CENSUS = (30, 28, 2)
EXPECTED_ROW_ORDER_SHA256 = "17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c"
EXPECTED_SCREEN = (284, 143, 141, 0)
EXPECTED_ORDER = (9, 9, 0, 0)
EXPECTED_FARKAS = (0, 141)
EXPECTED_SCREEN_SHA256 = "a86dc48861dbbf23d86cea24201a1b9634c118eec6e124effb5fa4e6aea688df"
EXPECTED_NEXT_CENSUS = (13, 13, 0)
EXPECTED_NEXT_ROW_ORDER_SHA256 = "ec2f1218d61dd69d3811d68eb87a23254d276bb75647be2d4c883affe53a520e"

LEDGER_BEFORE = 374215
LAYER_ROWS = 30
LEDGER_AFTER = 374185
NEXT_CAP = 226
EXPECTED_SEMANTIC_SHA256 = "17347efff69b101d47afc7d7097fcde498d210ed25737c1ef9935d6b3a5baacb"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    payload = path.read_bytes()
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_layers(base) -> tuple[tuple, tuple]:
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    selected = {LEVEL: [], NEXT_LEVEL: []}
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
        selected[level].append(
            (body, int(match.group(2)), high, level < high)
        )
    rows = tuple(selected[LEVEL])
    next_rows = tuple(selected[NEXT_LEVEL])
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == EXPECTED_NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256,
            "row-order digest")
    require(hashlib.sha256(repr(next_rows).encode()).hexdigest()
            == EXPECTED_NEXT_ROW_ORDER_SHA256, "next-row-order digest")
    return rows, next_rows


def screen_worker(task):
    predecessor = load(f"thm3114_predecessor_{task[1]}", SOURCE_3113)
    return predecessor.screen_worker(task)


def run_screen(rows: tuple, processes: int) -> tuple:
    tasks = tuple(
        (LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows
    )
    if processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost or duplicated task")
    screened = tuple(by_task[task] for task in tasks)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    order = tuple(row for row in screened if not row[5])
    order_totals = tuple(sum(row[index] for row in order) for index in (9, 10, 11, 12))
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(totals == EXPECTED_SCREEN, totals)
    require(order_totals == EXPECTED_ORDER, order_totals)
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "status verification")
    require(not any(row[12] for row in screened), "nonempty residual bank")
    screen_sha = hashlib.sha256(repr(screened).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    return screened, totals, order_totals, farkas, screen_sha


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3113) == SOURCE_3113_SHA256, "THM-3113 source changed")
    require(sha(OUTPUT_3113) == OUTPUT_3113_SHA256, "THM-3113 output changed")
    require(f"semantic_sha256={SEMANTIC_3113_SHA256}" in OUTPUT_3113.read_text(),
            "THM-3113 semantic changed")

    prior = load("thm3114_prior", SOURCE_3113)
    driver = prior.load("thm3114_driver")
    predecessor = driver.load("thm3114_screen")
    base = predecessor.load("thm3114_base")
    require(base.thm.ATLAS_SHA256 == ATLAS_SHA256, "atlas semantic changed")
    rows, next_rows = parse_layers(base)
    screened, totals, order_totals, farkas, screen_sha = run_screen(rows, args.processes)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    semantic_packet = (
        "lrc14-k3-z227-zero-residual-screen-v1",
        SOURCE_3113_SHA256,
        OUTPUT_3113_SHA256,
        SEMANTIC_3113_SHA256,
        ATLAS_SHA256,
        rows,
        screened,
        next_rows,
        (LEDGER_BEFORE, LAYER_ROWS, LEDGER_AFTER, NEXT_CAP),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected k3 z227 exact-screen zero-residual descent",
        f"dependency=THM3113_source:{SOURCE_3113_SHA256};output:{OUTPUT_3113_SHA256};semantic:{SEMANTIC_3113_SHA256}",
        f"atlas=sha256:{ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{EXPECTED_CENSUS[0]};wall:{EXPECTED_CENSUS[1]};order:{EXPECTED_CENSUS[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};order_control:{order_totals}",
    ]
    for row in screened:
        lines.append(
            f"ROW;z1={LEVEL};E={row[1]};L={row[2]};high={row[3]};"
            f"branch={'wall' if row[5] else 'order'};states={row[9]};"
            f"crude={row[10]};status={row[11]};residual={row[12]};"
            f"direct={row[19]};legacy={row[20]};stage_sha256={row[18]};"
            f"residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()}"
        )
    lines.extend([
        "direction=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;zero_residual_therefore_closes_the_full_layer",
        f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
        f"next_layer=z1:{NEXT_LEVEL};rows:{EXPECTED_NEXT_CENSUS[0]};wall:{EXPECTED_NEXT_CENSUS[1]};order:{EXPECTED_NEXT_CENSUS[2]};row_order_sha256:{EXPECTED_NEXT_ROW_ORDER_SHA256};status:occupied_unscouted_handoff",
        "scope=projected_k3_necessary_atlas_only;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ])
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
